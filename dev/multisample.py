"""Checks of the multi-sample significance search.

    python3 dev/multisample.py fwer 60 [min_support] [permutations]
                                              structureless samples: how often is anything reported?
    python3 dev/multisample.py replication    niches planted in all, one, or none of three samples

Each dataset is several independently generated tissues passed together as
repeated --input files. The family-wise guarantee is carried by the pooled
adjusted p-value; `replicated_in` counts samples whose own p-value reaches
alpha and is descriptive.
"""
import csv, os, random, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

import make_datasets as fixtures

DEV = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(DEV, "pam_st")
WORK = os.path.join(DEV, "multisample_out")
SAMPLES = 3


def structureless(seed):
    rng = random.Random(seed)
    return [[rng.uniform(0, fixtures.W), rng.uniform(0, fixtures.W),
             f"T{fixtures.pick(rng, fixtures.BASE)}"] for _ in range(fixtures.N)]


def write(rows, directory, name):
    path = os.path.join(directory, f"{name}.csv")
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["X_centroid", "Y_centroid", "Cell_Type"])
        writer.writerows(rows)
    return path


def run(inputs, out, permutations, min_support, seed, threads=1):
    args = [BINARY]
    for path in inputs:
        args += ["--input", path]
    args += ["--radius", "100", "--rho", "0.05", "--metric", "l2", "--statistic", "minp",
             "--min-support", str(min_support), "--split-size", "1000", "--max-motifs", "5",
             "--permutations", str(permutations), "--seed", str(seed), "--threads", str(threads),
             "--null-model", "block", "--block-size", "750", "--output-dir", out]
    subprocess.run(args, capture_output=True, text=True, check=True)
    return list(csv.DictReader(open(os.path.join(out, "motif_significance.csv"))))


def fwer(reps, permutations=2999, min_support=2):
    def one(rep):
        directory = tempfile.mkdtemp(dir=WORK)
        inputs = [write(structureless(10_000 * rep + s), directory, f"s{s}") for s in range(SAMPLES)]
        rows = run(inputs, os.path.join(directory, "out"), permutations, min_support, 1000 + rep)
        return (min((float(r["p_adjusted"]) for r in rows), default=1.0),
                min((float(r["p_raw"]) for r in rows), default=1.0),
                int(rows[0]["family_size"]) if rows else 0)
    with ThreadPoolExecutor(10) as pool:
        results = list(pool.map(one, range(reps)))
    families = [f for _, _, f in results if f > 0]
    print(f"{SAMPLES} structureless samples per dataset, {reps} datasets, B={permutations}; "
          f"{len(families)} formed a family (mean size {sum(families) / max(len(families), 1):.0f})")
    print(f"  pooled ADJUSTED p <= 0.05 in {sum(a <= 0.05 for a, _, _ in results)}/{reps}"
          f"  <- target about 5%")
    print(f"  pooled RAW      p <= 0.05 in {sum(p <= 0.05 for _, p, _ in results)}/{reps}")


def replication(permutations=2999):
    for planted_in in (3, 1, 0):
        directory = tempfile.mkdtemp(dir=WORK)
        inputs = []
        for s in range(SAMPLES):
            rows = fixtures.niches(seed=700 + s) if s < planted_in else structureless(900 + s)
            inputs.append(write(rows, directory, f"s{s}"))
        rows = run(inputs, os.path.join(directory, "out"), permutations, 5, 7, threads=4)
        planted = [r for r in rows if set(t.split(":")[0] for t in r["pattern"].split(";"))
                   == {"T3", "T5"}]
        significant = [r for r in rows if r["significant"] == "true"]
        top = planted[0] if planted else None
        summary = (f"top planted motif {top['pattern']}: pooled p_adj={top['p_adjusted']}, "
                   f"replicated_in={top['replicated_in']}/{SAMPLES}" if top else
                   "no T3+T5 motif reported")
        print(f"niches planted in {planted_in}/{SAMPLES} samples: "
              f"{len(significant)} significant motifs; {summary}")


if __name__ == "__main__":
    os.makedirs(WORK, exist_ok=True)
    if sys.argv[1] == "fwer":
        # Pooling samples grows the family, and the permutation budget must grow
        # with it (roughly 20x) or nothing can be rejected at all - which would
        # look like perfect calibration while measuring nothing.
        fwer(int(sys.argv[2]) if len(sys.argv) > 2 else 60,
             permutations=int(sys.argv[4]) if len(sys.argv) > 4 else 2999,
             min_support=int(sys.argv[3]) if len(sys.argv) > 3 else 2)
    else:
        replication()
