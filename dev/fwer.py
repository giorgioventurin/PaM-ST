"""False positives of the significance search where nothing should be found.

    python3 dev/fwer.py csr 40      uniform positions, i.i.d. labels: the null is exactly true
    python3 dev/fwer.py smooth 40   composition drifts across the tissue, still no local structure

At most about 5% of `csr` datasets should report anything at alpha = 0.05.
"""
import csv, os, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

import make_datasets as fixtures

DEV = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(DEV, "pam_st")
scenario, reps = sys.argv[1], int(sys.argv[2])


def dataset(rep):
    import math, random
    rng = random.Random(1000 + rep)
    rows = []
    for _ in range(fixtures.N):
        x, y = rng.uniform(0, fixtures.W), rng.uniform(0, fixtures.W)
        if scenario == "csr":
            label = fixtures.pick(rng, fixtures.BASE)
        else:
            t = 0.5 + 0.5 * math.sin(x / 800.0) * math.cos(y / 800.0)
            weights = [fixtures.BASE[i] * (t if i % 2 else 1 - t) + 1e-6 for i in range(fixtures.K)]
            total = sum(weights)
            label = fixtures.pick(rng, [w / total for w in weights])
        rows.append([x, y, f"T{label}"])
    return rows


def run(rep):
    work = os.path.join(DEV, "fwer_out")
    os.makedirs(work, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, dir=work) as f:
        writer = csv.writer(f)
        writer.writerow(["X_centroid", "Y_centroid", "Cell_Type"])
        writer.writerows(dataset(rep))
        path = f.name
    out = os.path.join(work, f"{scenario}_{rep}")
    try:
        done = subprocess.run([BINARY, "--input", path, "--radius", "100", "--rho", "0.05",
                               "--metric", "l2", "--statistic", "minp", "--min-support", "5",
                               "--split-size", "1000", "--max-motifs", "5", "--permutations", "199",
                               "--seed", str(1000 + rep), "--threads", "1", "--null-model", "block",
                               "--block-size", "750", "--output-dir", out],
                              capture_output=True, text=True)
        if done.returncode:
            return (1.0, 1.0, 0)
        rows = list(csv.DictReader(open(os.path.join(out, "motif_significance.csv"))))
        return (min((float(r["p_adjusted"]) for r in rows), default=1.0),
                min((float(r["p_raw"]) for r in rows), default=1.0),
                int(rows[0]["family_size"]) if rows else 0)
    finally:
        os.remove(path)


if __name__ == "__main__":
    with ThreadPoolExecutor(10) as pool:
        results = list(pool.map(run, range(reps)))
    families = [r for r in results if r[2] > 0]
    mean_family = sum(r[2] for r in families) / max(len(families), 1)
    print(f"{scenario}: {reps} datasets; {len(families)} produced a candidate family "
          f"(mean size {mean_family:.0f})")
    print(f"  ADJUSTED p <= 0.05 in {sum(1 for a, _, _ in results if a <= 0.05)}/{reps}"
          f"  <- target about 5%")
    print(f"  RAW      p <= 0.05 in {sum(1 for _, p, _ in results if p <= 0.05)}/{reps}")
