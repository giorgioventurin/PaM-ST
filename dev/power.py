"""Detection of planted two-type niches.

    python3 dev/power.py 12

Each dataset plants 40 niches of radius 150 filled with two rare types. The
search should recover that composition in essentially every dataset.
"""
import csv, os, subprocess, sys

import make_datasets as fixtures

DEV = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(DEV, "pam_st")
reps = int(sys.argv[1]) if len(sys.argv) > 1 else 12

if __name__ == "__main__":
    work = os.path.join(DEV, "power_out")
    os.makedirs(work, exist_ok=True)
    found = 0
    for rep in range(reps):
        path = os.path.join(work, f"niches_{rep}.csv")
        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["X_centroid", "Y_centroid", "Cell_Type"])
            writer.writerows(fixtures.niches(seed=500 + rep))
        out = os.path.join(work, str(rep))
        subprocess.run([BINARY, "--input", path, "--radius", "100", "--rho", "0.05", "--metric", "l2",
                        "--statistic", "minp", "--min-support", "5", "--split-size", "1000",
                        "--max-motifs", "5", "--permutations", "199", "--seed", str(500 + rep),
                        "--threads", "2", "--null-model", "block", "--block-size", "750",
                        "--output-dir", out], capture_output=True, text=True, check=True)
        reported = [r for r in csv.DictReader(open(os.path.join(out, "motif_significance.csv")))
                    if r["significant"] == "true"]
        if reported:
            found += 1
        if rep < 3 and reported:
            top = reported[0]
            print(f"  rep {rep}: {top['pattern']:28s} disjoint={top['disjoint_support']:>3} "
                  f"lift={float(top['lift']):8.1f} p_adj={top['p_adjusted']}")
    print(f"planted mixed niches: detected in {found}/{reps} datasets")
