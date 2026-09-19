"""Byte-for-byte output equivalence between two pam_st binaries.

    python3 dev/golden.py record BINARY   -> outputs of a known-good build into dev/golden/
    python3 dev/golden.py check  BINARY   -> outputs of a candidate build, diffed against it

Every output file is compared, along with stdout, stderr and exit codes. Timing
values are masked; nothing else may differ.
"""
import os, re, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

DEV = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(DEV)
FIXTURES = os.path.join(DEV, "data")
TESTS = os.path.join(ROOT, "tests", "data")
FULL = os.path.join(ROOT, "..", "data", "PDAC_202602", "annotated_cells",
                    "LSP31891_P181_rep1_annotated.csv")
DATA = {"niches": f"{FIXTURES}/niches.csv", "compart": f"{FIXTURES}/compartments.csv",
        "crop": f"{FIXTURES}/real_crop.csv", "sparse": f"{FIXTURES}/sparse.csv", "tiny": f"{TESTS}/block_null_cells.csv",
        "frozen": f"{TESTS}/frozen_cells.csv", "prop": f"{TESTS}/proportional_histograms.csv"}


def cases():
    c = {}
    base = ["--permutations", "7", "--max-motifs", "3", "--seed", "5", "--threads", "4"]
    for d in ["niches", "crop"]:
        for metric, rho in [("l2", "0"), ("l2", "0.05"), ("l2", "0.3"), ("l2", "1.5"),
                            ("js", "0"), ("js", "0.05"), ("js", "0.3"), ("js", "0.9")]:
            for mode in ["overlapping", "covering"]:
                for null in ["global", "block"]:
                    extra = (["--null-model", "block", "--block-size", "500",
                              "--block-origin-x", "37"] if null == "block" else [])
                    c[f"{d}_{metric}{rho}_{mode}_{null}"] = [
                        "--input", DATA[d], "--radius", "100", "--metric", metric,
                        "--rho", rho, "--neighborhood-mode", mode, *base, *extra]
        c[f"{d}_fewrs_l2"] = ["--input", DATA[d], "--radius", "100", "--rho", "0.05",
                              "--error-control", "fewrs-fdr", "--alpha", "0.2",
                              "--permutations", "auto", "--max-motifs", "12", "--threads", "3"]
        c[f"{d}_fewrs_js_block"] = ["--input", DATA[d], "--radius", "100", "--metric", "js",
                                    "--rho", "0.1", "--error-control", "fewrs-fdr", "--alpha", "0.2",
                                    "--fdr-failure-probability", "0.05", "--fdr-min-discoveries", "4",
                                    "--permutations", "80", "--max-motifs", "6",
                                    "--null-model", "block", "--block-size", "700"]
        for mode in ["overlapping", "covering"]:
            c[f"{d}_sweep_{mode}"] = ["--input", DATA[d], "--radius", "100", "--rho", "0.05",
                                      "--permutations", "9", "--max-motifs", "2",
                                      "--freeze-by-abundance", "--max-freeze-stages", "3",
                                      "--neighborhood-mode", mode, "--null-model", "block",
                                      "--block-size", "400"]
        c[f"{d}_sweep_none_global"] = ["--input", DATA[d], "--radius", "60", "--rho", "0",
                                       "--permutations", "5", "--freeze-by-abundance",
                                       "--stage-error-control", "none"]
        c[f"{d}_k1_threads1"] = ["--input", DATA[d], "--radius", "100", "--permutations", "19",
                                 "--threads", "1"]
        c[f"{d}_k1_threads10"] = ["--input", DATA[d], "--radius", "100", "--permutations", "19",
                                  "--threads", "10"]
        c[f"{d}_radius0"] = ["--input", DATA[d], "--radius", "0", "--permutations", "5",
                             "--max-motifs", "4"]
        c[f"{d}_radius300_js"] = ["--input", DATA[d], "--radius", "300", "--metric", "js",
                                  "--rho", "0.08", "--permutations", "4", "--max-motifs", "2",
                                  "--null-model", "block", "--block-size", "900"]
        # significance search: family selection, split, both nulls
        c[f"{d}_minp"] = ["--input", DATA[d], "--radius", "100", "--rho", "0.05", "--statistic",
                          "minp", "--min-support", "5", "--split-size", "1000", "--max-motifs", "5",
                          "--permutations", "49", "--seed", "5", "--threads", "4"]
        c[f"{d}_minp_block_types"] = ["--input", DATA[d], "--radius", "100", "--rho", "0.05",
                                      "--statistic", "minp", "--min-support", "5",
                                      "--split-size", "1000", "--min-types", "2",
                                      "--min-type-cells", "2", "--max-motifs", "5",
                                      "--permutations", "49", "--seed", "5", "--threads", "4",
                                      "--null-model", "block", "--block-size", "500"]
    # A family forms but nothing clears the reporting floor on the testing half,
    # which used to leave no rows at all and crash the writers.
    c["sparse_minp_nothing_reported"] = ["--input", DATA["sparse"], "--radius", "60", "--rho", "0.05",
                                         "--statistic", "minp", "--min-support", "2",
                                         "--split-size", "300", "--max-motifs", "4",
                                         "--permutations", "9", "--seed", "7", "--threads", "3",
                                         "--null-model", "block", "--block-size", "200"]
    c["sparse_minp_empty_family"] = ["--input", DATA["sparse"], "--radius", "60", "--rho", "0.05",
                                     "--statistic", "minp", "--min-support", "40",
                                     "--split-size", "300", "--max-motifs", "4",
                                     "--permutations", "9", "--seed", "7", "--threads", "2"]
    # Several samples in one run: pooled family, pooled test, per-sample breakdown.
    c["two_samples_minp_block"] = ["--input", DATA["niches"], "--input", DATA["compart"],
                                   "--radius", "100", "--rho", "0.05", "--statistic", "minp",
                                   "--min-support", "5", "--split-size", "1000", "--max-motifs", "5",
                                   "--permutations", "49", "--seed", "5", "--threads", "3",
                                   "--null-model", "block", "--block-size", "500"]
    c["three_samples_minp_global_js"] = ["--input", DATA["niches"], "--input", DATA["compart"],
                                         "--input", DATA["sparse"], "--radius", "80", "--metric", "js",
                                         "--rho", "0.1", "--statistic", "minp", "--min-support", "3",
                                         "--split-size", "700", "--max-motifs", "4",
                                         "--permutations", "29", "--seed", "9", "--threads", "2"]
    c["err_samples_without_minp"] = ["--input", DATA["tiny"], "--input", DATA["tiny"],
                                     "--permutations", "1", "--no-outdir"]
    c["niches_frozen2"] = ["--input", DATA["niches"], "--radius", "100", "--permutations", "9",
                           "--max-motifs", "3", "--freeze-cell-type", "T0", "--freeze-cell-type", "T1",
                           "--null-model", "block", "--block-size", "300"]
    c["crop_frozen_cov"] = ["--input", DATA["crop"], "--radius", "100", "--permutations", "9",
                            "--max-motifs", "2", "--freeze-cell-type", "Other_fibroblasts",
                            "--neighborhood-mode", "covering"]
    c["compart_k5_l2"] = ["--input", DATA["compart"], "--radius", "80", "--rho", "0.1",
                          "--permutations", "15", "--max-motifs", "5"]
    c["full_l2_global"] = ["--input", FULL, "--radius", "100", "--rho", "0.05",
                           "--permutations", "6", "--max-motifs", "2", "--threads", "10"]
    c["full_js_block"] = ["--input", FULL, "--radius", "100", "--metric", "js", "--rho", "0.05",
                          "--permutations", "4", "--null-model", "block", "--block-size", "500",
                          "--threads", "10"]
    for name, extra in [("tiny", []), ("tiny_block", ["--null-model", "block", "--block-size", "10"]),
                        ("tiny_single", ["--null-model", "block", "--block-size", "0.05"])]:
        c[name] = ["--input", DATA["tiny"], "--radius", "1", "--rho", "0", "--permutations", "19",
                   "--max-motifs", "3", *extra]
    c["frozen_sweep_all"] = ["--input", DATA["frozen"], "--radius", "1", "--rho", "0",
                             "--permutations", "5", "--max-motifs", "2", "--freeze-by-abundance"]
    c["prop_js"] = ["--input", DATA["prop"], "--radius", "1", "--metric", "js", "--rho", "0.05",
                    "--permutations", "5", "--max-motifs", "3"]
    c["report_file_only"] = ["--input", DATA["tiny"], "--radius", "1", "--permutations", "3",
                             "--output-file", "OUTDIR/custom/r.txt", "--no-outdir"]
    bad = {"missing_input": [], "unknown_arg": ["--input", DATA["tiny"], "--bogus"],
           "bad_metric": ["--input", DATA["tiny"], "--metric", "cos"],
           "bad_null": ["--input", DATA["tiny"], "--null-model", "grid"],
           "bad_mode": ["--input", DATA["tiny"], "--neighborhood-mode", "x"],
           "bad_ec": ["--input", DATA["tiny"], "--error-control", "bh"],
           "bad_statistic": ["--input", DATA["tiny"], "--statistic", "zscore"],
           "minp_sweep": ["--input", DATA["tiny"], "--statistic", "minp", "--freeze-by-abundance"],
           "family_without_minp": ["--input", DATA["tiny"], "--min-support", "5"],
           "split_too_small": ["--input", DATA["tiny"], "--statistic", "minp", "--radius", "100",
                               "--split-size", "50"],
           "bad_sec": ["--input", DATA["tiny"], "--freeze-by-abundance", "--stage-error-control", "x"],
           "block_wo_null": ["--input", DATA["tiny"], "--block-size", "5"],
           "block_nan": ["--input", DATA["tiny"], "--null-model", "block", "--block-size", "nan"],
           "alpha": ["--input", DATA["tiny"], "--alpha", "1"],
           "k0": ["--input", DATA["tiny"], "--max-motifs", "0"],
           "auto_pointwise": ["--input", DATA["tiny"], "--permutations", "auto"],
           "perm0": ["--input", DATA["tiny"], "--permutations", "0"],
           "fewrs_few": ["--input", DATA["tiny"], "--error-control", "fewrs-fdr", "--permutations", "3"],
           "sweep_and_manual": ["--input", DATA["tiny"], "--freeze-by-abundance",
                                "--freeze-cell-type", "A"],
           "stages_wo_sweep": ["--input", DATA["tiny"], "--max-freeze-stages", "2"],
           "sweep_fewrs": ["--input", DATA["tiny"], "--freeze-by-abundance", "--error-control",
                           "fewrs-fdr", "--permutations", "auto"],
           "unknown_frozen": ["--input", DATA["tiny"], "--freeze-cell-type", "Z"],
           "all_frozen": ["--input", DATA["frozen"], "--freeze-cell-type", "P",
                          "--freeze-cell-type", "A", "--freeze-cell-type", "B"],
           "dup_frozen": ["--input", DATA["tiny"], "--freeze-cell-type", "A",
                          "--freeze-cell-type", "A"],
           "missing_file": ["--input", "/nonexistent.csv"], "missing_value": ["--input"],
           "neg_radius": ["--input", DATA["tiny"], "--radius", "-1"],
           "neg_rho": ["--input", DATA["tiny"], "--rho", "-0.1"], "help": ["--help"]}
    for name, args in bad.items():
        c["err_" + name] = [*args, "--no-outdir"]
    return c


MASK_ROW = re.compile(r"^([a-z_]*seconds[=,]).*$", re.M)


def normalise(text, name):
    text = MASK_ROW.sub(r"\1<t>", text)
    if name.endswith("freeze_sweep_summary.csv"):
        lines = text.split("\n")
        column = lines[0].split(",").index("analysis_seconds")
        out = [lines[0]]
        for line in lines[1:]:
            parts = [p for p in re.findall(r'"[^"]*"|[^,]*', line) if p != ""] if line else []
            if len(parts) > column:
                parts[column] = "<t>"
            out.append(",".join(parts))
        text = "\n".join(out)
    return text


def run_case(binary, root, name, args):
    out = os.path.join(root, name)
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    args = [a.replace("OUTDIR", out) for a in args]
    if "--no-outdir" in args:
        args.remove("--no-outdir")
    else:
        args += ["--output-dir", os.path.join(out, "o")]
    result = subprocess.run([binary, *args], capture_output=True, text=True)
    open(os.path.join(out, "stdout"), "w").write(result.stdout.replace(binary, "PAM"))
    open(os.path.join(out, "stderr"), "w").write(result.stderr.replace(binary, "PAM"))
    open(os.path.join(out, "exit"), "w").write(str(result.returncode))


def collect(root):
    files = {}
    for directory, _, names in os.walk(root):
        for name in names:
            path = os.path.join(directory, name)
            files[os.path.relpath(path, root)] = normalise(open(path, errors="replace").read(), name)
    return files


if __name__ == "__main__":
    mode, binary = sys.argv[1], os.path.abspath(sys.argv[2])
    golden, root = os.path.join(DEV, "golden"), os.path.join(DEV, "golden" if mode == "record" else "candidate")
    shutil.rmtree(root, ignore_errors=True)
    all_cases = cases()
    heavy = {k: v for k, v in all_cases.items() if k.startswith("full")}
    light = {k: v for k, v in all_cases.items() if k not in heavy}
    with ThreadPoolExecutor(3) as pool:
        list(pool.map(lambda item: run_case(binary, root, *item), light.items()))
    for name, args in heavy.items():
        run_case(binary, root, name, args)
    if mode == "record":
        print(f"recorded {len(all_cases)} cases, {len(collect(root))} files")
        sys.exit(0)
    expected, actual = collect(golden), collect(root)
    actual = {k: v.replace(root, golden) for k, v in actual.items()}
    differing = sorted(set(expected) ^ set(actual)) + sorted(
        k for k in expected.keys() & actual.keys() if expected[k] != actual[k])
    print(f"{len(all_cases)} cases, {len(expected)} files: " +
          ("IDENTICAL" if not differing else f"{len(differing)} DIFFER"))
    for name in differing[:25]:
        print("  ", name)
    sys.exit(1 if differing else 0)
