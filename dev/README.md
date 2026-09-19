# Validation harness

Scripts that check PaM-ST beyond the CTest suite. They need only the Python
standard library and a built `pam_st`; nothing here is used at analysis time.

Build natively (the `cmake` on this Mac is an x86 binary and cannot run):

```sh
clang++ -std=c++20 -O3 -DNDEBUG -Iinclude src/*.cpp -o dev/pam_st
```

For the CTest suite, install a native CMake once:

```sh
python3 -m venv dev/.venv && dev/.venv/bin/pip install cmake
dev/.venv/bin/cmake -S . -B build-native -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++
dev/.venv/bin/cmake --build build-native -j 10 && dev/.venv/bin/ctest --test-dir build-native --output-on-failure
```

## make_datasets.py

Writes the fixtures the other scripts use into `dev/data/`: synthetic
tissues, and 8000 x 8000 crops of LSP31891, LSP31894, LSP31895 and LSP31896
centred on each sample's median cell. Deterministic.

## golden.py — output equivalence

Runs 131 configurations (both metrics, tolerances from 0 to above the metric
maximum, overlapping and covering neighbourhoods, both nulls, FewRS, frozen
types, freeze sweeps, thread counts, a full tissue, multi-sample runs,
invalid-argument cases)
and compares every output file, plus stdout, stderr and exit codes. Timing
fields are masked; nothing else may differ.

```sh
python3 dev/golden.py record dev/pam_st_reference   # baseline from a known-good binary
python3 dev/golden.py check  dev/pam_st             # after a change
```

Use it for any change that is meant to preserve behaviour. A refactor must come
back IDENTICAL; a deliberate change to the usage text will show up as
differences in the `err_*/stderr` cases only.

## fwer.py — false positives

Runs the significance search on tissues with no local structure, where nothing
should be reported. `csr` is uniform positions with i.i.d. labels (the null is
exactly true); `smooth` drifts composition across the tissue.

```sh
python3 dev/fwer.py csr 60 2 2999    # scenario, datasets, min-support, permutations
```

A correct search reports something in at most about 5% of `csr` datasets.
The check is only informative when the permutations can reach the level at
all: the smallest adjusted p is about family size / (permutations + 1), so use
roughly 20 permutations per candidate. This
is what caught the selection bias: choosing the candidate family and testing it
on the same neighbourhoods made every dataset with a candidate significant.

## power.py — detection

Plants 40 two-type niches of radius 150 and checks they are recovered.

```sh
python3 dev/power.py 12
```

## multisample.py — several samples in one run

```sh
python3 dev/multisample.py fwer 60 5 2999    # datasets, min-support, permutations
python3 dev/multisample.py replication
```

`fwer` passes three structureless tissues as repeated `--input` files and
counts datasets where anything is reported. `replication` plants niches in
three, one or none of three samples; the pooled test should find them in the
first two cases, with `replicated_in` equal to the number of planted samples.

The four PDAC crops are a realistic multi-sample case:

```sh
./dev/pam_st --input dev/data/LSP31891.csv --input dev/data/LSP31894.csv \
  --input dev/data/LSP31895.csv --input dev/data/LSP31896.csv \
  --radius 100 --rho 0.05 --metric l2 --statistic minp --min-support 10 \
  --min-types 2 --min-type-cells 2 --split-size 1500 --max-motifs 20 \
  --permutations 14999 --seed 37 --threads 10 --null-model block --block-size 500 \
  --output-dir dev/multisample_out/pdac4_crops
```

It forms a family of 668 and takes about 2.5 minutes; `docs/pipeline.tex`
discusses the result.
