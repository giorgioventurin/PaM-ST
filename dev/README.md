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

Writes the fixtures the other scripts use into `dev/data/`: two synthetic
tissues and an 8000 x 8000 crop of LSP31891. Deterministic.

## golden.py — output equivalence

Runs 118 configurations (both metrics, tolerances from 0 to above the metric
maximum, overlapping and covering neighbourhoods, both nulls, FewRS, frozen
types, freeze sweeps, thread counts, a full tissue, 25 invalid-argument cases)
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
python3 dev/fwer.py csr 40
```

A correct search reports something in at most about 5% of `csr` datasets. This
is what caught the selection bias: choosing the candidate family and testing it
on the same neighbourhoods made every dataset with a candidate significant.

## power.py — detection

Plants 40 two-type niches of radius 150 and checks they are recovered.

```sh
python3 dev/power.py 12
```
