# PaM-Mixed

A separate, dependency-free C++17 implementation for recurring mixed cell-type
groups, with Python plotting and reproducible synthetic examples. The original
PaM-ST implementation and outputs are not modified.

## Build and test

From this directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build --output-on-failure
python -m pip install -r requirements.txt
```

The C++ executable needs only the standard library and threads. Plotting needs
NumPy, pandas and Matplotlib. Synthetic generation and integration tests use
only Python's standard library. On this machine an existing working plotting
interpreter is `/Users/giorgio/anaconda3/envs/GNN4ST/bin/python`.

## Run on your data

```sh
./build/pam_mixed \
  --input ../../data/PDAC_202602/annotated_cells/LSP31891_P181_rep1_annotated.csv \
  --output results/LSP31891_blocks \
  --k 10 --max-radius 100 --min-types 2 --max-dominance 0.7 \
  --max-relabels 1 --min-support 10 --max-motifs 12 --candidate-pool 100 \
  --null blocks --block-size 1000 --region-size 1000 \
  --permutations 199 --seed 37 --threads 4
python scripts/plot_results.py --results results/LSP31891_blocks --top 6
```

This run is already included under `results/LSP31891_blocks`. Use a different
output directory for a new run, or explicitly pass `--overwrite` to replace
generated tables. When overwriting, regenerate plots to update their content.
The numerical radii above are in **input coordinate units**, not assumed microns.
The 1000-unit block size is an exploratory background scale, not a fitted or
biologically validated compartment definition.

## Input

CSV must contain `X_centroid,Y_centroid,Cell_Type`. Alternative names can be
specified with `--x-column`, `--y-column`, and `--type-column`.

| Optional column | Meaning |
| --- | --- |
| `Cell_ID` | Unique within a sample; otherwise a record-based ID is generated |
| `Sample` | Specimen/image ID; neighborhoods and shuffles never cross samples |
| `Component` | Disconnected tissue component; neighborhoods and shuffles never cross components |
| `Stratum` | Externally supplied compartment for stratified shuffling; neighborhoods may cross its boundary |

Without `Sample`, the input filename identifies the sample. Samples are not
automatically patients or independent biological replicates. `Unclassified`
labels are excluded case-insensitively and their number is recorded. Malformed
coordinates, empty labels/grouping values, and duplicate IDs within a sample
are rejected. Quoted CSV fields, CRLF and UTF-8 BOMs are supported.

`Component` is a hard categorical boundary, not a polygon mask. A physical
cutoff alone cannot recognize holes, vessel walls, or paths through background;
annotate components appropriately when those barriers matter.

## Method

1. Find the center and its `k-1` nearest cells in its sample/component. Reject
   centers with fewer than `k` cells within `max-radius`. The center is explicitly
   included even when coordinates coincide. The group diameter may be twice the
   center cutoff.
2. Deduplicate identical cell sets, retaining the lowest input-index center.
   Order groups by retained-center radius, then center index. All geometry is
   fixed independently of the labels and reused across permutations.
3. Calculate integer cell-type counts, Shannon entropy using natural logs,
   effective type count `exp(entropy)`, and dominant fraction. Require at least
   `min-types` and dominance no larger than `max-dominance`. These filters apply
   to **every matched occurrence** and are recomputed for each shuffled dataset.
4. Match equal-size counts with distance `sum(abs(c-c'))/2 <= max-relabels`.
   Supported tolerances are 0, 1 and 2. By default the set of represented types
   must be identical; `--allow-type-turnover` relaxes this requirement. Counts
   retain multiplicities; proportional vectors of different sizes are never
   pooled. An inverted subhistogram index avoids comparing every template with
   every neighborhood.
5. Count all unique matched cell sets (`raw_matches`) and greedily select
   cell-disjoint occurrences in the fixed geometry order (`support`). This is
   a reproducible maximal packing, not a maximum-cardinality optimization.
   Groups within a motif share no cells; different motifs can share cells.
   Spatial dependence remains, and occurrences are not biological replicates.
6. Search all eligible observed composition templates. The global statistic is
   their maximum disjoint support. The search is repeated in every permutation,
   including diversity and any contact filters. `min-support` screens reported
   candidates, but maxima retain their actual values even below that floor.
7. For descriptive prioritization, keep the top `candidate-pool` recurring
   templates by support (0 means all). Evaluate these **same templates** in each
   shuffled dataset and rank by observed minus mean null support. Report up to
   `max-motifs`, suppressing near variants at distance <= twice the matching
   tolerance when their type sets match. This suppression changes reporting
   only, not the search used to calibrate the global test. Setting pool=0 costs
   additional time and memory. The finite pool can omit rare enriched motifs.

The greedy support, radius, type restrictions, matching tolerance, and background
scale define the target. Choose them before confirmatory testing. Changing them
after inspecting p-values is exploratory unless the larger search is calibrated.

## Null models and inference

All implemented nulls condition on positions and shuffle labels, preserving type
counts in the indicated groups:

| `--null` | Shuffle groups |
| --- | --- |
| `global` | Sample x Component |
| `stratified` | Sample x Component x Stratum; requires `Stratum` column |
| `blocks` (default) | Sample x Component x Stratum x spatial tile |

Tiles use `floor(x/block-size), floor(y/block-size)`, anchored at zero. The
`region-size` tiles used to describe dispersion are separate from shuffle groups.
Strata with just one label cannot change; `run.json` reports how many cells and
strata admit label exchange. Shuffles include cells without an eligible
neighborhood, because their labels are still part of the tissue population.

With `B` permutations, the discovery test is

```text
T = maximum cell-disjoint support across all eligible observed templates
p_global = (1 + number of null maxima >= observed T) / (B + 1)
```

Each replicate starts from the original labels with its own deterministic random
seed. Worker count does not change results on the same C++ runtime. The p-value
resolution is `1/(B+1)`; 199 gives 0.005. More permutations improve resolution,
not the suitability of a null. C++ library differences may change generated
random permutations across platforms.

The global test is valid under the specified conditional random-label
exchangeability assumption. `p_search` compares a reported motif's support with
the same complete-search maxima. It describes evidence against that complete
null; it is **not a motif-specific FDR guarantee under partial alternatives**.
Discovery `p_value` and `p_value_holm` are deliberately blank. Enrichment is
descriptive after selection, with a recorded 0.5 pseudocount:
`(observed_support + 0.5)/(expected_support + 0.5)`.

Within-compartment shuffling retains coarse composition but destroys homotypic
clustering inside compartments. Thus none of these nulls establishes association
beyond arbitrary pre-existing clustering. A clustering-preserving spatial model
and irreducible higher-order interaction test are **not implemented**. Such a
model requires explicit assumptions, simulation calibration and uncertainty
treatment; substituting an unchecked MCMC shuffle would not provide exact tests.

### Fixed-template validation

Discover templates in separate specimens, freeze the geometry/filter/tolerance
settings, and test on held-out specimens:

```sh
./build/pam_mixed --input heldout.csv --output results/validation \
  --templates results/discovery/motif_composition.csv \
  --k 10 --max-radius 100 --min-types 2 --max-dominance 0.7 \
  --max-relabels 1 --null stratified --permutations 999 --threads 4
```

The template CSV contains `motif_id,cell_type,count`; absent validation types
receive observed count zero. All supplied templates, including zero-support
ones, are tested. Validation scores fixed templates in every shuffle, returns
their plus-one p-values, and applies Holm correction across **all** of them.
The minimum support requirement is an additional reporting decision, not a way
to exclude hypotheses from the correction. For the first Holm rejection among
M hypotheses, permutation resolution must reach alpha/M.

The program cannot verify that specimens are independent or templates/settings
were prespecified. Random cell splits leak through overlapping neighborhoods.
Holm handles dependence among valid p-values; it does not repair invalid null
assumptions or turn a random-label question into a mechanistic interaction test.

### Optional local contact feature

`--contact-radius D --min-cross-edge-fraction F` constructs all unordered pairs
within distance D among each group's cells and requires at least fraction F to
join different types. Edgeless groups fail a positive requirement. This feature
is recomputed after shuffling and can help distinguish mixing from segregation
at a chosen scale. It does not prove physical contact or higher-order interaction.

## Outputs and plots

| File | Contents |
| --- | --- |
| `run.json` | Complete settings, counts, assumptions, global p, timing |
| `cells.csv` | Classified cells and stable indices |
| `neighborhoods.csv` | Retained geometry, size, diversity, eligibility |
| `motifs.csv` | Support, spatial dispersion, descriptive enrichment and mode-specific tests |
| `motif_composition.csv` | Counts for reusable templates |
| `occurrences.csv`, `members.csv` | Selected disjoint groups and their actual cell identities |
| `null.csv` | Complete-search maxima and eligible counts for every permutation |
| `motif_null.csv` | Same-template support for each reported motif and permutation |

`scripts/plot_results.py` writes PNG and PDF versions of a composition/support
overview, null distributions, size/diversity diagnostics, specimen recurrence,
spatial maps, and an atlas of actual ten-cell occurrences. It also writes a
shared color palette, captions and plot manifest. See
[plotting details](scripts/PLOTTING.md). Zero-motif cases are handled explicitly.

## Reproduce the controls

```sh
python scripts/run_demo.py --binary build/pam_mixed \
  --output examples --permutations 199 --threads 4
```

This generates random labels, pure blobs, adjacent domains, and planted mixed
groups, then freezes discovery templates and validates them on a separately
generated planted dataset. All outputs and figures are retained. Use
`--skip-plots` if the plotting dependencies are unavailable. The included run
uses 1,800 cells per scenario and 72 planted 4A+3B+3C groups. These examples test
known behavior; they are not a repeated-simulation false-positive-rate study.

The test suite additionally compares indexed scoring and exact kNN against
brute-force references, including duplicate coordinates, negative coordinates,
contact filtering and sample/component boundaries. CLI tests exercise CSV
handling, no-motif cases, absent validation types, conditional nulls,
permutation arithmetic and thread-count reproducibility.

See [the completed run notes](RESULTS.md) for observed outcomes.
