# Completed runs

Run date: September 9, 2026. No parameters were changed in response to the
observed p-values. The original PaM-ST outputs remain untouched.

## Real tissue: LSP31891

Input: `LSP31891_P181_rep1_annotated.csv` in the original PDAC annotated-cell
directory. The implemented analysis retained 76,563 classified cells and excluded
18,775 Unclassified cells. It ran in 20.81 seconds using four permutation workers.

Settings: ten cells including the center, center distance at most 100 input
coordinate units, at least two types, dominant fraction <=0.7, one allowed label
replacement with the represented type set preserved, and at least ten
cell-disjoint occurrences for reported candidates. The null shuffled within
1000-unit spatial tiles (199 replicates, seed 37). Distances have not been
converted to microns. No contact-edge filter was used.

| Quantity | Result |
| --- | ---: |
| Centers rejected for insufficient cells within the cutoff | 29,726 |
| Duplicate cell sets removed | 3,414 |
| Distinct neighborhoods retained | 43,423 |
| Neighborhoods passing mixed-type eligibility | 28,528 |
| Observed eligible composition templates searched | 10,367 |
| Templates meeting minimum spatial support | 1,922 |
| Support-screened templates scored for descriptive enrichment | 100 |
| Reported templates after variant suppression | 12 |
| Maximum observed disjoint support | 217 |
| Global random-label search p-value | **0.89** |
| Reported motifs meeting the search threshold at alpha=0.05 | **0** |

The run finds mixed groups rather than homogeneous blobs, but their maximum
recurrence does **not** reject this particular background model. Descriptive
same-template enrichment is different from the complete-search significance
test. For example:

| Template | Disjoint support | Raw matched cell sets | Same-template mean null support |
| --- | ---: | ---: | ---: |
| 4 Epithelial_cancer + 6 Other_cancer | 184 | 725 | 11.73 |
| 6 Other_fibroblasts + 4 Other_immune | 217 | 658 | 103.71 |
| 7 Epithelial_cancer + 3 Other_cancer | 104 | 316 | 2.59 |

Matches allow one label replacement and must still meet all composition
requirements. These comparisons were selected from this tissue and are
descriptive; none is a confirmed motif-specific discovery. Unstandardized
maximum support may have low power for a comparatively rare but enriched
template because another composition can dominate the null maximum. A future
calibrated enrichment-based search would need its selection and score
normalization included in resampling, or independent template validation.

Explore the complete results in:

- [Overview](results/LSP31891_blocks/plots/overview.png)
- [Actual-cell occurrence atlas](results/LSP31891_blocks/plots/occurrence_atlas_01.png)
- [Spatial recurrence map](results/LSP31891_blocks/plots/spatial_01_LSP31891_P181_rep1_annotated.png)
- [Neighborhood diagnostics](results/LSP31891_blocks/plots/neighborhood_diagnostics.png)
- [Motif table](results/LSP31891_blocks/motifs.csv)
- [Settings and inference scope](results/LSP31891_blocks/run.json)

PDF equivalents and figure captions are included in the plots directory.

## Additional tissues: LSP31894, LSP31895 and LSP31896

All three requested datasets completed with exactly the LSP31891 settings above,
including 199 permutations, seed 37, and four workers. Each input was the
corresponding `LSP3189X_P181_rep1_annotated.csv` from the same PDAC annotated-cell
directory. Each run exported 12 descriptive templates and six PNG/PDF figure
pairs. The table includes LSP31891 as a reference.

| Dataset | Classified cells | Eligible mixed neighborhoods | Maximum disjoint support | Mean null maximum | Global p | Holm p across these four datasets |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| LSP31891 | 76,563 | 28,528 | 217 | 226.49 | 0.89 | 1.00 |
| LSP31894 | 89,301 | 32,292 | 341 | 525.88 | 1.00 | 1.00 |
| LSP31895 | 120,202 | 53,179 | 296 | 190.75 | **0.005** | **0.020** |
| LSP31896 | 93,877 | 59,296 | 211 | 170.21 | **0.005** | **0.020** |

LSP31895 and LSP31896 reject the specified complete blockwise random-label null;
LSP31894 does not. A p-value of 0.005 is the smallest attainable with 199
permutations: none of those null maxima reached the observed maximum in these
two datasets. The Holm correction covers these four dataset-level tests and
does not provide motif-specific error control. It does not require the files to
represent independent patients. Raw support should not be compared as a
normalized effect size across tissues with different cell counts and composition.

Examples from the reported templates:

| Dataset / motif | Template composition | Disjoint support | Same-template mean null support | Complete-null search p |
| --- | --- | ---: | ---: | ---: |
| LSP31894 / M001 | 6 Other_cancer + 4 Proliferating_Other_cancer | 186 | 11.60 | 1.00 |
| LSP31895 / M001 | 6 Epithelial_cancer + 4 Other_cancer | 252 | 10.29 | 0.005 |
| LSP31895 / M002 | 6 Endothelial_cells + 4 Other_fibroblasts | 273 | 44.46 | 0.005 |
| LSP31895 / M003 | 3 Endothelial_cells + 7 Other_fibroblasts | 251 | 66.60 | 0.005 |
| LSP31896 / M001 | 4 Epithelial_cancer + 6 Other_cancer | 211 | 48.87 | 0.005 |

Support includes the permitted one-label count variation while preserving the
represented type set and eligibility rules. Occurrences of the same motif are
cell-disjoint; different motifs may share cells. Five of LSP31895's 12 reported
templates, one of LSP31896's, and none of LSP31894's cross the complete-null
search threshold at 0.05. These are not counts of confirmed motif-specific
discoveries. Same-template expectations remain descriptive after selection.
Templates are reported by excess support after variant suppression, so the
leading reported template need not be the template with maximum raw support.

The endothelial–fibroblast examples in LSP31895 are particularly relevant to
the requested mixed-type patterns. Several other results distinguish cancer or
proliferation annotations; whether those annotations represent the desired
biological diversity requires interpreting the supplied cell-type taxonomy.
The test still does not establish association beyond within-tile homotypic
clustering or prove higher-order interactions.

Artifacts:

- [Four-dataset null comparison](results/comparison/null_comparison.png),
  [PDF](results/comparison/null_comparison.pdf),
  [summary CSV](results/comparison/dataset_summary.csv), and
  [settings and interpretation JSON](results/comparison/dataset_summary.json).
- LSP31894: [overview](results/LSP31894_blocks/plots/overview.png),
  [actual occurrences](results/LSP31894_blocks/plots/occurrence_atlas_01.png),
  [motif table](results/LSP31894_blocks/motifs.csv).
- LSP31895: [overview](results/LSP31895_blocks/plots/overview.png),
  [actual occurrences](results/LSP31895_blocks/plots/occurrence_atlas_01.png),
  [motif table](results/LSP31895_blocks/motifs.csv).
- LSP31896: [overview](results/LSP31896_blocks/plots/overview.png),
  [actual occurrences](results/LSP31896_blocks/plots/occurrence_atlas_01.png),
  [motif table](results/LSP31896_blocks/motifs.csv).

Run metadata and additional spatial maps and diagnostics are included in each
result directory. C++ analysis took approximately 25, 81 and 101 seconds for
LSP31894, LSP31895 and LSP31896 respectively; concurrent execution means these
times are not controlled performance comparisons. Unclassified cells excluded
were 30,550, 53,611 and 33,158 respectively.

## Synthetic controls

Each dataset has 1,800 cells. Discovery and validation data were generated with
different seeds (37 and 38). All runs used 199 shuffles, a 20-unit cutoff, and
the same two-type/dominance/matching criteria; the reporting support floor was
five. Synthetic coordinates are arbitrary units. The null shuffles labels
within sample/component.

| Scenario | Eligible groups | Maximum support | Statistical outcome |
| --- | ---: | ---: | --- |
| Random labels | 180 | 13 | Global p=0.805 |
| Homogeneous blobs | 0 | 0 | No eligible motifs; global p=1 |
| Adjacent homogeneous domains | 600 | 16 | Global p=1 |
| Planted 4A+3B+3C groups | 180 | 72 | Planted template recovered; global p=0.005 |
| Independent planted validation | 180 | 72 | Planted template p=0.005; Holm-adjusted p=0.025 across 5 fixed templates |

The planted template recovered all 72 planted groups in both generated datasets.
The adjacent-domain case supplies a useful boundary control; its non-rejection
here does not establish that all such boundaries will be rejected as candidates
or that the method distinguishes arbitrary mixing from segregation.

- [Benchmark table](examples/benchmark_summary.csv)
- [Synthetic discovery overview](examples/results/planted_mixed/plots/overview.png)
- [Synthetic validation overview](examples/results/validation_planted/plots/overview.png)
- [Generated ground truth](examples/data/manifest.json)

These are individual controlled examples, not a simulation estimate of the
false-positive rate across many independently generated datasets.

## Verification

The Release build passed both CTest suites: 18 CLI integration cases and a C++
reference-comparison suite containing 60 randomized geometry cases and 360
matching/mining comparisons. The core comparisons cover all supported matching
tolerances, type-set modes, duplicate coordinates and spatial boundaries.
Generated figures were inspected visually, including real-tissue composition,
null distribution, spatial occurrence maps, and actual group membership.

The main remaining methodological limitation is the random-label null: spatial
tiles preserve coarse composition, not within-tile homotypic clustering. These
results should not be presented as association beyond all background clustering
or as irreducible higher-order interactions.
