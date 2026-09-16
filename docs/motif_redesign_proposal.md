# Why the top motifs are blobs, and what to change

Status: **proposal only — no code has been modified.**
Reviewed: `src/`, `include/`, `docs/`, `output/` (LSP31891, LSP31894), `data/PDAC_202602`.
Diagnostics: run read-only, in a separate Python reimplementation of the
neighborhood/histogram step (see *Caveats* at the end).

## Summary

The blobs are not a tuning problem. Three parts of the design each push the
answer toward a single cell type, and they compound.

1. **The null model** (global label shuffle) tests whether the tissue has *any*
   spatial structure. That answer is always yes, which is why every stage of the
   freeze sweep reports `z` between 10 and 680 with `p_raw` pinned at the
   permutation floor.
2. **The pattern language** (a rho-ball around a direction in 20-dimensional
   composition space) has the simplex vertices as attractors. Six of the eight
   motifs the sweep reported sit *inside* the rho = 0.05 ball of a pure-type
   vertex.
3. **Ranking by frequency** selects against mixtures by construction. On
   LSP31891 at r = 100 the most frequent exact composition has support 3671 and
   is literally `unclassified:1`; the best composition with three types at >= 2
   cells each has support 53.

The encouraging part: with a null that keeps each cell type's own clumping and
only decouples types from each other, the top-ranked patterns in the same data
became `Bcells + CD8_Tcells + Other_immune` (lift 9.7x), `CD8 + Macrophages +
Mem_CD4` (8.2x), `CD8 + Macrophages + Treg` (5.8x). The signal is there; the
current statistic cannot see it.

---

## Diagnosis

### A. The global label shuffle tests the wrong hypothesis

`run_pattern_test_prepared` shuffles labels over all permutable cells, so every
cell type becomes completely spatially random. But every cell type in a tumor is
clustered — that is what tissue is. Any neighborhood statistic that responds to
clustering is therefore astronomically significant, and the fibroblast blob wins
by the largest margin because fibroblasts are both the most abundant and the
most clustered type.

The sweep is the proof: 20 stages, 20 motifs, all at `p_raw = 1e-4`, `z` from 10
to 680. A test that rejects everything carries no ranking information, and the
Holm correction across stages is then decorative — it is adjusting p-values that
are all at the resampling floor.

Measured (LSP31891, r = 100, exact composition matching, 95338 cells, 20 types),
support of the most frequent composition:

| condition                          | top support | note                          |
|------------------------------------|-------------|-------------------------------|
| observed                           | 3671        | top class = `unclassified:1`  |
| global shuffle (current null)      | ~240        | 3 resamples: 218, 241, 260    |
| block shuffle, 500-unit tiles      | ~2318       | 3 resamples: 2270, 2326, 2357 |
| block shuffle, 1000-unit tiles     | ~1743       | 3 resamples: 1719, 1746, 1763 |

Under the current null the observed blob looks impossible. Under a shuffle that
merely preserves composition inside 500-unit tiles — still destroying all
fine-scale structure — the same blob is within 1.6x of chance. Most of the
"significance" is tissue-scale composition heterogeneity, not a discovery.

### B. The rho-ball has the simplex vertices as attractors

Histograms are L2-normalized before matching, so a motif is a *direction*, and
rho = 0.05 is a chord distance — an angle of 2.87 degrees. Two consequences:

- **Size is discarded.** A neighborhood of 1 fibroblast and one of 60
  fibroblasts are the same point.
- **Mass concentrates at the vertices.** All near-pure neighborhoods pile onto
  the same tiny region of the sphere, while a genuine 5-type mixture has four
  free ratios that scatter its instances far apart in 20 dimensions. The
  rho-ball is a kernel density estimate with a 2.87-degree bandwidth in a
  19-dimensional space; the only places with enough mass are the vertices.

The argmax then drifts to a direction just inside the vertex ball, so that it
captures the whole vertex plus a rim. Chord distance from each reported motif to
its nearest pure vertex:

| freeze-sweep winner (LSP31894, r = 100)                                      | dominant | dist to vertex | inside rho=0.05 |
|------------------------------------------------------------------------------|----------|----------------|-----------------|
| `Other_fibroblasts:20; Prolif_Other_fibroblasts:1`                            | 95.2%    | 0.0500         | **yes**         |
| `Other_immune:64; CD4:1; Mem_CD4:2; Mem_CD8:1; Prolif_Mem_CD8:1`              | 92.8%    | 0.0413         | **yes**         |
| `Other_cancer:50; Other_fibroblasts:2; Prolif_Other_cancer:1`                 | 94.3%    | 0.0447         | **yes**         |
| `Macrophages:68; CD4_Tcells:1; Tcells:2`                                      | 95.8%    | 0.0329         | **yes**         |
| `CD4_Tcells:61; Endothelial:1; Macrophages:2; Tcells:1`                       | 93.8%    | 0.0401         | **yes**         |
| `Prolif_Bcells:74; Prolif_Mem_CD8:3; Other_immune:1`                          | 94.9%    | 0.0427         | **yes**         |
| `Bcells:58; Prolif_Mem_CD8:5; Other_immune:3`                                 | 87.9%    | 0.1002         | no              |
| `Other_cancer:58; Other_fibroblasts:14; Epi_cancer:1; Other_immune:1` (r=300) | 78.4%    | 0.2374         | no              |

Row 1 is the clearest case: at 0.049954 it sits just inside the 0.05 ball of the
pure-fibroblast vertex, so its match set contains every purely-fibroblast
neighborhood in the sample regardless of size, plus everything >= 95% fibroblast.
The r = 300 run escapes the vertex trap — its winner is a genuine 78/19
cancer-fibroblast mixture — which is the "blob with a rim".

Increasing rho does not fix this. It grows the vertex balls faster than any
interior mode, because that is where the mass is.

### C. Frequency and interest are anti-correlated in this data

At r = 100 on LSP31891: 9.4% of neighborhoods are 100% one type, 13.4% are >= 90%
one type, mean neighborhood size 16.3 (median 12, max 90), mean 3.56 distinct
types. Of 50162 distinct exact compositions, the three most frequent are
`unclassified:1` (3671), `Other_fibroblasts:1` (3087), `Other_cancer:1` (1622) —
single cells. The most frequent composition containing three types with >= 2
cells each has support 53. No re-ranking inside a frequency-first search closes a
70x gap.

It goes deeper. Scoring presence-based type triples against a *block-permutation*
null, **every** observed triple came back depleted (lift 0.31-0.70). That is not
a bug: real tissue is more segregated than a reshuffle, so "mixed neighborhoods
are more frequent than chance" is a hypothesis this tissue actively contradicts.
Frequent-mixture mining is fighting the dominant signal. The question that does
have an answer is co-occurrence *given* each type's own clustering.

### D. Freeze-by-abundance is whack-a-mole

Freezing the top type removes it as a *center* but leaves it in every histogram,
and the next-most-abundant type immediately wins the same way: fibroblast blob ->
immune blob -> cancer blob -> macrophage blob -> CD4 blob -> proliferating-B blob
-> B blob. Stages 3-6, 7-12 and 13-17 return literally the same motif, because
freezing a type that was not in the winning pattern changes nothing about that
pattern. Cost: 8852 s for 20 stages, for information stage 0 already contained.

---

## Proposed changes

### P0 — stop the artifacts (small edits, current architecture)

**01. Drop `unclassified` from the pattern alphabet.**
19.7% of LSP31891 cells are `unclassified`, and it is the single most frequent
motif in the data. It is not a cell type and should not be able to define, or sit
inside, a motif. Keep those cells as spatial context (they occupy space and
displace neighbors) but exclude their column from the composition vector. The
same argument applies more weakly to the `Other_*` catch-alls — worth a
sensitivity run with and without.
*Where:* new `--exclude-label NAME` in `cli.cpp`; mask the column in
`histograms_*_into`.

**02. Minimum neighborhood size, and switch to k-nearest neighbours.**
The target is groups "in the order of 10". A fixed radius does not deliver that:
r = 300 averages 55 cells, r = 100 averages 16 but ranges 1-90, and the small
ones are pure by construction and feed the vertex classes. Use kNN with k ~ 10-15
(optionally capped by a max radius so sparse regions do not produce absurd
neighborhoods). Every neighborhood then has the same denominator, removing size
as a confounder in a scale-invariant statistic.
*Where:* `build_neighbors()` — add a kNN mode; reuse the existing grid with a
bounded heap per center.

**03. Constrain the candidate family to multi-type patterns — in the null too.**
Require >= m types with >= c cells each (start at m = 3, c = 2). This must be
applied identically when mining the observed data *and* inside every permutation,
or the max-statistic calibration is wrong: the filter goes inside `top_motifs()`,
not after it. This alone changes what wins; it does not fix why the winner is
uninteresting (that is P1).
*Where:* `top_motifs()` / `select_top_motif_indices()`, mirrored in the
`null_motifs` path of `run_pattern_test_prepared`.

**04. Report the class centroid, not an arbitrary member.**
`composition_classes()` stores the raw histogram of the *first* row it
encounters as the class representative, and that is what `pattern.csv` prints.
`Macrophages:68; Tcells:2` is not the motif — it is one neighborhood that
happened to come first in index order, and a different center ordering would
print a different-looking motif for the same class. Report the mean proportion
vector over matched neighborhoods, the number of matches, and the distribution of
matched neighborhood sizes. A reader seeing `Prolif_Bcells:74` will reasonably
think you found 74 proliferating B cells inside one 100-unit disc.
*Where:* `composition_classes()` in `utils.cpp`, `write_outputs()` / `pattern.csv`.

**05. Rank by surprise, not by count.**
The selected motif is currently `argmax` of raw frequency, with the permutation
test applied afterwards. Invert it: compute a per-pattern null distribution, rank
by standardized excess (or observed/expected lift), then take the max of the
*standardized* statistic across patterns within each permutation for the
family-wise calibration. That is a Westfall-Young max-T, it is a small change to
the loop that already exists, and it is what lets a rare-but-consistent triad
beat an abundant blob.
*Where:* `run_pattern_test_prepared()` — keep per-pattern null vectors,
standardize before the max.

**06. Retire freeze-by-abundance, or redefine it.**
As implemented it costs 20x the runtime and returns the same handful of motifs.
If you want the concept — "what is left once the dominant architecture is
accounted for" — the right version is to *condition*, not to freeze centers:
restrict to neighborhoods whose dominant type is X and ask which minority types
co-occur beyond expectation within that stratum. That directly extracts the "rim"
biology the r = 300 blob is hinting at.
*Where:* `run_freeze_sweep()`, `docs/cumulative_freeze_sweep.md`.

### P1 — replace the null (the decision that matters)

**07. Use a null that preserves each cell type's own spatial clustering.**
The scientific question is "do these types co-occur beyond what each type's own
distribution implies?", so the null must keep each type clumped and destroy only
the registration between types. Three options, increasing fidelity:

- *Block permutation.* Shuffle labels only within tiles of side T. One
  parameter, ~15 lines, preserves composition above scale T. This is what
  CRAWDAD does, swept over T so the scale at which a relationship appears becomes
  part of the result rather than a hidden assumption. Caveat: within a tile it
  spreads a clumped type out, so it *over-predicts* mixing — which is why the
  triples all came back depleted against it.
- *Smoothed-intensity resampling (recommended).* Estimate each type's intensity
  lambda_t(x) with a kernel of bandwidth h, then resample labels for the fixed
  cell positions with P(label = t) proportional to lambda_t(x_i), using a
  Gumbel-top-k assignment so exact per-type totals are preserved. No tile
  artifacts, works with an irregular tissue mask, and h becomes the explicit
  statement of "structure coarser than this is given, finer than this is a
  discovery".
- *Per-type random shift.* The classical Lotwick-Silverman construction:
  translate each type's point set by an independent random vector. Highest
  fidelity to within-type structure, but needs care with the tissue mask, since
  shifted points can land off-tissue.

*Where:* replaces `std::shuffle(shuffled_labels...)` in
`run_pattern_test_prepared`; everything downstream is unchanged.

What the same data yields with a prototype of the second option (LSP31891,
r = 100, presence-based triples, `unclassified` excluded, 10 resamples):

| triple                                    | support | null mean | lift  | z    |
|-------------------------------------------|---------|-----------|-------|------|
| Bcells + CD8_T + Other_immune             | 5230    | 538       | 9.7x  | 19.5 |
| CD8_T + Macrophages + Mem_CD4             | 5041    | 613       | 8.2x  | 18.2 |
| CD8_T + Mem_CD4 + Other_immune            | 8084    | 984       | 8.2x  | 20.8 |
| Bcells + CD8_T + Mem_CD4                  | 4188    | 535       | 7.8x  | 10.5 |
| CD8_T + Macrophages + Other_immune        | 7536    | 1115      | 6.8x  | 16.4 |
| Macrophages + Mem_CD4 + Other_immune      | 7859    | 1195      | 6.6x  | 12.8 |
| CD8_T + Macrophages + Treg                | 2081    | 357       | 5.8x  |  9.7 |
| CD8_T + Mem_CD4 + Treg                    | 2358    | 473       | 5.0x  | 10.6 |
| *Other_fibroblasts + Other_immune + CD8_T*| *12760* | *5001*    | *2.6x*| *8.3*|

74 of 461 triples with support >= 50 were never exceeded in any resample; 154 were
significantly depleted. Note the last row: the *most frequent* triple ranks 9th
by lift — exactly the inversion that makes frequency-first search fail. The top
triads are recognisable biology (B/CD8/CD4 co-occurrence is a tertiary lymphoid
structure; CD8 + macrophage + Treg is a suppressed myeloid-T niche).

### P2 — change the pattern language, and get exact multiple-testing control

**08. Patterns are type *sets*, not points in 20-dimensional histogram space.**
"Groups of different cell types frequently close to each other" is naturally an
itemset, not an exact count vector. Define a neighborhood's pattern as the set of
types present at >= c cells; support of a set S is the number of neighborhoods
whose present-set contains S. Anti-monotone, interpretable, robust to exact
ratios, and it removes the curse-of-dimensionality problem: supports go from 53
to the thousands. Run it at two or three values of c and include them in the
tested family.

**09. With K = 20, enumerate the whole lattice with a zeta transform.**
No FP-growth or KD-tree needed. Represent each neighborhood's present-set as a
20-bit mask, histogram the masks into a dense array of 2^20 counters (4 MB), then
run the superset zeta transform:

```
for each bit b:
    for each mask m without bit b:
        cnt[m] += cnt[m | b]
```

20 * 2^20 ~ 21M operations gives the exact support of **every** one of the
1048576 possible type sets in about 20 ms. Repeat per permutation and you have
the complete null distribution for the entire pattern lattice, which makes
Westfall-Young max-T over all patterns exact rather than restricted to a top-K
shortlist — and removes the need for the 2*rho packing rule and the FewRS
order-statistic machinery.
*Where:* replaces `KDTree.cpp`, `VPTreeJS.cpp`, `composition_classes()`;
~150 lines total.

**10. Layer Tarone pruning on top for power.**
Once patterns are itemsets with a support-based test, Tarone's minimum attainable
p-value prunes untestable patterns before correction, so the effective family is
far smaller than 2^20 and the corrected threshold is much less punishing than
Bonferroni. (Significant-pattern-mining literature: LAMP; Westfall-Young light.)
Plugs directly into the zeta-transform structure.

**11. Make cross-sample replication a first-class criterion.**
There are four samples. A pattern enriched in 4/4 is worth more than any p-value
from one slide, and reviewers will ask. Report per-sample lift and a combined
statistic, and treat sample-level consistency as part of the discovery rule
rather than a post-hoc check. With the zeta transform this is essentially free —
one 2^20 array per sample.

### P3 — code-level notes worth fixing regardless

- **Threads are created and joined inside the permutation loop.**
  `histograms_for_rows_into` and `top_motifs` each spawn `hardware_concurrency`
  threads per call, so a 1000-permutation run creates ~1e5 threads. Permutations
  are independent — parallelize the *outer* loop with per-permutation seeds
  (`seed + b`, or a counter-based RNG) and keep reproducibility.
- **`composition_classes()` allocates a `std::vector<int>` key per row** — 95k
  heap allocations per permutation. Hash the row in place, or use the bitmask
  representation from change 09.
- **Neighbours are a vector-of-vectors.** A CSR layout (offsets + one flat index
  array) removes 95k allocations and makes the permutation-loop histogram pass
  cache-friendly.
- **The gcd equivalence is arithmetically arbitrary.** `composition_key` merges
  (2,4) with (1,2) but not (3,5) with (2,3); with rho > 0 the ball query re-merges
  everything anyway, so the gcd step mostly deduplicates identical rows at the
  cost of a gcd per element. It also makes size-invariance a property of the
  *exact*-match path only.
- **KD-tree in 20 dimensions.** Bounding-box pruning degrades badly at this
  dimensionality, and the tree is rebuilt from scratch every permutation.
  Sparsity would help (mean 3.6 non-zero columns), but change 09 removes the need.
- **Covering mode is recomputed per freeze stage**, so the number of tested
  neighborhoods differs between stages and frequencies are not comparable across
  the sweep.
- **Every center is its own neighbour**, and centers are restricted to non-frozen
  types, so each candidate neighborhood is guaranteed >= 1 cell of its own
  center's type. Combined with small neighborhoods this is a direct push toward
  pure-center-type motifs. Consider excluding the center, or conditioning on it
  explicitly.
- **`Result::null_max` holds per-rank counts, not maxima** — a naming trap.
- **Source and outputs have diverged.** `output/*/report.txt` reports
  `density_sampling_enabled`, `density_knn_k`, `preprocess_seconds` and 19 labels
  / 63892 cells for LSP31891; there is no density sampling in `src/`, and the raw
  file has 20 types and 95338 cells. Worth reconciling before any of these
  numbers go in the paper. `code/PaM-ST` is also not under version control,
  unlike `code/MiMo-ST`.

---

## Positioning relative to published work

- **SpatialQuery** (An et al.) already does frequent-pattern mining of cell-type
  itemsets over kNN / radius neighborhoods with FP-growth — i.e. change 08,
  published. The itemset framing is not where the novelty can live.
- **CRAWDAD** (*Nat Commun*, 2025) tests cell-type relationships with grid-based
  label shuffling swept across length scales, with offset grids to suppress tile
  artifacts — change 07 option one, published; a good citation and baseline.
- What is *not* settled is the combination this project is placed to deliver: a
  null that preserves each type's own spatial process rather than shuffling
  labels globally (squidpy's `nhood_enrichment`, imcRtools' `testInteractions`
  and most niche-enrichment tools use a global or near-global permutation — the
  same weakness as section A), combined with *exact* family-wise control over the
  entire pattern lattice rather than a top-K shortlist, plus multi-sample
  replication as part of the discovery rule. The zeta-transform trick is what
  makes that computationally trivial where others needed heuristics.

One reframing worth considering for the paper: given the depletion result in
section C, "frequent co-occurrence" may be the wrong headline. *Which cell types
are found together more often than their individual spatial distributions
predict, and at what length scale* is a sharper question, it is what the data
supports, and it turns the segregation findings (154 depleted triples) into a
result rather than a nuisance.

---

## Caveats on these diagnostics

All diagnostics were run read-only on LSP31891 at r = 100, using a separate
Python reimplementation of the neighborhood/histogram step — not through
`pam_st` itself. Exact composition matching (rho = 0) was used for the
class-support numbers, so they are not identical to what rho = 0.05 would give.
The vertex-distance table is computed directly from the motifs in the existing
`report.txt` files and is exact. The field-shift null in P1 is a crude prototype
(100-unit binned, lightly smoothed, toroidal, no tissue mask), so the absolute
lifts are inflated; the ranking and the direction of the effect are the parts to
trust. Nothing in `code/PaM-ST` was modified.

---

**Implemented.** A ground-up C++ implementation of P0-P2 lives in
`code/CoMo-ST` (see its `README.md` and `docs/method.md`). First run and
validated results: `code/CoMo-ST/docs/first_run_2026_09_09.md`.
