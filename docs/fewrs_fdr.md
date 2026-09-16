# Few-shot FDR control for PaM-ST

## Goal

`fewrs-fdr` extends the FewRS maximum-statistic idea from family-wise error
control to a conservative false-discovery-rate guarantee for the ranked motif
frequencies mined by PaM-ST. It does not apply BH to the rank-specific
permutation p-values: lower motif ranks are shifted when a higher-ranked motif
is non-null, so that shortcut does not provide motif-wise FDR control under
mixed alternatives.

## Construction

Let `T_1 >= T_2 >= ...` be the observed motif frequencies. Fix, before
resampling:

- a target FDR `q` (`--alpha`);
- a resampling failure probability `beta < q`;
- a minimum useful discovery count `R0`;
- `K >= R0` motifs to mine (`--max-motifs K`).

Set

```
q0 = q - beta
k  = 1 + floor(q0 * R0).
```

For every null label permutation, PaM-ST mines the top `K` motif frequencies
and retains the `k`-th largest value. The rejection threshold is the maximum
of that order statistic over all resamples. Observed motifs are candidate
discoveries only when their frequency is strictly larger than this threshold.

In FDR mode, motif representatives are selected greedily in decreasing
frequency with pairwise distance strictly greater than `2 * rho`. Their
radius-`rho` match balls therefore cannot overlap. The same packing rule is
applied to the observed data and every null resample. This prevents many tiny
variations of one compositional mode from artificially increasing the number
of discoveries in the denominator of the FDP bound.

The automatic FewRS budget is

```
m = ceil(log(1 / beta) / log(1 / (1 - beta))).
```

This is independent of the number of motifs. For example, `beta=0.05` needs
59 resamples and `beta=0.025` needs 146.

## Guarantee

Under the same subset-pivotality assumption used by FewRS, calibrating the
`k`-th complete-null order statistic controls

```
P(V >= k) <= beta,
```

where `V` is the number of false candidate motifs. Let `R` be the number of
candidates. PaM-ST returns no motifs unless both `R >= R0` and

```
beta + (k - 1) / R <= q.
```

On the event `V < k`, every certified output has FDP at most
`(k - 1) / R <= q - beta`. On the complementary event the FDP is at most one,
and that event has probability at most `beta`. Therefore the final reporting
rule has FDR at most `q`. The value written as `fdr_bound` is this conservative
certification bound, not an estimate of the realized FDP. The guarantee does
not require independence among motif scores; dependence is handled by the
order-statistic exceedance event.

The case `k=1` reduces to a FewRS-style no-false-discovery event and is as
conservative as FWER control. FDR gains require enough anticipated discoveries
that `floor((q - beta) * R0) >= 1`.

## Assumptions and scope

- Label permutations must be independent draws from the intended null. The
  default global model shuffles throughout the tissue; `--null-model block`
  conditions on the cell-type counts in each fixed grid tile and shuffles
  within tiles. See [block null settings](block_null.md).
- With `--freeze-cell-type NAME`, the null is conditional on the locations of
  `NAME` cells: those labels stay fixed, all other labels are shuffled, and
  only non-`NAME` cells are candidate motif centers.
- Subset pivotality must hold for true-null motif statistics.
- Radius, distance metric, `rho`, `q`, `beta`, `R0`, and `K` must be fixed
  before inspecting the resampled thresholds. Searching across configurations
  requires including those configurations in the tested family or using an
  independent tuning split.
- The guarantee applies to the motif-frequency analyses returned by this run.
  It is not a statement about cell-type annotations or biological replication.
- The `2 * rho` packing rule makes the reported motifs non-overlapping in
  composition space. Changing that rule changes the tested family and must be
  mirrored in every null resample.

## Usage

```bash
./pam_st \
  --input cells.csv \
  --output-dir output/fewrs_fdr \
  --metric l2 \
  --radius 300 \
  --rho 0.05 \
  --error-control fewrs-fdr \
  --alpha 0.10 \
  --fdr-failure-probability 0.05 \
  --max-motifs 40 \
  --fdr-min-discoveries 40 \
  --permutations auto
```

With these settings, `k=3`, 59 label permutations are used, and a nonempty
result is certified only when `0.05 + 2/R <= 0.10`.

## Outputs

The existing `pattern.csv` and `matches.csv` continue to describe the first
ranked motif for compatibility. In `fewrs-fdr` mode the run also writes:

- `motif_tests.csv`: all tested ranks and their certification decisions;
- `motif_patterns.csv`: long-form compositions for every tested motif;
- `motif_matches.csv`: matches for every certified motif, keyed by rank;
- `fewrs_fdr_null_order.csv`: the calibrated null order statistic;
- `motif_separation` in `summary.csv` and `report.txt`: the minimum distance
  (`2 * rho`) enforced between reported motif centers;
- `summary.csv`: the target FDR, failure budget, order `k`, threshold, achieved
  bound, candidate and certified discovery counts, and subset-pivotality
  assumption.
