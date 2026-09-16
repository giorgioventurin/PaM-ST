# Cumulative cell-type freeze sweep

## Purpose

The cumulative freeze sweep measures how motif significance changes after the
spatial locations of increasingly many abundant cell types are treated as a
fixed background. The cell-type order is computed once from the observed
counts, with ties resolved by cell-type name.

For stage `j`, the frozen set contains the `j` most abundant cell types. Stage
zero freezes no types. A frozen cell keeps its observed label and remains in
all neighborhood count vectors, but it cannot be a motif center. Labels are
shuffled only among non-frozen cell locations.

With `--null-model block --block-size L`, that shuffle is further restricted to
non-frozen locations within each grid tile. The same grid is used at every
stage. See [the block null model](block_null.md) for grid options and diagnostics.

The sweep stops while at least two cell types remain active. With only one
active type the conditional permutation distribution is degenerate.

## Stage-level null hypothesis

Let `F_j` be the frozen types and let

```text
I_j = {i : label(i) is not in F_j}
```

be the eligible motif centers. Conditional on coordinates, frozen labels and
locations, and the active label counts, the stage null is exchangeability of
the labels over `I_j`.

For the block model, exchangeability is within each tile's subset of `I_j`,
conditional also on the active label counts in each tile.

The primary statistic is the largest motif frequency among eligible centers.
This maximum includes the motif search in the permutation test. With `B`
conditional permutations, the raw stage p-value is

```text
p_j = (1 + number of null maxima >= observed maximum) / (B + 1).
```

The abundance order is valid under this conditional test because the label
counts are fixed by every permutation. The order must not be changed using
intermediate motif results or p-values.

## Whole-sweep error control

The stage hypotheses are dependent and are not logically nested. By default,
PaM-ST applies Holm correction to all completed stage p-values. For sorted raw
p-values `p_(1) <= ... <= p_(S)`, the adjusted values are

```text
p_adj_(r) = min(1, max over s <= r of (S - s + 1) * p_(s)).
```

Holm controls the probability of at least one false stage-level rejection at
the requested `alpha` under arbitrary dependence. Use
`--stage-error-control none` only for explicitly exploratory output.

The current `fewrs-fdr` mode controls motifs for one frozen configuration. It
is intentionally disallowed in cumulative-sweep mode because running it at
the same target independently in every stage does not provide pooled FDR
control across the sweep.

## Usage

```bash
./pam_st \
  --input cells.csv \
  --output-dir output/freeze_sweep \
  --metric l2 \
  --radius 300 \
  --rho 0.05 \
  --permutations 10000 \
  --max-motifs 3 \
  --freeze-by-abundance \
  --stage-error-control holm
```

To stop after freezing at most three types, add:

```text
--max-freeze-stages 3
```

The existing single-configuration option remains available and is repeatable:

```text
--freeze-cell-type Other_fibroblasts --freeze-cell-type Other_cancer
```

Manual frozen types and `--freeze-by-abundance` are mutually exclusive.

## Outputs

The sweep output root contains:

- `freeze_order.csv`: abundance order, counts, fractions, and included stages;
- `freeze_sweep_summary.csv`: rank-one statistics, raw and adjusted p-values,
  candidate fractions, effect sizes, patterns, and timings;
- `report.txt`: a compact text summary of every stage;
- `stage_00_none`, `stage_01_top_1`, and subsequent stage directories using
  the existing PaM-ST per-run CSV contract.

Every stage `summary.csv` records the cumulative frozen types, stage index,
active label count, raw stage p-value, adjusted stage p-value, and whole-sweep
significance decision. Ranked motifs beyond rank one are retained as
descriptive results; the Holm inference is based on the rank-one maximum.

## Plotting one dataset

From the repository's `code/` directory, plot the pattern and significance
trajectory for one sweep with:

```bash
python scripts/plot_freeze_sweep.py ../PaM-ST/output/path/to/freeze_sweep
```

The figure contains the mean raw count vector across all matched neighborhoods
at every stage, the raw and stage-adjusted p-values, and the z-score. Its
default destination is
`plots/freeze_sweep/<dataset>/freeze_sweep_pattern_significance.png`.

Use `--pattern-scale log1p` to make low-abundance components more visible, or
`--pattern-scale proportions` for column-normalized composition. Use
`--output-dir PATH` to select a different destination. Relative input and output
paths are resolved from `scripts/`, consistently with the other plotting
scripts.

## Interpretation

A loss of significance after freezing a type means that the remaining motif
is not unusual under a null conditioned on the frozen spatial scaffold. It
does not establish that the newly frozen type caused the original pattern.

Raw motif frequencies are not directly comparable between stages because the
candidate-center population shrinks. Prefer observed frequency divided by
candidate centers, observed-to-null ratio, z-score, and adjusted p-value.

## Computational reuse

The all-cell spatial neighbor graph and observed all-cell histograms are built
once. Each stage selects eligible rows from these cached structures. Only the
conditional label permutations and their histograms are recomputed per stage.
