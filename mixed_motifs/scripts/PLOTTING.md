# Plotting result tables

Run the C++ analysis first, then generate the complete figure set:

```sh
python -m pip install -r mixed_motifs/requirements.txt
python mixed_motifs/scripts/plot_results.py --results path/to/results --top 8
```

`--output` chooses a different output directory (default: `RESULTS/plots`).
`--atlas-matches 3` selects the number of actual occurrences shown per motif;
`--dpi 300` raises PNG resolution. Every figure is also saved as a PDF with
editable text and rasterized dense tissue point clouds.

The script reads exported tables without recomputing neighborhoods or changing
the analysis. It preserves the C++ motif reporting order and shares one cell-type
palette across all figures. Labels in legends are wrapped rather than truncated.

| Output | What to inspect |
| --- | --- |
| `overview` | Template counts, matching centers versus cell-disjoint support, support against same-template null realizations, and statistical evidence. |
| `neighborhood_diagnostics` | Physical radii, eligibility, diversity, dominance, and cross-type fractions among member pairs within `contact_radius`. The edge panel is explicitly disabled when that radius is zero. |
| `sample_recurrence` | Counts in each specimen and descriptive rates per 1,000 input cells. |
| `spatial_XX_SAMPLE` | One specimen per figure, with a cell-type reference and maps of actual cells belonging to selected occurrences. An outlined point marks each occurrence center. |
| `occurrence_atlas_XX` | Actual member cells for representative groups, covering different specimens first and then separated positions. Coordinates are translated to the group center; all panels on one page retain the same scale and correct aspect ratio. |
| `cell_type_palette.csv` | Reusable mapping from full cell-type labels to hexadecimal colors. |
| `plots_manifest.json`, `figure_captions.md` | File list, figure descriptions, and interpretation notes. |

Discovery and validation are deliberately distinguished. In discovery, the
null histogram compares the maximum observed support with maxima from complete
null searches. Same-template null means and ranges are descriptive after
selection and are not motif-specific significance tests. In validation, the
statistical panel shows Holm-adjusted fixed-template p-values; the templates
and settings must have been chosen independently of the validation tissue.

Cell-disjoint occurrences are not necessarily independent biological
replicates. Counts from distinct motifs can include the same cells. Recurrence
rates per 1,000 cells do not adjust for cell-type abundance or compartment
composition. Coordinate distances retain the input data's units. The atlas
shows membership and positions, without inventing cell boundaries or biological
contacts.

For a run with no reportable motifs, the script still generates an overview,
geometric diagnostics, and tissue reference maps with explicit empty-result
annotations.

Compare completed discovery runs with identical settings:

```sh
python mixed_motifs/scripts/compare_runs.py \
  --results mixed_motifs/results/LSP31891_blocks \
            mixed_motifs/results/LSP31894_blocks \
            mixed_motifs/results/LSP31895_blocks \
            mixed_motifs/results/LSP31896_blocks \
  --output mixed_motifs/results/comparison
```

This writes a CSV/JSON summary and `null_comparison.png` / `.pdf`. It checks
analysis settings, permutation completeness, and recorded global p-values
against the exported null maxima before plotting. Holm adjustment covers all
supplied dataset-level tests; it is not motif-level error control. Include the
full dataset family you intend to compare, rather than only significant runs.
