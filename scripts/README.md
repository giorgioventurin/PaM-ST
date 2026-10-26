# Plots of PaM-ST results

Each script reads one `pam_st --output-dir` (`summary.csv`, `patterns.csv`,
`labellings.csv`, `occurrences.csv`) and writes its figure to `plots/` inside
it. Run them from `code/`; `--help` lists every option.

| Script | Shows |
|---|---|
| `plot_occurrences.py RUN --rank 1` | every occurrence of a pattern, circled in red, in every sample, and a few up close |
| `plot_statistics.py RUN` | observed z-scores, the null of best scores, p-values, recurrence against frequency, patterns per labelling |
| `plot_runtime.py RUN [RUN ...]` | time by phase, and seconds per labelling against cells and universe size; a folder of runs works too |
| `plot_compositions.py RUN` | the mean cell-type shares of the top significant patterns, grouped by main cell type |
| `plot_hotspots.py RUN` | where the significant patterns concentrate, and of which kind |

The scripts that draw the tissue find the input tables from `summary.csv`. If
they have moved, give them again with `--input`, once per sample, in the order
given to `pam_st`.

`pam_results.py` holds what the scripts share: reading a run and the plot style.

```bash
python3 -m pip install -r scripts/requirements.txt
```
