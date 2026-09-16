# Block null model for the original PaM-ST

The original `pam_st` executable now supports shuffling labels independently
inside square tiles of a fixed grid. Enable it with `--null-model block`.
`--null-model global` remains the default and retains the original global
shuffle sequence for the same seed and build environment.

This option changes label resampling. The existing radius-based neighborhoods,
normalized-composition matching, candidate ranking, and statistical calculations
are used as before. Neighborhoods can still have variable numbers of cells and
contain only one cell type. The separate `mixed_motifs` implementation is not
used by these commands.

## Build

Run from the PaM-ST project directory:

```sh
cmake -S . -B build-block -DCMAKE_BUILD_TYPE=Release
cmake --build build-block -j 4
ctest --test-dir build-block --output-on-failure
```

The executable is `build-block/pam_st`. This separate build directory leaves
previously compiled executables available for comparison.

## Run one tissue

```sh
./build-block/pam_st \
  --input ../data/PDAC_202602/annotated_cells/LSP31895_P181_rep1_annotated.csv \
  --output-dir output/block_LSP31895_tile1000 \
  --radius 300 --rho 0.05 --metric l2 \
  --null-model block --block-size 1000 \
  --permutations 199 --max-motifs 3 --seed 37 --threads 4
```

The radius, tolerance, and permutation count here are examples; keep your
existing experiment settings when comparing null models. Tile size is measured
in the input coordinate units. The value 1000 is an initial configurable value,
not a calibrated biological distance.

To run the three datasets with the same settings:

```sh
for sample in LSP31894 LSP31895 LSP31896; do
  ./build-block/pam_st \
    --input "../data/PDAC_202602/annotated_cells/${sample}_P181_rep1_annotated.csv" \
    --output-dir "output/block_${sample}_tile1000" \
    --radius 300 --rho 0.05 --metric l2 \
    --null-model block --block-size 1000 \
    --permutations 199 --max-motifs 3 --seed 37 --threads 4 || break
done
```

Use a different output directory when changing settings; the original program
writes result files into the specified directory.

## Grid definition

| Option | Meaning | Default |
| --- | --- | --- |
| `--null-model global\|block` | Label resampling model | `global` |
| `--block-size L` | Positive, finite square tile side length | `1000` |
| `--block-origin-x X` | Grid x origin in input coordinate units | `0` |
| `--block-origin-y Y` | Grid y origin in input coordinate units | `0` |

Block size and origin options require `--null-model block`. A cell at `(x,y)`
belongs to tile

```text
tile_x = floor((x - origin_x) / block_size)
tile_y = floor((y - origin_y) / block_size)
```

Tiles include their lower boundary and exclude their upper boundary. With side
1000 and origin `(0,0)`, x=999 belongs to column 0, x=1000 to column 1, and
x=-1 to column -1. The grid is fixed throughout the run; it is not fitted to the
label distribution, shifted between permutations, or reanchored at each freeze
stage. Negative coordinates and translated grid origins are supported.

All retained cells are assigned to the grid. For every permutation, the labels
of non-frozen cells are randomly reordered separately in each tile. This
preserves each tile's count of every cell type exactly. Coordinates and
neighborhood membership remain fixed; neighborhoods may cross tile boundaries.
Cells not chosen as covering-neighborhood centers still participate in the
appropriate tile's label shuffle.

The existing input reader expects `X_centroid,Y_centroid,Cell_Type` columns and
excludes `Unclassified` cells. Run one tissue coordinate system per input file:
the original reader does not separate specimens by a sample identifier.

## Frozen cells and existing modes

Both manual `--freeze-cell-type NAME` and cumulative `--freeze-by-abundance`
are supported. Frozen cells keep their labels and positions, remain in
neighborhood count vectors, and follow the original center-exclusion rules.
Other labels shuffle only among the non-frozen locations of the same tile.
The grid stays identical across cumulative stages, while shuffle groups shrink
as additional types are frozen.

The block null works with `l2` and `js`, overlapping and covering neighborhoods,
and the existing pointwise and FewRS options. The preexisting restriction on
combining cumulative freeze sweeps with FewRS remains in effect. The existing
statistical guarantees and assumptions of each scoring mode are unchanged; a
new shuffle model alone does not establish motif-wise error control.

## Outputs and diagnostics

Existing motif, match, and null CSV formats are preserved. `report.txt` and
`summary.csv` additionally record the null model and the grid settings, along
with:

| Field | Meaning |
| --- | --- |
| `null_blocks` | Occupied grid tiles, including tiles with only frozen cells; zero for the global model |
| `null_shuffle_groups` | Tiles with at least one non-frozen cell; one group for a nonempty global shuffle |
| `null_mixed_groups` | Shuffle groups with at least two different non-frozen labels |
| `null_exchangeable_cells` | Non-frozen cells in those mixed groups, whose labels can actually change |

`permutable_cells` retains its original meaning: all non-frozen cells, including
cells in pure or singleton tiles. Those tiles are valid but produce no label
changes. A run with no exchangeable cells has an unchanged null tissue on every
permutation. The program reports this explicitly.

Each freeze stage records its own diagnostics. The sweep summary also includes
the null model, grid settings, and per-stage group counts.

## Interpretation

The block null asks whether the existing motif statistic is unusually large
when labels are exchangeable inside each tile, conditional on the observed
tile compositions and any frozen cells. Coarse spatial differences in cell-type
abundance are retained. Fine spatial organization within a tile is randomized.

Tile size therefore defines part of the hypothesis: small tiles can leave very
little randomization, while one tile containing the whole tissue reduces to the
global null. Prespecify a size and origin for a confirmatory comparison. A sweep
of sizes or offsets can be useful exploration, but choosing the smallest
p-value from that sweep requires accounting for the additional comparisons.
