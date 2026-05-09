# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository overview

GeneTrajectory is an R package (R >= 4.0) for inferring gene trajectories in single-cell RNA-seq data using graph-based optimal transport. It integrates with Seurat objects and delegates the most expensive step — pairwise gene-gene Wasserstein distance computation — to a companion Python package (`gene-trajectory`) called via `reticulate`.

A standalone Python implementation lives at https://github.com/KlugerLab/GeneTrajectory-python.

## Common commands

This is a roxygen2-documented R package with no test suite. Typical workflow from the package root:

```r
# Regenerate man/*.Rd files and NAMESPACE from roxygen comments in R/*.R.
# Required after any change to roxygen comments or @export tags.
devtools::document()

# Build the pkgdown site under docs/
pkgdown::build_site()

# Build vignettes (uses knitr; vignettes assume access to a Seurat .rds
# downloaded from figshare and are not self-contained).
devtools::build_vignettes()

# Install locally from source
devtools::install()

# Full CRAN-style check
devtools::check()              # or: R CMD check .
```

There are no unit tests (no `tests/` directory) and no linter config. Verify changes by running the example scripts in `scripts/` against a downloaded Seurat object, or by reproducing a vignette.

## Architecture: the four-step pipeline

The whole package is organized around one pipeline. Understanding the data flow between functions is more useful than reading files individually.

**Step 1 — Cell graph + graph distances (R, on Seurat object).**
`RunDM(seurat)` adds a Diffusion Map embedding to `seurat@reductions$dm`, computed from a precomputed pairwise distance matrix on PCA (or any other reduction) via `diffusion.map()`. `GetGraphDistance(seurat)` then builds a kNN graph on the DM embedding (via `FNN::get.knn`), turns it into an `igraph` graph, and returns the cell-cell shortest-path distance matrix. If the kNN graph has disconnected components, this errors out and asks the user to increase `K`.

**Step 2 — Coarse-grain, then compute gene-gene Wasserstein distances (R → Python).**
`CoarseGrain(seurat, graph.dist, features, N)` runs `kmeans` on the DM embedding to collapse cells into `N` meta-cells, returning a list with `gene.expression` (meta-cells × genes), `graph.dist` (meta-cell × meta-cell), and `features`. The Wasserstein computation itself is **not in R**; it runs in Python via `reticulate::import('gene_trajectory.compute_gene_distance_cmd')$cal_ot_mat_from_numpy(...)`. The R package can also write the coarse-grained matrices to CSV/MTX files (see `scripts/human_myeloid_example.R`) for an offline Python run. `LoadGeneDistMat(dir.path)` reads back the resulting `emd.csv` plus `gene_names.csv` into a labeled gene-gene distance matrix.

**Step 3 — Gene embedding + sequential trajectory extraction (R).**
`GetGeneEmbedding(dist.mat, K)` runs `diffusion.map` on the gene-gene Wasserstein matrix and drops the trivial first eigenvector. `ExtractGeneTrajectory(gene.embedding, dist.mat, N, t.list, K)` then iterates `N` times: pick the gene farthest from the origin in the DM embedding as a terminus, run a `t`-step random walk on the gene graph (via `get.RW.matrix`) seeded at that gene, threshold by `quantile * max(diffused)`, label those genes as "Trajectory-i", and exclude them from subsequent rounds. Early stops if no unassigned genes remain.

**Step 4 — Pseudo-order genes within each trajectory and score cells (R).**
`GetGenePseudoorder(dist.mat, subset, max.id)` runs Diffusion Map on the within-trajectory submatrix and ranks genes by the second eigenvector (flipping sign so the chosen terminus sits at the high end). `ExtractGeneTrajectory` already calls this and stores the result as `Pseudoorder1`, `Pseudoorder2`, ... columns in the returned data frame. Finally, `AddGeneBinScore(seurat, gene.trajectory, N.bin, trajectories, assay, reverse)` splits each trajectory into `N.bin` ordered gene bins and writes per-cell bin scores back into `seurat@meta.data` as `Trajectory{i}_genes{1..N.bin}`.

## Naming conventions

Two distinct API layers, distinguished by case:

- **`PascalCase` functions** (`RunDM`, `GetGraphDistance`, `CoarseGrain`, `GetGeneEmbedding`, `ExtractGeneTrajectory`, `GetGenePseudoorder`, `AddGeneBinScore`, `LoadGeneDistMat`) are the user-facing API. They take/return Seurat objects or labeled matrices.
- **`dot.case` functions** (`diffusion.map`, `coarse.grain`, `get.RW.matrix`) are lower-level primitives operating on plain matrices. `CoarseGrain` is a thin Seurat wrapper over `coarse.grain`; `RunDM` and `GetGeneEmbedding` are wrappers over `diffusion.map`.

Both layers are exported (see `NAMESPACE`). When adding a new pipeline step, follow the same split if the underlying logic is matrix-only and may be reused.

## Parameter consistency rules

Several parameters must match across calls or results will be silently inconsistent:

- `K` in `GetGeneEmbedding` and `K` in `ExtractGeneTrajectory` must be the same — both build adaptive kernels on the same gene-gene distance matrix.
- `dims` and `reduction` in `GetGraphDistance` and `CoarseGrain` should refer to the same cell embedding used in `RunDM`.
- `N` in `CoarseGrain` (number of meta-cells) trades accuracy for speed; see `vignettes/fast_computation.Rmd` for the documented strategies (cell-graph coarse-graining and gene-graph sparsification via `FNN::get.knn` on a low-N preliminary distance matrix).

## Conventions for editing R source

- Every exported function has roxygen2 comments above it. Editing the signature or `@export`/`@param` tags requires running `devtools::document()` so `man/*.Rd` and `NAMESPACE` stay in sync — both are committed and machine-generated.
- `RoxygenNote: 7.2.3` in `DESCRIPTION` pins the roxygen2 version used to generate the docs; bumping it produces large diffs.
- The `simulate_*` functions (`simulate_linear`, `simulate_bifurcation`, `simulate_cyclic`, `simulate_cylinder`, `simulate_coral`) generate synthetic gene-by-cell count matrices for testing/benchmarking the pipeline. They use `MASS::rnegbin` for the negbin model — note `MASS` is not in `Imports`, so the negbin path will fail unless `MASS` is loaded.
- `data/test/` contains a small fixture (`gene_expression.mtx`, `ot_cost.csv`, `gene_names.csv`, `emd.csv`) representing a coarse-grained input + its expected Wasserstein output; useful as a smoke-test fixture for `LoadGeneDistMat` and downstream steps without needing Python.

## Python interop

The Python `gene-trajectory` package is required at runtime for Step 2. The recommended setup (per `README.md`) is:

```r
if (!reticulate::virtualenv_exists('gene_trajectory')) {
  reticulate::virtualenv_create('gene_trajectory', packages = c('gene_trajectory'))
}
reticulate::use_virtualenv('gene_trajectory')
```

When calling `cal_ot_mat_from_numpy` with a precomputed `gene_pairs` matrix (sparsification path), subtract `1L` from R indices to convert to Python 0-based indexing — this is easy to forget and silently produces wrong distances.
