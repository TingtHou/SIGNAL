# Creates data/SIGNAL_example.rda. Run from the package root with SIGNAL installed.
library(SIGNAL)

SIGNAL_example <- simulate_signal_data(
  n_cells = 400,
  n_genes = c(null = 40, celltype = 10, spatial = 10, ctSVG = 10),
  n_celltypes = 4, fold = 4, theta = 10, seed = 2025
)
storage.mode(SIGNAL_example$counts) <- "integer"
save(SIGNAL_example, file = "data/SIGNAL_example.rda", compress = "xz")
