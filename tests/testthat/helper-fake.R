# A stand-in "working object" for tests.
#
# The real object is a Seurat object, which cannot be built here: Seurat lives in
# Suggests and is absent on CI / on a bare dev box. A named counts matrix is
# enough for what these tests assert -- that every module's outputs *evaluate*.
# `obj_dims()` reads its dims, `obj_meta()` / `obj_reductions()` fall through
# their tryCatch guards to empty results, and any module that actually needs
# Seurat stops at its `require_pkgs()` gate, which is exactly the degraded path
# we want covered.
fake_obj <- function(n_genes = 40, n_cells = 20) {
  m <- matrix(stats::rpois(n_genes * n_cells, lambda = 3),
              nrow = n_genes, ncol = n_cells)
  dimnames(m) <- list(paste0("GENE", seq_len(n_genes)),
                      paste0("Cell", seq_len(n_cells)))
  m
}
