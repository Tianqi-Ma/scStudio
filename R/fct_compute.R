#' Compute wrappers around Seurat / Bioconductor
#'
#' Each function performs one analysis step on a Seurat object and returns the
#' updated object (or a result). Kept separate from the UI so the science is
#' testable and the modules stay thin. All heavy dependencies are checked by the
#' caller via [require_pkgs()].
#'
#' Methods reflect current (2023-2025) best practice; see README for references.
#'
#' @name fct_compute
#' @keywords internal
NULL

# ---- QC ---------------------------------------------------------------------

#' Compute QC metrics (mitochondrial / ribosomal / hemoglobin percentages)
#' @param obj Seurat object.
#' @param species "human" or "mouse" (controls gene-name patterns).
#' @return Seurat object with `percent.mt`, `percent.ribo`, `percent.hb` in meta.
#' @keywords internal
qc_add_metrics <- function(obj, species = c("human", "mouse")) {
  species <- match.arg(species)
  mt   <- if (species == "human") "^MT-"   else "^mt-"
  ribo <- if (species == "human") "^RP[SL]" else "^Rp[sl]"
  hb   <- if (species == "human") "^HB[^P]" else "^Hb[^p]"
  obj[["percent.mt"]]   <- Seurat::PercentageFeatureSet(obj, pattern = mt)
  obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(obj, pattern = ribo)
  obj[["percent.hb"]]   <- Seurat::PercentageFeatureSet(obj, pattern = hb)
  # Dissociation/stress-gene percentage (borrowed from scCancer): flags cells
  # stressed during tissue dissociation. Uses genes present in the object.
  diss <- dissociation_genes(species)
  diss <- intersect(diss, rownames(obj))
  if (length(diss) > 0) {
    obj[["percent.diss"]] <- Seurat::PercentageFeatureSet(obj, features = diss)
  }
  obj
}

#' Dissociation/stress gene set (van den Brink et al.), human or mouse
#' @keywords internal
dissociation_genes <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  g <- c("FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "ATF3", "HSPA1A",
         "HSPA1B", "HSP90AB1", "HSPB1", "DNAJB1", "DNAJA1", "DUSP1", "IER2",
         "NR4A1", "PPP1R15A", "SOCS3", "ZFP36", "UBC", "HSPA8", "JUN")
  if (species == "mouse") g <- paste0(substr(g, 1, 1),
                                      tolower(substring(g, 2)))
  unique(g)
}

#' Adaptive (MAD-based) outlier flags from a metadata data.frame (pure)
#'
#' The testable core of the QC rule: flags cells more than `nmads`
#' median-absolute-deviations from the median on log1p(counts), log1p(genes),
#' and (upper only) percent.mt -- the modern alternative to fixed cutoffs.
#'
#' @param md data.frame with `nCount_RNA`, `nFeature_RNA`, optional `percent.mt`.
#' @param nmads_lib,nmads_mt MAD multipliers.
#' @return Logical vector: TRUE = keep, FALSE = flagged outlier.
#' @keywords internal
qc_mad_keep_from_meta <- function(md, nmads_lib = 5, nmads_mt = 3) {
  is_out <- function(x, nmads, type = "both", log = FALSE) {
    if (log) x <- log1p(x)
    med <- stats::median(x, na.rm = TRUE)
    dev <- stats::mad(x, center = med, na.rm = TRUE)
    lower <- med - nmads * dev
    upper <- med + nmads * dev
    switch(type,
           both   = x < lower | x > upper,
           higher = x > upper,
           lower  = x < lower)
  }
  out_counts <- is_out(md$nCount_RNA,   nmads_lib, "both",   log = TRUE)
  out_genes  <- is_out(md$nFeature_RNA, nmads_lib, "both",   log = TRUE)
  out_mt     <- if (!is.null(md$percent.mt)) is_out(md$percent.mt, nmads_mt, "higher") else FALSE
  !(out_counts | out_genes | out_mt)
}

#' Adaptive (MAD-based) outlier flags for a Seurat object
#' @param obj Seurat object with QC metrics.
#' @param nmads_lib,nmads_mt MAD multipliers.
#' @return Logical keep vector.
#' @keywords internal
qc_mad_keep <- function(obj, nmads_lib = 5, nmads_mt = 3) {
  qc_mad_keep_from_meta(obj_meta(obj), nmads_lib, nmads_mt)
}

#' Apply manual QC thresholds to a metadata data.frame (pure)
#' @param md data.frame with `nFeature_RNA`, optional `percent.mt`.
#' @keywords internal
qc_manual_keep_from_meta <- function(md, min_genes, max_genes, max_mt) {
  keep <- md$nFeature_RNA >= min_genes & md$nFeature_RNA <= max_genes
  if (!is.null(md$percent.mt)) keep <- keep & md$percent.mt <= max_mt
  keep
}

#' Apply manual QC thresholds to a Seurat object
#' @keywords internal
qc_manual_keep <- function(obj, min_genes, max_genes, max_mt) {
  qc_manual_keep_from_meta(obj_meta(obj), min_genes, max_genes, max_mt)
}

# ---- Doublets ---------------------------------------------------------------

#' Run doublet detection (scDblFinder by default) and add a call to metadata
#' @param obj Seurat object; @param method "scDblFinder" or "DoubletFinder".
#' @return Seurat object with `doublet_score` and `doublet_class` in meta.
#' @keywords internal
run_doublets <- function(obj, method = c("scDblFinder", "DoubletFinder")) {
  method <- match.arg(method)
  if (method == "scDblFinder") {
    sce <- Seurat::as.SingleCellExperiment(obj)
    sce <- scDblFinder::scDblFinder(sce)
    obj[["doublet_score"]] <- SummarizedExperiment::colData(sce)$scDblFinder.score
    obj[["doublet_class"]] <- SummarizedExperiment::colData(sce)$scDblFinder.class
  } else {
    stop("DoubletFinder path requires per-dataset pK sweep; use scDblFinder for the one-click default.")
  }
  obj
}

# ---- Normalization ----------------------------------------------------------

#' Normalize counts
#' @param obj Seurat object; @param method "LogNormalize","SCT","scran".
#' @keywords internal
normalize_obj <- function(obj, method = c("LogNormalize", "SCT"), scale_factor = 1e4) {
  method <- match.arg(method)
  if (method == "LogNormalize") {
    obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize",
                                 scale.factor = scale_factor, verbose = FALSE)
  } else {
    obj <- Seurat::SCTransform(obj, verbose = FALSE)
  }
  obj
}

# ---- Features + PCA ---------------------------------------------------------

#' HVG selection, scaling, and PCA
#' @keywords internal
reduce_obj <- function(obj, n_hvg = 2000, npcs = 50, hvg_method = "vst") {
  obj <- Seurat::FindVariableFeatures(obj, selection.method = hvg_method,
                                      nfeatures = n_hvg, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = npcs, verbose = FALSE)
  obj
}

# ---- Integration ------------------------------------------------------------

#' Batch integration (Harmony by default; CCA/RPCA via Seurat v5)
#'
#' @param obj Seurat object (normalized, with PCA for harmony).
#' @param batch Metadata column identifying batch/sample.
#' @param method "harmony" (default), "none", "CCA", or "RPCA".
#' @param dims Dimensions to use for anchor-based methods.
#' @return Seurat object with an integrated reduction ("harmony" or
#'   "integrated.dr"); downstream steps can point their `reduction` at it.
#' @keywords internal
integrate_obj <- function(obj, batch, method = c("harmony", "none", "CCA", "RPCA"),
                          dims = 30) {
  method <- match.arg(method)
  if (method == "none") return(obj)

  if (method == "harmony") {
    obj <- harmony::RunHarmony(obj, group.by.vars = batch)
    return(obj)
  }

  # CCA / RPCA use Seurat v5 IntegrateLayers. Split layers by batch first.
  if (!"IntegrateLayers" %in% getNamespaceExports("Seurat")) {
    stop("CCA/RPCA integration requires Seurat v5 (IntegrateLayers). ",
         "Use Harmony, or upgrade Seurat.")
  }
  obj[["RNA"]] <- base::split(obj[["RNA"]], f = obj_meta(obj)[[batch]])
  obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, verbose = FALSE)
  m <- if (method == "CCA") "CCAIntegration" else "RPCAIntegration"
  obj <- Seurat::IntegrateLayers(obj, method = get(m, envir = asNamespace("Seurat")),
                                 orig.reduction = "pca",
                                 new.reduction = "integrated.dr", verbose = FALSE)
  obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
  obj
}

# ---- Clustering -------------------------------------------------------------

#' Neighbors + clustering at one or more resolutions
#' @keywords internal
cluster_obj <- function(obj, reduction = "pca", dims = 30,
                        resolutions = 0.5, algorithm = 4) {
  obj <- Seurat::FindNeighbors(obj, reduction = reduction, dims = seq_len(dims),
                               verbose = FALSE)
  for (res in resolutions) {
    obj <- Seurat::FindClusters(obj, resolution = res, algorithm = algorithm,
                                verbose = FALSE)
  }
  obj
}

# ---- Embedding --------------------------------------------------------------

#' Non-linear embedding for visualization (UMAP / t-SNE)
#' @keywords internal
embed_obj <- function(obj, method = c("umap", "tsne"), reduction = "pca",
                      dims = 30, n_neighbors = 30, min_dist = 0.3, perplexity = 30) {
  method <- match.arg(method)
  d <- seq_len(dims)
  if (method == "umap") {
    obj <- Seurat::RunUMAP(obj, reduction = reduction, dims = d,
                           n.neighbors = n_neighbors, min.dist = min_dist,
                           verbose = FALSE)
  } else {
    obj <- Seurat::RunTSNE(obj, reduction = reduction, dims = d,
                           perplexity = perplexity)
  }
  obj
}

# ---- Markers ----------------------------------------------------------------

#' Marker genes per cluster
#' @keywords internal
markers_obj <- function(obj, test = "wilcox", logfc = 0.25, min_pct = 0.1,
                        only_pos = TRUE) {
  Seurat::FindAllMarkers(obj, test.use = test, logfc.threshold = logfc,
                         min.pct = min_pct, only.pos = only_pos, verbose = FALSE)
}

# ---- Annotation -------------------------------------------------------------

#' Reference-based annotation with SingleR
#' @keywords internal
annotate_singler <- function(obj, ref) {
  sce <- Seurat::as.SingleCellExperiment(obj)
  pred <- SingleR::SingleR(test = sce, ref = ref$data, labels = ref$labels)
  obj[["SingleR"]] <- pred$labels
  obj
}

# ============================================================================
# scop engine wrappers.
# Thin wrappers over scop's Run* API so the modules stay declarative. All are
# gated by require_pkgs("scop") and wrapped in tryCatch. scop signatures are
# verified at runtime on the user's machine; a few arg names may need tweaks
# against the installed scop version (flagged in CHANGELOG/README).
# ============================================================================

#' Standard preprocessing (normalize + HVG + scale + linear reduction)
#' @param method normalization: "LogNormalize"/"SCT"/"TFIDF"; hvf: vst/mvp/disp.
#' @keywords internal
sc_standard <- function(srt, normalization = "LogNormalize", hvf = "vst",
                        nHVF = 2000, linear = "pca", npcs = 50) {
  scop::Standard_SCP(srt, normalization_method = normalization,
                     HVF_method = hvf, nHVF = nHVF,
                     linear_reduction = linear, linear_reduction_dims = npcs)
}

#' Linear dimensionality reduction (PCA/ICA/NMF/MDS)
#' @keywords internal
sc_reduce <- function(srt, method = "pca", ndim = 50) {
  scop::RunDimReduction(srt, prefix = method, reduction_method = method,
                        ndim = ndim)
}

#' Batch integration via scop::Integration_SCP (many methods)
#' @param method Uncorrected/Seurat/Harmony/scVI/scanorama/BBKNN/fastMNN/LIGER/CSS/Conos/ComBat
#' @keywords internal
sc_integrate <- function(srt, batch, method = "Harmony", nHVF = 2000) {
  scop::Integration_SCP(srt, batch = batch, integration_method = method,
                        nHVF = nHVF)
}

#' Non-linear embedding (UMAP/tSNE/DM/PHATE/FR/PaCMAP)
#' @keywords internal
sc_embed <- function(srt, method = "umap", reduction = "pca", dims = 30,
                     n_neighbors = 30, min_dist = 0.3) {
  scop::RunDimReduction(srt, prefix = method, reduction_method = method,
                        reduction_use = reduction, dims_use = seq_len(dims),
                        n.neighbors = n_neighbors, min.dist = min_dist)
}

#' Differential expression / markers via scop::RunDEtest
#' @keywords internal
sc_detest <- function(srt, group_by, test = "wilcox") {
  scop::RunDEtest(srt, group_by = group_by, test.use = test)
}

#' Reference/label transfer annotation via scop::RunKNNPredict
#' @keywords internal
sc_annotate_knn <- function(srt, ref, ref_group, distance = "cosine") {
  scop::RunKNNPredict(srt_query = srt, srt_ref = ref, ref_group = ref_group,
                      distance_metric = distance)
}

#' Over-representation enrichment via scop::RunEnrichment
#' @keywords internal
sc_enrichment <- function(srt, group_by, db = "GO", species = "Homo_sapiens") {
  scop::RunEnrichment(srt, group_by = group_by, db = db, species = species)
}

#' GSEA via scop::RunGSEA
#' @keywords internal
sc_gsea <- function(srt, group_by, db = "GO", species = "Homo_sapiens") {
  scop::RunGSEA(srt, group_by = group_by, db = db, species = species)
}

#' Trajectory / pseudotime (Slingshot/Monocle2/Monocle3/PAGA/Palantir/WOT)
#' @keywords internal
sc_trajectory <- function(srt, method = "slingshot", group_by = NULL, ...) {
  fun <- switch(tolower(method),
                slingshot = scop::RunSlingshot,
                monocle2  = scop::RunMonocle2,
                monocle3  = scop::RunMonocle3,
                paga      = scop::RunPAGA,
                palantir  = scop::RunPalantir,
                wot       = scop::RunWOT,
                stop("Unknown trajectory method: ", method))
  fun(srt, group.by = group_by, ...)
}

#' RNA velocity via scop::RunSCVELO (Python)
#' @keywords internal
sc_velocity <- function(srt, group_by = NULL, mode = "dynamical", ...) {
  scop::RunSCVELO(srt, group_by = group_by, mode = mode, ...)
}

#' Dynamic features along pseudotime via scop::RunDynamicFeatures
#' @keywords internal
sc_dynamic <- function(srt, lineages, n_candidates = 200) {
  scop::RunDynamicFeatures(srt, lineages = lineages, n_candidates = n_candidates)
}

#' Cell-cycle scoring (Seurat) + module scores (UCell/AddModuleScore)
#' @keywords internal
sc_cellcycle <- function(srt) {
  cc <- Seurat::cc.genes.updated.2019
  Seurat::CellCycleScoring(srt, s.features = cc$s.genes,
                           g2m.features = cc$g2m.genes, set.ident = FALSE)
}

#' Signature module scoring; method "UCell" or "AddModuleScore"
#' @param features Named list of gene sets.
#' @keywords internal
sc_modulescore <- function(srt, features, method = "UCell") {
  if (method == "UCell" && has_pkg("UCell")) {
    UCell::AddModuleScore_UCell(srt, features = features)
  } else {
    Seurat::AddModuleScore(srt, features = features, name = names(features))
  }
}

#' Cell-cell communication (LIANA/CellChat/CellPhoneDB/NicheNet)
#' @keywords internal
sc_cellcomm <- function(srt, group_by, method = "liana", ...) {
  method <- tolower(method)
  if (method == "liana") {
    if (!require_pkgs("liana", "Cell-cell communication")) return(NULL)
    liana::liana_wrap(srt, idents_col = group_by, ...)
  } else if (method == "cellchat") {
    if (!require_pkgs("CellChat", "Cell-cell communication")) return(NULL)
    cc <- CellChat::createCellChat(object = srt, group.by = group_by)
    cc
  } else {
    stop("Cell-cell communication method '", method,
         "' runs via Python/other setup; see docs.")
  }
}

#' Malignant-cell / CNV estimation (CopyKAT / inferCNV / Numbat)
#' @keywords internal
sc_cnv <- function(srt, method = "copykat", ref_cells = NULL, ...) {
  method <- tolower(method)
  if (method == "copykat") {
    if (!require_pkgs("copykat", "CNV / malignant cells")) return(NULL)
    mat <- as.matrix(SeuratObject::LayerData(srt, layer = "counts"))
    norm_names <- if (!is.null(ref_cells)) ref_cells else ""
    copykat::copykat(rawmat = mat, norm.cell.names = norm_names, ...)
  } else {
    stop("CNV method '", method, "' requires its own setup (inferCNVpy/Numbat).")
  }
}

#' Stemness score (mRNAsi-style signature)
#' @keywords internal
sc_stemness <- function(srt, features = NULL) {
  if (is.null(features)) {
    stop("Provide a stemness gene signature (features).")
  }
  sc_modulescore(srt, list(stemness = features))
}
