#' Plain-language cell-type explanation dictionary
#'
#' A lookup used by the annotation module's hover tooltips: given a predicted
#' cell-type label, show a one-line, beginner-friendly explanation of what that
#' cell type is. Matching is case-insensitive and tolerant of common synonyms.
#'
#' This is intentionally a small, general starter set covering frequent immune,
#' stromal, epithelial, and developmental types. Extend `celltype_dictionary()`
#' as needed for your tissue of interest.
#'
#' @name data_celltype_dict
#' @keywords internal
NULL

#' The cell-type explanation table
#'
#' @return A data.frame with columns `key` (lowercase canonical token) and
#'   `explanation` (one-line plain-language description).
#' @keywords internal
celltype_dictionary <- function() {
  d <- c(
    "t cell"            = "Immune cell that coordinates and executes adaptive immune responses.",
    "cd4 t cell"        = "Helper T cell; directs other immune cells (a coordinator).",
    "cd8 t cell"        = "Cytotoxic T cell; kills infected or abnormal cells.",
    "regulatory t cell" = "Treg; dampens immune responses to prevent over-reaction.",
    "b cell"            = "Immune cell that produces antibodies.",
    "plasma cell"       = "Mature B cell specialised for mass antibody production.",
    "nk cell"           = "Natural killer cell; kills stressed/infected cells without prior priming.",
    "monocyte"          = "Circulating immune cell that becomes a macrophage in tissue.",
    "macrophage"        = "Tissue immune cell that engulfs debris and pathogens.",
    "dendritic cell"    = "Antigen-presenting cell that activates T cells.",
    "neutrophil"        = "Fast-responding immune cell against bacterial infection.",
    "mast cell"         = "Immune cell releasing histamine; involved in allergy.",
    "erythrocyte"       = "Red blood cell; carries oxygen (often a contaminant in scRNA-seq).",
    "platelet"          = "Cell fragment involved in blood clotting.",
    "epithelial cell"   = "Cell forming the lining of organs, glands and surfaces.",
    "endothelial cell"  = "Cell lining blood and lymphatic vessels.",
    "fibroblast"        = "Structural cell producing extracellular matrix (connective tissue).",
    "smooth muscle cell"= "Involuntary muscle cell in vessels and organs.",
    "pericyte"          = "Cell wrapping capillaries to regulate blood flow.",
    "cardiomyocyte"     = "Heart muscle cell responsible for contraction.",
    "neuron"            = "Nerve cell that transmits electrical/chemical signals.",
    "astrocyte"         = "Support (glial) cell in the brain maintaining neurons.",
    "oligodendrocyte"   = "Glial cell that myelinates neurons in the CNS.",
    "microglia"         = "Resident immune cell of the brain.",
    "hepatocyte"        = "Main functional cell of the liver.",
    "stem cell"         = "Undifferentiated cell that can self-renew and specialise.",
    "progenitor cell"   = "Partially committed cell on its way to a mature type.",
    "proliferating cell"= "Actively dividing cell (high cell-cycle gene expression)."
  )
  data.frame(
    key = names(d),
    explanation = unname(d),
    stringsAsFactors = FALSE
  )
}

#' Look up a plain-language explanation for a cell-type label
#'
#' @param label Character vector of predicted cell-type labels.
#' @return Character vector of explanations (empty string if not found).
#' @keywords internal
explain_celltype <- function(label) {
  dict <- celltype_dictionary()
  norm <- tolower(trimws(as.character(label)))
  # strip common suffixes/markers so "CD8+ T cells" matches "cd8 t cell"
  norm <- gsub("\\+", "", norm)
  norm <- gsub("s$", "", norm)          # naive plural -> singular
  norm <- gsub("_", " ", norm)
  out <- vapply(norm, function(x) {
    hit <- which(dict$key == x)
    if (length(hit)) return(dict$explanation[hit[1]])
    # partial contains match as fallback
    part <- which(vapply(dict$key, function(k) grepl(k, x, fixed = TRUE), logical(1)))
    if (length(part)) return(dict$explanation[part[1]])
    ""
  }, character(1))
  unname(out)
}
