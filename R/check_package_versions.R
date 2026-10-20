# =============================================================================
# Report the R / Bioconductor / package versions used in this project
#
# Run from the repository root:   source("R/check_package_versions.R")
#
# Output:
#   - a table in the console
#   - output/package_versions.csv     (for the Reporting Summary)
#   - output/sessionInfo.txt          (full session information)
#   - a ready-to-paste "Software" sentence for the Methods
# =============================================================================

packages <- list(
  "snRNA-seq" = c("Seurat", "SeuratObject", "DoubletFinder", "harmony", "clustree",
                  "SingleR", "celldex", "scuttle", "MAST"),
  "Bulk RNA-seq" = c("DESeq2", "SummarizedExperiment", "limma", "org.Hs.eg.db"),
  "Spatial transcriptomics" = c("zellkonverter", "SingleCellExperiment", "RANN"),
  "Enrichment" = c("enrichR"),
  "Epigenetic age" = c("cgageR"),
  "Visualisation and data handling" = c("ggplot2", "ggpubr", "EnhancedVolcano", "gghalves",
                                        "ggalluvial", "VennDiagram", "patchwork", "cowplot",
                                        "dplyr", "tidyr", "openxlsx", "tableone", "rmarkdown")
)

get_version <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) as.character(packageVersion(pkg)) else NA_character_
}

versions <- do.call(rbind, lapply(names(packages), function(group) {
  data.frame(Analysis = group, Package = packages[[group]],
             Version = vapply(packages[[group]], get_version, character(1)),
             row.names = NULL)
}))

r_version    <- paste(R.version$major, R.version$minor, sep = ".")
bioc_version <- if (requireNamespace("BiocManager", quietly = TRUE)) {
  as.character(BiocManager::version())
} else NA_character_

cat("\nR version:            ", r_version, "\n")
cat("Bioconductor version: ", bioc_version, "\n")
cat("Platform:             ", R.version$platform, "\n\n")
print(versions, row.names = FALSE)

missing <- versions$Package[is.na(versions$Version)]
if (length(missing)) {
  cat("\nNot installed in this R library:", paste(missing, collapse = ", "), "\n")
}

# ---- Save files ----
dir.create("output", showWarnings = FALSE)
write.csv(rbind(data.frame(Analysis = "Base", Package = c("R", "Bioconductor"),
                           Version = c(r_version, bioc_version)),
                versions),
          "output/package_versions.csv", row.names = FALSE)
for (p in versions$Package[!is.na(versions$Version)]) {
  suppressPackageStartupMessages(try(library(p, character.only = TRUE), silent = TRUE))
}
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")

# ---- Methods sentence ----
fmt <- function(pkgs) {
  v <- versions[versions$Package %in% pkgs & !is.na(versions$Version), ]
  paste0(v$Package, " (v", v$Version, ")", collapse = ", ")
}
cat("\n---- Methods: Software ----\n\n")
cat(sprintf(paste(
  "All analyses were performed in R (v%s) with Bioconductor (v%s).",
  "snRNA-seq data were analysed with %s.",
  "Differential expression of the endothelial cell stimulation experiment was performed with %s.",
  "Spatial transcriptomics data were analysed with %s.",
  "Functional enrichment analyses were performed with %s.",
  "Figures were generated with %s.\n"),
  r_version, bioc_version,
  fmt(packages[["snRNA-seq"]]), fmt(packages[["Bulk RNA-seq"]]),
  fmt(packages[["Spatial transcriptomics"]]), fmt(packages[["Enrichment"]]),
  fmt(packages[["Visualisation and data handling"]])))
cat("\nFiles written: output/package_versions.csv, output/sessionInfo.txt\n")
