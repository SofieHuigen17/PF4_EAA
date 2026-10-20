# =============================================================================
# One-off export of the manual lumen annotation used in the manuscript.
#
# The lumen boundaries were traced interactively (see 04_spatial_Xenium.Rmd,
# chunk `draw-lumen`). This script exports the saved boundaries and the
# per-sample distance cut-offs to plain files in this folder, so that the
# luminal-adjacent EC classification can be reproduced without redrawing.
#
# Input : lumen_boundaries.rds            (list of boundary coordinates per sample)
#         panel1_with_manual_lumen.rds    (Xenium SCE with column lumen_cutoff_used)
# Output: lumen_annotation/lumen_boundaries.rds
#         lumen_annotation/lumen_boundaries.csv
#         lumen_annotation/lumen_cutoffs.csv
# =============================================================================

library(SingleCellExperiment)

boundary_file <- "path/to/lumen_boundaries.rds"
panel_file    <- "path/to/panel1_with_manual_lumen.rds"
out_dir       <- "lumen_annotation"

lumen_boundaries <- readRDS(boundary_file)
panel1 <- readRDS(panel_file)

cutoffs <- unique(data.frame(sample = as.character(panel1$lumen_sample),
                             cutoff = panel1$lumen_cutoff_used))
cutoffs <- cutoffs[!is.na(cutoffs$sample), ]
stopifnot(!anyDuplicated(cutoffs$sample))

saveRDS(lumen_boundaries, file.path(out_dir, "lumen_boundaries.rds"))
write.csv(cutoffs, file.path(out_dir, "lumen_cutoffs.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(names(lumen_boundaries), function(s) {
  data.frame(sample = s, point = seq_len(nrow(lumen_boundaries[[s]])),
             x_um = lumen_boundaries[[s]][, 1], y_um = lumen_boundaries[[s]][, 2])
})), file.path(out_dir, "lumen_boundaries.csv"), row.names = FALSE)
