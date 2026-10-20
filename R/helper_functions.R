# =============================================================================
# Shared helper functions
# PF4 / epigenetic age acceleration (EAA) in human carotid plaques
#
# Sourced by all four analysis notebooks (01-04).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Colour schemes used throughout the manuscript figures
# -----------------------------------------------------------------------------
col_EAA <- c("Decelerated" = "#1b9e77", "Accelerated" = "#8152A2")
col_condition <- c("Control" = "#E69F00", "PF4" = "#004C99", "TGF" = "#CC79A7")
col_location <- c("Lumen-adjacent ECs" = "#DD8452", "Non-lumen-adjacent ECs" = "#4C72B0")

enrichment_palettes <- list(
  ACC  = c("#C9A9D8", "#8152A2", "#5B3775"),
  pf4  = c("#D6EAF8", "#5B9BD5", "#004C99", "#00264D"),
  tgf  = c("#E6F5FC", "#91D3F2", "#56B4E9", "#0077B6"),
  purp = c("#DCCBE8", "#7E4D9E", "#3F274F")
)

# -----------------------------------------------------------------------------
# Output helpers
# -----------------------------------------------------------------------------
save_fig <- function(plot, filename, width, height, out_dir, units = "cm", dpi = 300) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(out_dir, filename), plot = plot,
         width = width, height = height, units = units, dpi = dpi)
  invisible(plot)
}

# -----------------------------------------------------------------------------
# Enrichr
# -----------------------------------------------------------------------------
enrichr_dbs_GO <- c("GO_Molecular_Function_2025",
                    "GO_Cellular_Component_2025",
                    "GO_Biological_Process_2025")
enrichr_dbs_pathways <- c("Reactome_Pathways_2024",
                          "KEGG_2021_Human",
                          "MSigDB_Hallmark_2020")

init_enrichr <- function() {
  enrichR::setEnrichrSite("Enrichr")
  invisible(TRUE)
}

#' Run Enrichr against GO and pathway libraries with a custom background
run_enrichr <- function(genes, background) {
  genes <- unique(na.omit(genes))
  background <- unique(na.omit(background))
  c(enrichR::enrichr(genes, enrichr_dbs_GO, background = background),
    enrichR::enrichr(genes, enrichr_dbs_pathways, background = background))
}

#' Horizontal bar plot of selected Enrichr terms (gene count on x, colour =
#' -log10 adjusted P). Term labels are drawn inside (white) or next to the bar.
#'
#' @param enrich_table  Enrichr result table (already subset to the terms shown)
#' @param inside_n      number of top terms whose label is drawn inside the bar
#' @param inside_terms  alternatively, the exact terms to label inside the bar
plot_enrichment <- function(enrich_table, title = "", n_terms = 12, num_char = 50,
                            palette = "ACC", inside_n = 3, inside_terms = NULL) {
  df <- enrich_table %>%
    as_tibble() %>%
    arrange(Adjusted.P.value) %>%
    distinct(Term, .keep_all = TRUE) %>%
    slice_head(n = n_terms) %>%
    mutate(Term_short = str_trunc(Term, width = num_char),
           neglogP    = -log10(Adjusted.P.value),
           gene_count = str_count(Genes, ";") + 1,
           label_mode = if (!is.null(inside_terms)) {
             ifelse(Term %in% inside_terms, "inside", "outside")
           } else {
             ifelse(row_number() <= inside_n, "inside", "outside")
           })
  df$Term_short <- factor(df$Term_short, levels = rev(df$Term_short))

  ggplot(df, aes(x = gene_count, y = Term_short, fill = neglogP)) +
    geom_col(width = 0.75) +
    geom_text(data = filter(df, label_mode == "inside"),
              aes(label = Term_short), hjust = 1.05, size = 6,
              fontface = "bold", colour = "white") +
    geom_text(data = filter(df, label_mode == "outside"),
              aes(label = Term_short), hjust = -0.05, size = 6,
              fontface = "bold", colour = "black") +
    scale_fill_gradientn(colours = enrichment_palettes[[palette]]) +
    guides(fill = guide_colorbar(title = expression(-Log[10]~adj.~italic(P)))) +
    labs(title = title, x = "Gene count", y = NULL) +
    theme_classic(base_size = 16) +
    theme(plot.title = element_text(size = 15),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),
          legend.position = "right", legend.justification = "top",
          plot.margin = margin(t = 10, r = 40, b = 10, l = 10))
}

#' Strip the "(GO:0000000)" suffix from Enrichr GO term names
strip_go_id <- function(df) mutate(df, Term = str_replace(Term, "\\s*\\(GO:\\d+\\)$", ""))

# -----------------------------------------------------------------------------
# Volcano plots
# -----------------------------------------------------------------------------

#' Volcano plot with custom colours for up / down / not significant
#'
#' @param res   data frame with log fold-change and P value columns
#' @param x,y   column names for effect size and (adjusted) P value
#' @param lab   vector of labels (gene/protein names)
#' @param fc_cut absolute effect-size cut-off used for colouring (0 for slopes)
plot_volcano <- function(res, x, y, lab, select_lab = NULL, p_cut = 0.05, fc_cut = 0.25,
                         up_col = "#8152A2", down_col = "#1b9e77", ns_col = "grey80",
                         xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
                         label_size = 7, point_size = 2) {
  up   <- res[[x]] >  fc_cut & res[[y]] < p_cut
  down <- res[[x]] < -fc_cut & res[[y]] < p_cut
  keyvals <- ifelse(up, up_col, ifelse(down, down_col, ns_col))
  keyvals[is.na(keyvals)] <- ns_col
  names(keyvals) <- ifelse(keyvals == up_col, "Up",
                           ifelse(keyvals == down_col, "Down", "NS"))

  args <- list(toptable = res, lab = lab, x = x, y = y, selectLab = select_lab,
               pCutoff = p_cut, FCcutoff = fc_cut, colCustom = keyvals,
               title = "", subtitle = "", caption = "",
               pointSize = point_size, labSize = label_size,
               boxedLabels = TRUE, drawConnectors = TRUE, widthConnectors = 0.5,
               colConnectors = "black", max.overlaps = Inf,
               legendPosition = "none", gridlines.major = FALSE, gridlines.minor = FALSE)
  if (!is.null(xlim)) args$xlim <- xlim
  if (!is.null(ylim)) args$ylim <- ylim
  if (!is.null(xlab)) args$xlab <- xlab
  if (!is.null(ylab)) args$ylab <- ylab
  do.call(EnhancedVolcano::EnhancedVolcano, args)
}

#' Top n up- and downregulated genes (by adjusted P) of a Seurat FindMarkers table
top_genes_sn <- function(de, n = 3, p_cut = 0.05) {
  sig <- de[de$p_val_adj < p_cut, ]
  up   <- rownames(sig[sig$avg_log2FC > 0, ])[order(sig$p_val_adj[sig$avg_log2FC > 0])]
  down <- rownames(sig[sig$avg_log2FC < 0, ])[order(sig$p_val_adj[sig$avg_log2FC < 0])]
  c(head(up, n), head(down, n))
}

# -----------------------------------------------------------------------------
# Gene-set overlap
# -----------------------------------------------------------------------------

#' One-sided Fisher's exact test for enrichment of overlap between two gene sets
overlap_fisher <- function(set_a, set_b, universe) {
  set_a <- unique(set_a); set_b <- unique(set_b)
  n_univ <- length(unique(c(set_a, set_b, universe)))
  k <- length(intersect(set_a, set_b))
  tab <- matrix(c(k, length(set_a) - k,
                  length(set_b) - k, n_univ - (length(set_a) + length(set_b) - k)),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("In_A", "Not_in_A"), c("In_B", "Not_in_B")))
  list(table = tab, overlap = intersect(set_a, set_b),
       test = fisher.test(tab, alternative = "greater"))
}
