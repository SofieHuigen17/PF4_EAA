# PF4 and epigenetic age acceleration in human carotid atherosclerotic plaques

This repository contains the R Markdown notebooks used to generate the results,
figures and supplementary tables of the manuscript. **Patient-level data are not
included** 

---

## Overview

Plaque epigenetic age acceleration (EAA; Horvath clock on HM450k DNA methylation)
was linked to plaque cell states, proteome and endothelial PF4 signalling using
four complementary datasets:

| # | Notebook | Data | Main outputs |
|---|---|---|---|
| 01 | [`01_snRNAseq_EAA.Rmd`](01_snRNAseq_EAA.Rmd) | snRNA-seq of 8 carotid plaques (Athero-Express) | Fig. 1, Fig. 2, Fig. S1-S9, Tables S1-S8 |
| 02 | [`02_proteomics_EAA.Rmd`](02_proteomics_EAA.Rmd) | ECM-enriched plaque proteomics (n = 108) and plaque bulk RNA-seq | Fig. 3, Tables S9-S10 |
| 03 | [`03_EC_PF4_TGFb_stimulation_and_overlap.Rmd`](03_EC_PF4_TGFb_stimulation_and_overlap.Rmd) | Bulk RNA-seq of primary plaque ECs stimulated with PF4 / TGF-β1 (n = 5 donors); integration with 01 | Fig. 4, Fig. 5, Tables S11-S13 |
| 04 | [`04_spatial_Xenium.Rmd`](04_spatial_Xenium.Rmd) | Public Xenium dataset (Pauli *et al.*) | Fig. 6, Fig. S10-S16 |

Notebooks 01 and 02 are independent. Notebook 03 uses outputs of 01, and
notebook 04 uses the responsive gene sets from 03:

```
01_snRNAseq_EAA ──────────────┐
                              ├──► 03_EC_PF4_TGFb_stimulation_and_overlap ──► 04_spatial_Xenium
02_proteomics_EAA (standalone)│
```

### Repository structure

```
.
├── README.md
├── 01_snRNAseq_EAA.Rmd
├── 02_proteomics_EAA.Rmd
├── 03_EC_PF4_TGFb_stimulation_and_overlap.Rmd
├── 04_spatial_Xenium.Rmd
└── R/
    └── helper_functions.R         # shared plotting / enrichment / overlap functions
```

Each notebook writes its figures and tables to `output/<notebook>/figures` and
`output/<notebook>/tables`. Input locations are set in the `params:` block at the
top of every notebook.

---

## System requirements

- **R** 4.4.x (snRNA-seq QC, integration and clustering) and R 4.5.1 / Bioconductor 3.21 (all other analyses); run on Windows 11 and macOS
- RAM: ≥ 16 GB recommended for notebooks 01 and 04
- No non-standard hardware is required

### Main R packages

| Package | Version | Used in |
|---|---|---|
| Seurat / SeuratObject | 5.1–5.2 / 5.0.2 (QC, integration, clustering); 5.3.0 / 5.2.0 (downstream) | 01, 03 |
| DoubletFinder | 2.0.6 | 01 |
| harmony | 1.2.3 | 01 |
| clustree | 0.5.1 | 01 |
| SingleR / celldex | [to confirm] / 1.16.0 (annotation); SingleR 2.10.0 (downstream) | 01 |
| scuttle | 1.18.0 | 01 |
| MAST | 1.33.0 | 01 |
| DESeq2 | 1.48.1 | 03 |
| SummarizedExperiment / limma / org.Hs.eg.db | 1.38.1 / 3.64.3 / 3.21.0 | 03 |
| enrichR | 3.4 | 01-03 |
| EnhancedVolcano | 1.26.0 | 01-03 |
| zellkonverter / SingleCellExperiment / RANN | 1.17.4 / 1.30.1 / 2.6.2 | 04 |
| ggplot2 / ggpubr / gghalves / ggalluvial / VennDiagram / patchwork / cowplot | 4.0.1 / 1.0.0 / 0.1.4 / 0.12.6 / 1.8.2 / 1.3.2 / 1.2.0 | all |
| dplyr / tidyr / openxlsx / tableone / rmarkdown | 1.2.1 / 1.3.2 / 4.2.8 / 0.13.2 / 2.32 | all |

The exact versions of all packages are reported in the `sessionInfo()` output at
the end of each rendered notebook.

Upstream processing (not part of this repository): CeleScope v2.9.0 (snRNA-seq
alignment, GRCh38/Ensembl 99), BWA v0.7.13 (CEL-Seq2 bulk RNA-seq, Ensembl v84),
Proteome Discoverer v2.5 (proteomics), `cgageR` v0.1.0 (Horvath epigenetic age).

## Installation

```r
install.packages(c("tidyverse", "ggpubr", "gghalves", "ggalluvial", "enrichR", "openxlsx", "tableone",
                   "VennDiagram", "RANN", "clustree", "harmony", "Seurat", "BiocManager", "remotes"))
BiocManager::install(c("DESeq2", "limma", "org.Hs.eg.db", "SingleR", "celldex", "scuttle",
                       "MAST", "EnhancedVolcano", "zellkonverter", "SingleCellExperiment"))
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
```

> **Reproducing the clustering.** QC, Harmony integration and clustering of the
> snRNA-seq data were run in February 2025 with R 4.4, Seurat 5.1–5.2
> (SeuratObject 5.0.2) and harmony 1.2.3. Later versions (in particular
> harmony ≥ 2.0) may give slightly different clusters, which changes the
> cluster-to-cell-type mapping in notebook 01 (section 7d). To reproduce the
> published clusters, install the original harmony version:
>
> ```r
> remotes::install_version("harmony", version = "1.2.3")
> ```
>
> Downstream analyses (differential expression, enrichment, figures) were run
> with R 4.5.1 / Bioconductor 3.21.

Typical installation time: ~30-60 min on a standard desktop computer.

## Instructions for use

1. Clone this repository and open it as the working directory in R / RStudio.
2. Place the input data in `data/` (or adapt the paths in the `params:` block of each notebook).
3. Render the notebooks in order:

```r
rmarkdown::render("01_snRNAseq_EAA.Rmd")
rmarkdown::render("02_proteomics_EAA.Rmd")
rmarkdown::render("03_EC_PF4_TGFb_stimulation_and_overlap.Rmd")
rmarkdown::render("04_spatial_Xenium.Rmd")
```

Approximate run times on a standard desktop: 01 ~2-4 h (JackStraw, DoubletFinder
and MAST are the slow steps; set `run_slow_steps: false` to skip the diagnostic
steps), 02 < 10 min, 03 ~15 min, 04 ~20 min. Enrichr analyses require an internet
connection; results may differ slightly if the Enrichr libraries are updated.

### Expected input files

| Notebook | File | Content |
|---|---|---|
| 01 | `data/snRNAseq/<sample>/` | CeleScope filtered matrices (`matrix.mtx(.gz)`, `features/genes.tsv(.gz)`, `barcodes.tsv(.gz)`) |
| 01 | `data/snRNAseq/sample_metadata.csv` | `sample_id, label, folder, sex, batch, age, age_horvath, EAA, EAA_dich, bmi, diabetes` (`label` = anonymous P1-P8) |
| 01 | `data/reference/Pan_et_al_human.rds` | Pan *et al.* human plaque scRNA-seq reference with `Author_Provided` labels |
| 02 | `data/proteomics/plaque_proteomics_EAA.rds` | Patients × (`STUDY_NUMBER, Age, Sex, EAA, EAA.dich`, proteins) |
| 02 | `data/proteomics/plaque_bulkRNA_expression.rds` | Patients × genes, log2 normalised plaque bulk RNA-seq |
| 03 | `data/EC_stimulation/EC_PF4_TGF_raw_counts.txt` | Ensembl gene × sample UMI counts (`<Condition>_<Timepoint>_<Donor>.sam.counts`) |
| 03 | `data/EC_stimulation/donor_metadata.csv` | `Donor, Sex` |
| 04 | `data/spatial/Panel1_cells_scale_10_asg_conf_0.75_final_data.h5ad` | Xenium panel 1 from Pauli *et al.* |

## Key analysis parameters

| Analysis | Parameters |
|---|---|
| snRNA-seq QC | 200 < nFeature < 3,000; mitochondrial reads < 0.2%; *MALAT1* and MT genes removed |
| Doublets | DoubletFinder, pN = 0.25, pK = 0.02, expected rate 4%, PCs 1-10 |
| Clustering | LogNormalize, 2,000 VST HVGs, Harmony (batch = sample), 20 PCs, resolution 0.9 |
| snRNA-seq DE | MAST, accelerated vs decelerated, latent variables sex + sequencing batch; adj. P < 0.01, \|log2FC\| ≥ 0.25 (figures, tables and enrichment) |
| snRNA-seq enrichment | Enrichr, genes with adj. P < 0.01; background = genes tested per cell group |
| Proteomics | `lm(EAA ~ protein + age + sex)` per protein; BH correction; nominal P < 0.05 for enrichment |
| EC stimulation DE | DESeq2, `~ Donor + Condition`, genes with total count > 30; FDR < 0.05, \|log2FC\| ≥ 0.25 |
| Spatial scores | Per-patient z-score of area-normalised expression, mean per cell; 14 PF4 / 26 TGF-β panel genes |
| Spatial statistics | Wilcoxon rank-sum (cell level, per sample); paired Wilcoxon signed-rank (sample-level medians) |

---

## Data availability

- **snRNA-seq, plaque proteomics, plaque bulk RNA-seq and EC stimulation RNA-seq**
  are derived from the Athero-Express Biobank. Because of patient privacy and the
  informed consent, individual-level data cannot be shared publicly. Data are
  available upon reasonable request to.
- **Xenium spatial transcriptomics**: publicly available from Pauli *et al.*
- **Pan *et al.* scRNA-seq reference**: publicly available 

