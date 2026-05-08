############################################################
# GSE120575 — VALIDATION PIPELINE
# Melanoma TIL single-cell RNA-seq
# State-focused validation analysis
#
# Dataset:
# GSE120575
#
# Main goals:
# - Validate recurrent cytotoxic/transitional states
# - Validate overlap of activation/exhaustion programs
# - Validate recurrent genes:
#   CCL4, CCL4L2, GZMK, GNLY, FGFBP2,
#   PRF1, NKG7, IFNG
#
# Outputs:
# - UMAPs
# - Heatmaps
# - Dotplots
# - Feature plots
# - Module scoring
# - Correlation analyses
# - Patient-level analyses
# - Volcano plots
# - DE tables
# - Publication-ready figures
############################################################

############################################################
# STEP 1 — Libraries
############################################################

pkgs <- c(
  "Seurat",
  "Matrix",
  "data.table",
  "dplyr",
  "ggplot2",
  "patchwork",
  "pheatmap",
  "RColorBrewer",
  "reshape2",
  "future"
)

to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

if(length(to_install) > 0){
  install.packages(to_install)
}

invisible(lapply(pkgs, library, character.only = TRUE))

options(future.globals.maxSize = 50 * 1024^3)

set.seed(1234)

############################################################
# STEP 2 — Paths
############################################################

BASE_DIR <- "D:/Papers/2026/In Vitro Paper/GSE120575_melanoma_checkpoint_TIL"

RAW_DIR <- file.path(BASE_DIR, "raw_data")

OUT_DIR <- file.path(BASE_DIR, "analysis_outputs")
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDS_DIR <- file.path(OUT_DIR, "rds")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

expr_file <- file.path(
  RAW_DIR,
  "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt"
)

############################################################
# STEP 3 — Theme
############################################################

theme_pub <- function(){

  theme_minimal(base_size = 14) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5
      ),

      axis.title = element_text(face = "bold"),

      axis.text = element_text(color = "black"),

      legend.title = element_text(face = "bold"),

      strip.text = element_text(face = "bold"),

      panel.grid.minor = element_blank()
    )
}

save_plot <- function(
    plot,
    filename,
    width = 8,
    height = 6,
    dpi = 300
){

  ggsave(
    filename = file.path(FIG_DIR, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

############################################################
# STEP 4 — Read TPM matrix
############################################################

cat("\nReading header lines...\n")

header_lines <- readLines(expr_file, n = 2)

cell_ids <- strsplit(
  header_lines[1],
  "\t",
  fixed = TRUE
)[[1]][-1]

sample_ids <- strsplit(
  header_lines[2],
  "\t",
  fixed = TRUE
)[[1]][-1]

cat("\nCells detected:", length(cell_ids), "\n")
cat("Unique samples:", length(unique(sample_ids)), "\n")

############################################################
# STEP 5 — Read expression matrix
############################################################

cat("\nReading expression matrix...\n")

expr_dt <- fread(
  expr_file,
  sep = "\t",
  skip = 2,
  header = FALSE,
  data.table = TRUE,
  showProgress = TRUE
)

gene_names <- expr_dt[[1]]

expr_dt[[1]] <- NULL

############################################################
# IMPORTANT FIX
############################################################

expr_dt <- expr_dt[, 1:length(cell_ids), with = FALSE]

expr_mat <- as.matrix(expr_dt)

storage.mode(expr_mat) <- "numeric"

rownames(expr_mat) <- make.unique(gene_names)
colnames(expr_mat) <- cell_ids

rm(expr_dt)

gc()

############################################################
# STEP 6 — Metadata
############################################################

meta <- data.frame(
  cell_barcode = cell_ids,
  sample_label = sample_ids,
  stringsAsFactors = FALSE
)

meta$treatment <- ifelse(
  grepl("^Pre_", meta$sample_label),
  "Pre",
  "Post"
)

meta$patient_id <- sub(
  "^(Pre|Post)_",
  "",
  meta$sample_label
)

rownames(meta) <- meta$cell_barcode

############################################################
# STEP 7 — Create Seurat object
############################################################

cat("\nCreating Seurat object...\n")

seu <- CreateSeuratObject(
  counts = expr_mat,
  meta.data = meta,
  project = "GSE120575"
)

rm(expr_mat)

gc()

############################################################
# STEP 8 — Basic preprocessing
############################################################

DefaultAssay(seu) <- "RNA"

seu <- NormalizeData(seu)

seu <- FindVariableFeatures(
  seu,
  selection.method = "vst",
  nfeatures = 3000
)

seu <- ScaleData(seu)

############################################################
# STEP 9 — PCA / UMAP / Clustering
############################################################

seu <- RunPCA(seu)

seu <- FindNeighbors(
  seu,
  dims = 1:30
)

seu <- FindClusters(
  seu,
  resolution = 0.6
)

seu <- RunUMAP(
  seu,
  dims = 1:30
)

############################################################
# STEP 10 — Initial UMAPs
############################################################

p1 <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) +
  ggtitle("GSE120575 Clusters") +
  theme_pub()

save_plot(
  p1,
  "Figure_3A_UMAP_clusters.png",
  8,
  6
)

p2 <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "treatment"
) +
  ggtitle("Pre vs Post Treatment") +
  theme_pub()

save_plot(
  p2,
  "Figure_3B_UMAP_treatment.png",
  8,
  6
)

############################################################
# STEP 11 — Canonical markers
############################################################

marker_genes <- c(
  "CD3D",
  "CD3E",
  "IL7R",
  "CCR7",
  "LTB",
  "GZMK",
  "NKG7",
  "GNLY",
  "PRF1",
  "GZMB",
  "FGFBP2",
  "CCL4",
  "CCL4L2",
  "PDCD1",
  "TIGIT",
  "LAG3",
  "TOX",
  "CXCL13",
  "MKI67",
  "TOP2A"
)

marker_genes <- marker_genes[
  marker_genes %in% rownames(seu)
]

############################################################
# STEP 12 — Dotplot
############################################################

p_dot <- DotPlot(
  seu,
  features = marker_genes
) +
  RotatedAxis() +
  ggtitle("Canonical T-cell markers") +
  theme_pub()

save_plot(
  p_dot,
  "Figure_3C_Dotplot_markers.png",
  14,
  7
)

############################################################
# STEP 13 — Heatmap
############################################################

top_markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

fwrite(
  top_markers,
  file.path(
    TAB_DIR,
    "Table_top_markers.csv"
  )
)

top10 <- top_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  )

heatmap_genes <- unique(top10$gene)

heatmap_genes <- heatmap_genes[
  heatmap_genes %in% rownames(seu)
]

p_heat <- DoHeatmap(
  seu,
  features = heatmap_genes,
  size = 3
) +
  ggtitle("Top cluster markers")

save_plot(
  p_heat,
  "Figure_3D_Heatmap_top_markers.png",
  12,
  10
)

############################################################
# STEP 14 — Functional modules
############################################################

cytotoxicity_genes <- c(
  "NKG7",
  "GNLY",
  "PRF1",
  "GZMB",
  "GZMK",
  "IFNG",
  "FGFBP2",
  "CCL4",
  "CCL4L2"
)

exhaustion_genes <- c(
  "PDCD1",
  "LAG3",
  "TIGIT",
  "HAVCR2",
  "TOX",
  "CXCL13"
)

prolif_genes <- c(
  "MKI67",
  "TOP2A",
  "TYMS"
)

activation_genes <- c(
  "IFNG",
  "TNF",
  "NFKBIA",
  "JUN",
  "FOS"
)

interferon_genes <- c(
  "IFITM1",
  "IFITM2",
  "IFITM3",
  "STAT1",
  "ISG15"
)

############################################################
# STEP 15 — Module scoring
############################################################

score_sets <- list(
  cytotoxicity_genes,
  exhaustion_genes,
  prolif_genes,
  activation_genes,
  interferon_genes
)

score_names <- c(
  "Cytotoxicity",
  "Exhaustion",
  "Proliferation",
  "Activation",
  "Interferon"
)

for(i in 1:length(score_sets)){

  genes <- score_sets[[i]]

  genes <- genes[
    genes %in% rownames(seu)
  ]

  seu <- AddModuleScore(
    seu,
    features = list(genes),
    name = paste0(score_names[i], "_Score")
  )
}

############################################################
# STEP 16 — FeaturePlots
############################################################

score_features <- c(
  "Cytotoxicity_Score1",
  "Exhaustion_Score1",
  "Proliferation_Score1",
  "Activation_Score1",
  "Interferon_Score1"
)

p_scores <- FeaturePlot(
  seu,
  features = score_features,
  ncol = 3
) &
  theme_pub()

save_plot(
  p_scores,
  "Figure_4A_Module_featureplots.png",
  14,
  8
)

############################################################
# STEP 17 — Violin plots
############################################################

p_vln <- VlnPlot(
  seu,
  features = score_features,
  group.by = "seurat_clusters",
  pt.size = 0
) &
  theme_pub()

save_plot(
  p_vln,
  "Figure_4B_Module_violinplots.png",
  16,
  8
)

############################################################
# STEP 18 — Correlation matrix
############################################################

score_df <- FetchData(
  seu,
  vars = score_features
)

cor_mat <- cor(score_df)

write.csv(
  cor_mat,
  file.path(
    TAB_DIR,
    "Module_correlations.csv"
  )
)

png(
  file.path(
    FIG_DIR,
    "Figure_4C_Module_correlations.png"
  ),
  width = 1800,
  height = 1600,
  res = 300
)

pheatmap(
  cor_mat,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  display_numbers = TRUE
)

dev.off()

############################################################
# STEP 19 — Scatterplots
############################################################

p_scatter1 <- ggplot(
  score_df,
  aes(
    Cytotoxicity_Score1,
    Exhaustion_Score1
  )
) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm") +
  theme_pub() +
  ggtitle(
    "Cytotoxicity vs Exhaustion"
  )

save_plot(
  p_scatter1,
  "Figure_4D_Cyto_vs_Exhaustion.png",
  7,
  6
)

p_scatter2 <- ggplot(
  score_df,
  aes(
    Cytotoxicity_Score1,
    Activation_Score1
  )
) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm") +
  theme_pub() +
  ggtitle(
    "Cytotoxicity vs Activation"
  )

save_plot(
  p_scatter2,
  "Figure_4E_Cyto_vs_Activation.png",
  7,
  6
)

############################################################
# STEP 20 — Patient-level analyses
############################################################

patient_summary <- seu@meta.data %>%
  group_by(
    patient_id,
    treatment,
    seurat_clusters
  ) %>%
  summarise(
    cells = n()
  ) %>%
  ungroup()

write.csv(
  patient_summary,
  file.path(
    TAB_DIR,
    "Patient_cluster_abundance.csv"
  ),
  row.names = FALSE
)

p_patient <- ggplot(
  patient_summary,
  aes(
    patient_id,
    cells,
    fill = seurat_clusters
  )
) +
  geom_bar(
    stat = "identity"
  ) +
  facet_wrap(~treatment) +
  theme_pub() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  ) +
  ggtitle(
    "Patient-level cluster abundance"
  )

save_plot(
  p_patient,
  "Figure_5A_Patient_cluster_abundance.png",
  14,
  8
)

############################################################
# STEP 21 — Differential expression
############################################################

Idents(seu) <- "treatment"

de_treatment <- FindMarkers(
  seu,
  ident.1 = "Post",
  ident.2 = "Pre",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

de_treatment$gene <- rownames(de_treatment)

write.csv(
  de_treatment,
  file.path(
    TAB_DIR,
    "DE_Post_vs_Pre.csv"
  ),
  row.names = FALSE
)

############################################################
# STEP 22 — Volcano plot
############################################################

de_treatment$significant <- ifelse(
  de_treatment$p_val_adj < 0.05 &
    abs(de_treatment$avg_log2FC) > 0.5,
  "yes",
  "no"
)

p_volcano <- ggplot(
  de_treatment,
  aes(
    avg_log2FC,
    -log10(p_val_adj),
    color = significant
  )
) +
  geom_point(alpha = 0.7) +
  theme_pub() +
  ggtitle(
    "Post vs Pre differential expression"
  )

save_plot(
  p_volcano,
  "Figure_6A_Volcano_Post_vs_Pre.png",
  8,
  6
)

############################################################
# STEP 23 — Recurrent genes
############################################################

recurrent_genes <- c(
  "CCL4",
  "CCL4L2",
  "GZMK",
  "GNLY",
  "FGFBP2",
  "PRF1",
  "NKG7"
)

recurrent_genes <- recurrent_genes[
  recurrent_genes %in% rownames(seu)
]

############################################################
# STEP 24 — FeaturePlots recurrent genes
############################################################

p_recurrent <- FeaturePlot(
  seu,
  features = recurrent_genes,
  ncol = 3
) &
  theme_pub()

save_plot(
  p_recurrent,
  "Figure_6B_Recurrent_gene_featureplots.png",
  14,
  10
)

############################################################
# STEP 25 — Recurrent heatmap
############################################################

avg_exp <- AverageExpression(
  seu,
  features = recurrent_genes
)

avg_exp <- avg_exp$RNA

png(
  file.path(
    FIG_DIR,
    "Figure_6C_Recurrent_gene_heatmap.png"
  ),
  width = 1800,
  height = 1200,
  res = 300
)

pheatmap(
  avg_exp,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)

dev.off()

############################################################
# STEP 26 — Save object
############################################################

saveRDS(
  seu,
  file.path(
    RDS_DIR,
    "GSE120575_validation_seurat.rds"
  )
)

############################################################
# STEP 27 — Session info
############################################################

writeLines(
  capture.output(sessionInfo()),
  file.path(
    OUT_DIR,
    "sessionInfo.txt"
  )
)

############################################################
# FINISHED
############################################################

cat("\nVALIDATION PIPELINE FINISHED.\n")
cat("\nOutputs saved in:\n")
cat(OUT_DIR, "\n")
