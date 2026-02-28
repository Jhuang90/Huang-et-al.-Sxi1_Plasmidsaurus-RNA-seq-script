## Run in RStudio
## Input: featureCounts gene-level count table
## Adjust paths below before running
infile <- "data/featureCounts_gene_level.txt"   # <-- set your input file path here

## Output folder
outdir <- "results/sample_cor_pca"              # <-- set your output folder path here
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## Install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(edgeR)
library(ggplot2)

message("Reading: ", infile)
fc <- read.delim(infile, comment.char = "#", check.names = FALSE)

meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
stopifnot(all(meta_cols %in% colnames(fc)))

counts <- fc[, setdiff(colnames(fc), meta_cols)]
rownames(counts) <- fc$Geneid
counts <- as.matrix(counts)
mode(counts) <- "numeric"

message("Counts: ", nrow(counts), " genes x ", ncol(counts), " samples")

## TMM normalization
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")

## log2 CPM on TMM normalized library sizes
logCPM <- cpm(dge, log = TRUE, prior.count = 1)

## Optional filter to reduce noise
keep <- rowSums(cpm(dge) > 1) >= 2
logCPM_f <- logCPM[keep, , drop = FALSE]
message("After filtering: ", nrow(logCPM_f), " genes kept")

## Pearson correlation among samples
cor_mat <- cor(logCPM_f, method = "pearson")
write.table(
  cor_mat,
  file = file.path(outdir, "TMM_logCPM_pearson_sample_cor.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

## Heatmap (base R)
pdf(file.path(outdir, "sample_sample_correlation_heatmap.pdf"), width = 7, height = 6)
heatmap(
  cor_mat,
  symm = TRUE,
  margins = c(10, 10),
  main = "Sample sample Pearson correlation (TMM logCPM)"
)
dev.off()

## PCA on normalized counts
pca <- prcomp(t(logCPM_f), center = TRUE, scale. = FALSE)

pve <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_lab <- paste0("PC1 (", round(100 * pve[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(100 * pve[2], 1), "%)")

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)
write.table(
  pca_df,
  file = file.path(outdir, "PCA_TMM_logCPM_scores.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, size = 3) +
  xlab(pc1_lab) + ylab(pc2_lab) +
  ggtitle("PCA on TMM normalized logCPM") +
  theme_bw()

ggsave(file.path(outdir, "PCA_TMM_logCPM_PC1_PC2.pdf"), plot = p, width = 7, height = 5)

message("Done. Outputs in: ", outdir)
