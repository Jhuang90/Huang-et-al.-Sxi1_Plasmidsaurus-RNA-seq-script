## edgeR differential expression + EnhancedVolcano
## only label genes that are BOTH selected by label_patterns AND are DEG
## shorten y axis by limiting to the data range with a small headroom
## remove outlier outline layer
## keep your previous color code
## keep drop and label partial string match
## keep debug prints
## keep DEG tables, summary, gene lists

## Input
infile <- "data/featureCounts_gene_level.txt"   # <-- set your input file path here

## Output
outdir <- "results/DE_edgeR"                    # <-- set your output folder path here
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## Install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_bioc <- c("edgeR", "EnhancedVolcano")
for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}

library(edgeR)
library(EnhancedVolcano)

## Parameters
lfc_cut <- 1
fdr_cut <- 0.05

## Helper for partial matching, literal substring match
match_any_fixed <- function(x, patterns) {
  if (length(patterns) == 0) return(rep(FALSE, length(x)))
  Reduce(`|`, lapply(patterns, function(p) grepl(p, x, fixed = TRUE)))
}

## Sample naming convention:
## Column names in the featureCounts table are expected to contain either
## "_WT-" (for wild-type samples) or "_S-" (for mutant/treatment samples).
## Example: "sample1_WT-rep1", "sample2_S-rep2"
## Adjust the grepl() patterns in section 3 if your naming differs.

## 1 Patterns to drop if Geneid contains any of these substrings.
## Add gene IDs or partial strings you want to exclude (e.g. known artefacts,
## mitochondrial genes, or annotation noise specific to your reference genome).
drop_patterns <- c(
  ## "GENE001",   # example: a known artefact gene
  ## "GENE002"    # example: another gene to exclude
)

## 2 Patterns to label on the volcano plot if Geneid contains any of these
## substrings AND the gene passes DEG cutoffs.
## Add gene IDs or partial strings for genes of biological interest.
label_patterns <- c(
  ## "GENE003",   # example: gene of interest
  ## "GENE004"    # example: another gene of interest
)



## 1 Read featureCounts table
fc <- read.delim(infile, comment.char = "#", check.names = FALSE)
meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
stopifnot(all(meta_cols %in% colnames(fc)))

counts <- fc[, setdiff(colnames(fc), meta_cols)]
rownames(counts) <- fc$Geneid
counts <- as.matrix(counts)
mode(counts) <- "numeric"

## 2 Drop genes by partial match
rn <- rownames(counts)
drop_flag <- match_any_fixed(rn, drop_patterns)
present_drop <- rn[drop_flag]

message("Drop patterns: ", paste(drop_patterns, collapse = ", "))
message("Rows matched for dropping: ", length(present_drop))
if (length(present_drop) > 0) {
  message(paste(head(present_drop, 50), collapse = "\n"))
  if (length(present_drop) > 50) message("More omitted")
  counts <- counts[!drop_flag, , drop = FALSE]
} else {
  message("NOTE no rows matched drop_patterns")
}

samples <- colnames(counts)

## 3 Define groups from sample names
group <- ifelse(grepl("_WT-", samples), "WT",
                ifelse(grepl("_S-", samples), "S", NA))

if (any(is.na(group))) {
  stop(
    "Some samples could not be assigned to WT or S. Edit the group definition. Samples:\n",
    paste(samples, collapse = "\n")
  )
}

group <- factor(group, levels = c("WT", "S"))

write.table(
  data.frame(sample = samples, group = group),
  file.path(outdir, "sample_info.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## 4 edgeR standard practice with filterByExpr default
y <- DGEList(counts = counts, group = group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ group)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "groupS")

## 5 Results table
res <- topTags(qlf, n = Inf)$table
res$FDR <- p.adjust(res$PValue, method = "BH")

write.table(
  res,
  file.path(outdir, "edgeR_DE_S_vs_WT_all_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = TRUE
)

## 6 DEG classification using your cutoffs
deg_flag <- (res$FDR <= fdr_cut) & (abs(res$logFC) >= lfc_cut)

deg_tbl <- data.frame(
  gene = rownames(res),
  logFC = res$logFC,
  logCPM = res$logCPM,
  F = res$F,
  PValue = res$PValue,
  FDR = res$FDR,
  direction = ifelse(res$logFC >= lfc_cut, "Up_in_S",
                     ifelse(res$logFC <= -lfc_cut, "Down_in_S", "Not_DEG")),
  pass_cutoff = deg_flag,
  stringsAsFactors = FALSE
)

deg_pass <- deg_tbl[deg_tbl$pass_cutoff, ]
deg_pass <- deg_pass[order(deg_pass$FDR, -abs(deg_pass$logFC)), ]

write.table(
  deg_pass,
  file.path(outdir, paste0("DEG_pass_logFC", lfc_cut, "_FDR", fdr_cut, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

n_up <- sum(deg_pass$direction == "Up_in_S")
n_down <- sum(deg_pass$direction == "Down_in_S")

deg_summary <- data.frame(
  cutoff_logFC = lfc_cut,
  cutoff_FDR = fdr_cut,
  Up_in_S = n_up,
  Down_in_S = n_down,
  Total_DEG = n_up + n_down,
  stringsAsFactors = FALSE
)

write.table(
  deg_summary,
  file.path(outdir, "DEG_summary_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

up_genes <- deg_pass$gene[deg_pass$direction == "Up_in_S"]
down_genes <- deg_pass$gene[deg_pass$direction == "Down_in_S"]

write.table(
  data.frame(gene = up_genes),
  file.path(outdir, "DEG_genes_Up_in_S.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

write.table(
  data.frame(gene = down_genes),
  file.path(outdir, "DEG_genes_Down_in_S.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

write.table(
  data.frame(
    category = c(rep("Up_in_S", length(up_genes)), rep("Down_in_S", length(down_genes))),
    gene = c(up_genes, down_genes),
    stringsAsFactors = FALSE
  ),
  file.path(outdir, "DEG_gene_lists_by_category.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## 7 EnhancedVolcano
## label only genes that are DEG and match label_patterns
lab_match <- match_any_fixed(rownames(res), label_patterns)
sel_deg <- rownames(res)[lab_match & deg_flag]

message("Label patterns: ", paste(label_patterns, collapse = ", "))
message("Rows matched label_patterns: ", sum(lab_match))
message("Rows that are DEG: ", sum(deg_flag))
message("Rows matched for labeling (label_patterns AND DEG): ", length(sel_deg))
if (length(sel_deg) > 0) {
  message("Matched rownames for labeling:")
  message(paste(sel_deg, collapse = "\n"))
} else {
  message("No selected DEG labels. Showing first 30 rownames(res):")
  message(paste(head(rownames(res), 30), collapse = "\n"))
}

## tighter y axis limit using plotted y = -log10(FDR)
yvals <- -log10(res$FDR)
yvals <- yvals[is.finite(yvals)]
ymax <- max(yvals, na.rm = TRUE)
ymax <- ymax * 1.05

pdf(file.path(outdir, "EnhancedVolcano_S_vs_WT_selected_DEG_labels_tightY.pdf"), width = 7.5, height = 6.5)

p <- EnhancedVolcano(
  res,
  lab = rownames(res),
  selectLab = sel_deg,
  x = "logFC",
  y = "FDR",
  xlab = "log2 Fold Change (S vs WT)",
  ylab = "FDR",
  pCutoff = fdr_cut,
  FCcutoff = lfc_cut,
  pointSize = 1.6,
  labSize = 4,
  drawConnectors = TRUE,
  widthConnectors = 0.4,
  max.overlaps = Inf,
  col = c("grey70", "grey70", "grey70", "#E0A8A9"),
  colAlpha = 0.95,
  legendPosition = "right",
  legendLabSize = 12,
  legendIconSize = 4
)

## tighten y axis to remove empty space
p <- p + ggplot2::coord_cartesian(ylim = c(0, ymax))

print(p)
dev.off()

## 8 Save normalized logCPM matrix used for QC and checks
logCPM <- cpm(y, log = TRUE, prior.count = 1)
write.table(
  logCPM,
  file.path(outdir, "TMM_logCPM_filtered.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

## 9 Messages
message("Done. Outputs in: ", outdir)
message("Removed rows matching drop_patterns count: ", length(present_drop))
message("DEG passing cutoffs. Up in S: ", n_up, " Down in S: ", n_down, " Total: ", n_up + n_down)
message("Label rows matched count: ", length(sel_deg))
message("Tight y axis ymax used: ", signif(ymax, 4))
