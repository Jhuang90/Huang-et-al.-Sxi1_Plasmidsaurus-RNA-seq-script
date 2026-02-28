## =========================
## TMM logCPM + grouped bar plot - 3 SEPARATE PANELS
## Updated: Custom colors, clearer axis labels, split into 3 plots
## =========================

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
})

## Input
infile <- "data/featureCounts_gene_level.txt"   # <-- set your input file path here

## Output
outdir <- "results/IndividualPLOT"              # <-- set your output folder path here
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(infile))

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

## Optional: drop patterns
## Add Geneid substrings you want to exclude (e.g. known artefacts or
## annotation noise specific to your reference genome).
drop_patterns <- c(
  ## "GENE001",   # example: a known artefact gene
  ## "GENE002"    # example: another gene to exclude
)

## Genes you want to plot (edit as needed).
## Use the exact Geneid strings as they appear in your featureCounts table.
## Genes are split across 3 panels: panel1 = genes 1-5, panel2 = 6-10, panel3 = 11+
plot_genes <- c(
  ## "gene-GENE001",   # example gene of interest
  ## "gene-GENE002"    # example gene of interest
)

## 1 Read featureCounts table
fc <- read.delim(infile, comment.char = "#", check.names = FALSE)
meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
stopifnot(all(meta_cols %in% colnames(fc)))

counts <- fc[, setdiff(colnames(fc), meta_cols), drop = FALSE]
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

## 4 edgeR standard practice + TMM
y <- DGEList(counts = counts, group = group)
keep <- filterByExpr(y)
message("Rows before filterByExpr: ", nrow(y))
message("Rows kept after filterByExpr: ", sum(keep))

y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

## Show TMM normalization factors
message("\n=== TMM Normalization Factors ===")
print(y$samples)
message("Effective library sizes = lib.size × norm.factors")
message("These factors are used by cpm() to calculate normalized CPM values\n")

## 5 logCPM matrix (TMM-adjusted)
logCPM <- cpm(y, log = TRUE, prior.count = 1)

write.table(
  logCPM,
  file.path(outdir, "TMM_logCPM_filtered.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

## 6 Subset genes and build long table for plotting
present <- intersect(plot_genes, rownames(logCPM))
missing <- setdiff(plot_genes, rownames(logCPM))

message("Plot genes requested: ", length(plot_genes))
message("Plot genes present: ", length(present))

if (length(missing) > 0) {
  message("Missing genes (not found in rownames(logCPM)):")
  message(paste(missing, collapse = "\n"))
}

sub <- logCPM[present, , drop = FALSE]

df <- data.frame(
  gene = rep(rownames(sub), times = ncol(sub)),
  sample = rep(colnames(sub), each = nrow(sub)),
  logCPM = as.vector(sub),
  stringsAsFactors = FALSE
)

df$group <- ifelse(grepl("_WT-", df$sample), "WT",
                   ifelse(grepl("_S-", df$sample), "S", NA))
df$group <- factor(df$group, levels = c("WT", "S"))
df$gene <- factor(df$gene, levels = present)

## 7 Perform statistical tests (t-test for each gene)
stat_results <- data.frame()

for (g in present) {
  wt_vals <- df$logCPM[df$gene == g & df$group == "WT"]
  s_vals <- df$logCPM[df$gene == g & df$group == "S"]
  
  if (length(wt_vals) >= 2 && length(s_vals) >= 2) {
    t_result <- t.test(wt_vals, s_vals)
    stat_results <- rbind(stat_results, data.frame(
      gene = g,
      WT_mean = mean(wt_vals),
      S_mean = mean(s_vals),
      log2FC = mean(s_vals) - mean(wt_vals),
      p_value = t_result$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

# Adjust p-values for multiple testing (FDR)
if (nrow(stat_results) > 0) {
  stat_results$p_adj <- p.adjust(stat_results$p_value, method = "BH")
  stat_results$significance <- ifelse(stat_results$p_adj < 0.001, "***",
                                      ifelse(stat_results$p_adj < 0.01, "**",
                                             ifelse(stat_results$p_adj < 0.05, "*", "ns")))
  
  # Save statistics
  write.table(
    stat_results,
    file.path(outdir, "gene_statistics.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  message("\n=== Statistical Testing Results ===")
  print(stat_results)
}

## 8 Create THREE SEPARATE PLOTS
## Split genes into 3 groups:
## Panel 1: genes 1-5
## Panel 2: genes 6-10
## Panel 3: genes 11-end

custom_colors <- c("WT" = "#f1e2cc", "S" = "#cbd5e8")
dot_colors <- c("WT" = "#8B7355", "S" = "#4682B4")

# Define gene groups
n_present <- length(present)
gene_groups <- list(
  panel1 = present[1:min(5, n_present)],
  panel2 = present[6:min(10, n_present)],
  panel3 = present[11:n_present]
)

# Remove empty groups
gene_groups <- gene_groups[sapply(gene_groups, length) > 0]

# Create plot for each group
for (panel_name in names(gene_groups)) {
  panel_genes <- gene_groups[[panel_name]]
  
  # Subset data for this panel
  df_panel <- df[df$gene %in% panel_genes, ]
  df_panel$gene <- factor(df_panel$gene, levels = panel_genes)
  
  # Create plot
  p <- ggplot(df_panel, aes(x = gene, y = logCPM, fill = group)) +
    stat_summary(
      fun = mean,
      geom = "col",
      position = position_dodge(width = 0.82),
      width = 0.72,
      alpha = 0.8
    ) +
    geom_point(
      aes(color = group),
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.82),
      size = 2.5,
      alpha = 1,
      shape = 21,
      fill = NA,
      stroke = 0.5
    ) +
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = dot_colors) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
      panel.grid.major.x = element_blank(),
      legend.position = "top"
    ) +
    labs(
      x = "Gene",
      y = "log₂ CPM (TMM normalized)"
    )
  
  # Determine appropriate width based on number of genes
  plot_width <- 2 + (length(panel_genes) * 1.5)
  
  # Save plot
  ggsave(
    filename = file.path(outdir, paste0("BarPlot_TMM_logCPM_", panel_name, ".png")),
    plot = p,
    width = plot_width,
    height = 4.5,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(outdir, paste0("BarPlot_TMM_logCPM_", panel_name, ".pdf")),
    plot = p,
    width = plot_width,
    height = 4.5
  )
  
  message("Created plot for ", panel_name, " with ", length(panel_genes), " genes")
}

message("\nDone. Outputs in: ", outdir)
message("Removed rows matching drop_patterns count: ", length(present_drop))
message("Created ", length(gene_groups), " separate plots:")
for (panel_name in names(gene_groups)) {
  message("  - ", panel_name, ": ", length(gene_groups[[panel_name]]), " genes")
}
