# Differential gene expression for one contrast, from Salmon quantifications.
#
# Adapted from the lab's DESeq2.R. Behaviour preserved:
#
#   * tximport with countsFromAbundance = "lengthScaledTPM", ignoreAfterBar = TRUE
#   * transcript -> gene mapping parsed out of the GENCODE-style quant.sf
#     "Name" field: field 2 of the "|"-separated header is the gene ID and
#     the second-to-last field is the gene symbol
#   * output table carries gene_id and gene_name alongside the DESeq2
#     statistics, written tab-separated without quoting or row names
#   * PCA built from vst(dds, blind = FALSE) -> plotPCA(returnData = TRUE)
#
# Deliberate differences from the original script:
#
#   * samples are taken from Snakemake (params$sample_names / sample_groups)
#     instead of list.dirs() on the working directory, and only the samples
#     belonging to this contrast are loaded
#   * the tx2gene / transcript tables are held in memory rather than written
#     to CSV and immediately read back
#   * the condition is read from the design file rather than inferred with
#     grepl() on the sample name
#   * txdbmaker / GenomicFeatures are not loaded; they were only used by the
#     commented-out GTF section
#   * the PCA is written to a PDF instead of being drawn interactively

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

# Restrict the library search path to the activated conda environment.
#
# R_LIBS_USER is commonly exported in HPC shell profiles and points at a
# personal library built against a different R version. R prepends it to
# .libPaths(), so it takes precedence over the conda environment and packages
# get loaded from the wrong place -- e.g. DESeq2 from the conda env but
# MatrixGenerics from the personal library, which fails with:
#     object 'colMeans' is not exported by 'namespace:MatrixGenerics'
# Rscript --vanilla does NOT prevent this: it ignores .Renviron files but
# still honours R_LIBS_USER inherited from the environment.
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (nzchar(conda_prefix)) {
  conda_lib <- file.path(conda_prefix, "lib", "R", "library")
  if (dir.exists(conda_lib)) {
    .libPaths(conda_lib)
  }
}
cat("Library paths in use:\n")
print(.libPaths())

suppressPackageStartupMessages({
  library("DESeq2")
  library("tximport")
  library("ggplot2")
})

cat("\nPackage versions:\n")
for (pkg in c("DESeq2", "tximport", "ggplot2", "MatrixGenerics", "BiocGenerics")) {
  cat(" ", pkg, as.character(packageVersion(pkg)), "\n")
}

files <- unlist(snakemake@input[["quant"]])
sample_names <- unlist(snakemake@params[["sample_names"]])
sample_groups <- unlist(snakemake@params[["sample_groups"]])
reference <- snakemake@params[["reference"]]
treatment <- snakemake@params[["treatment"]]
cfa <- snakemake@params[["counts_from_abundance"]]

stopifnot(length(files) == length(sample_names))
stopifnot(length(files) == length(sample_groups))
names(files) <- sample_names

cat("Contrast:", treatment, "vs", reference, "(reference =", reference, ")\n")
cat("Samples:\n")
print(data.frame(sample = sample_names, group = sample_groups))

####
## transcript -> gene mapping, parsed from the quant.sf header
####

quant <- read.table(files[1], header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)
parts <- strsplit(as.character(quant$Name), "\\|")

if (max(lengths(parts)) < 3) {
  stop(
    "The 'Name' column in ", files[1], " does not look like a GENCODE-style ",
    "header (expected fields separated by '|'). The gene ID and gene symbol ",
    "are parsed out of that header, so the Salmon index must have been built ",
    "from a FASTA with full GENCODE headers."
  )
}

tx2gene <- data.frame(
  TXNAME = vapply(parts, `[`, character(1), 1),
  GENEID = vapply(parts, `[`, character(1), 2),
  stringsAsFactors = FALSE
)

# gene_id -> gene_name lookup; the symbol is the second-to-last field
gene_map <- data.frame(
  GENEID = vapply(parts, `[`, character(1), 2),
  gene_name = vapply(parts, function(x) x[length(x) - 2], character(1)),
  stringsAsFactors = FALSE
)
gene_map <- gene_map[!duplicated(gene_map$GENEID), ]

cat("Parsed", nrow(tx2gene), "transcripts spanning",
    nrow(gene_map), "genes\n")

####
## import quantifications and run DESeq2
####

txi <- tximport(files, type = "salmon", tx2gene = tx2gene,
                countsFromAbundance = cfa, ignoreAfterBar = TRUE)

# Setting `reference` as the first factor level makes it the baseline, so
# positive log2FoldChange means higher in `treatment`.
sampleTable <- data.frame(
  sampleName = sample_names,
  group = factor(sample_groups, levels = c(reference, treatment)),
  row.names = sample_names,
  stringsAsFactors = FALSE
)

dds <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", treatment, reference))

cat("\nDESeq2 results summary:\n")
summary(res)

####
## result table: gene_id, gene_name, then the DESeq2 statistics
####

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df$gene_name <- gene_map$gene_name[match(res_df$gene_id, gene_map$GENEID)]

stat_cols <- setdiff(colnames(res_df), c("gene_id", "gene_name"))
res_df <- res_df[, c("gene_id", "gene_name", stat_cols)]

write.table(res_df, file = snakemake@output[["table"]],
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

cat("\nWrote", nrow(res_df), "genes to", snakemake@output[["table"]], "\n")

####
## PCA
####

# vst() needs enough genes to fit its parametric model (nsub = 1000 by
# default); fall back to the full-model transform on small test datasets.
vsd <- tryCatch(
  vst(dds, blind = FALSE),
  error = function(e) {
    message("vst() failed (", conditionMessage(e),
            "); falling back to varianceStabilizingTransformation()")
    varianceStabilizingTransformation(dds, blind = FALSE)
  }
)

pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(paste0(treatment, " vs ", reference)) +
  coord_fixed() +
  theme_bw()

ggsave(snakemake@output[["pca"]], plot = p, width = 7, height = 6)

# Same plot with each point annotated by sample name. show.legend = FALSE
# keeps the text glyph from replacing the point markers in the group legend.
p_labels <- p +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE)

ggsave(snakemake@output[["pca_labels"]], plot = p_labels, width = 7, height = 6)

cat("Wrote PCA to", snakemake@output[["pca"]], "\n")
cat("Wrote labelled PCA to", snakemake@output[["pca_labels"]], "\n")

sink(type = "message")
sink(type = "output")
close(log_con)
