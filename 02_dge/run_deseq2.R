library(tidyverse)
library(readxl)
library(tximport)
library(DESeq2)

ref  <- "../01_expression/0_index/data/reference"
outd <- "output/deseq2"

dir.create(outd, showWarnings = FALSE, recursive = TRUE)

# metadata ---------------------------------------------------------------------
mut <-
    "../00_fastq_qc/data/Broad_SK-691I_RLB for Vitor.xls" |>
    read_excel(sheet = "Sample Info", skip = 2, col_names = FALSE) |>
    select(sample_name = 3, genotype = 5, rin = 9) |>
    filter(!is.na(genotype))

metrics <-
    "../01_expression/1_quantification/output/mapping_metrics.tsv" |>
    read_tsv() |>
    select(experiment, run) |>
    mutate(experiment = str_replace(experiment, "Experiment ", "exp"))

meta <-
    read_tsv("../00_fastq_qc/data/metadata_fastp.tsv") |>
    select(1:3) |>
    mutate(run = str_c(well, sample_name, sep = "_")) |>
    left_join(mut, join_by(sample_name)) |>
    left_join(metrics, join_by(run)) |>
    mutate(rep_half = if_else(str_detect(well, "^[A-D]"), "rep1", "rep2"),
	   mouse = str_c(experiment, str_extract(sample_name, "[0-9]+"), sep = "_"),
	   condition = case_when(
				 str_detect(sample_name, "^Tregs-")    ~ "Treg",
				 str_detect(sample_name, "^Teff-")     ~ "Teff",
				 str_detect(sample_name, "_IL-6$")     ~ "IL-6",
				 str_detect(sample_name, "_IL-1b$")    ~ "IL-1b",
				 str_detect(sample_name, "_RANKLneg$") ~ "RANKLneg",
				 str_detect(sample_name, "_RANKL$")    ~ "RANKL",
				 str_detect(sample_name, "^iYFP-")     ~ "iTregs",
				 str_detect(sample_name, "^Thymus-")   ~ "Thymus"))

# run-level depth filter -------------------------------------------------------
read_meta_info <- function(run) {

    f <- file.path("../01_expression/1_quantification/output/salmon", run, "aux_info/meta_info.json")

    j <- jsonlite::fromJSON(f)

    tibble(run = run, num_processed = j$num_processed, percent_mapped = j$percent_mapped)
}

qc <- map_dfr(meta$run, read_meta_info)

meta_filt <-
    meta |>
    left_join(qc, by = "run") |>
    filter(!is.na(num_processed), num_processed >= 1e6)

# --------------------------------------------- import + collapse tech reps ---
# Same input as the edgeR pipeline: lengthScaledTPM counts (length correction
# baked in, so no length offset is used here either) - this keeps DESeq2 and
# edgeR maximally comparable. DESeq2 then applies its own median-of-ratios size
# factors on top. (The DESeq2-native alternative is DESeqDataSetFromTximport with
# countsFromAbundance="no", which builds a per-sample length offset instead; it is
# equally valid but would make the two pipelines differ in normalisation.)
files <-
    file.path("../01_expression/1_quantification/output/salmon", meta_filt$run, "quant.sf") |>
    set_names(meta_filt$run)

tx2gene <-
    file.path(ref, "tx2gene.tsv") |>
    read_tsv()

txi <-
    tximport(files, type = "salmon", tx2gene = tx2gene,
	     countsFromAbundance = "lengthScaledTPM")

samples <-
    meta_filt |>
    as.data.frame() |>
    column_to_rownames("run")

samples <- samples[colnames(txi$counts), ]
identical(rownames(samples), colnames(txi$counts))

# Collapse the two library reps per biological sample. collapseReplicates() sums
# the count columns sharing an ID and reduces colData to one row per sample - the
# idiomatic DESeq2 way. It needs an integer dds first, so counts are rounded before
# summing; vs rounding after the sum this changes results by well under one count
# per gene, i.e. not at all. (edgeR did not require integers, so this rounding is
# the only numerical difference between the two pipelines' inputs.)
dds_all <-
    DESeqDataSetFromMatrix(round(txi$counts), colData = samples, design = ~ 1) |>
    collapseReplicates(groupby = samples$sample_name)

cts     <- counts(dds_all)
coldata <- as.data.frame(colData(dds_all))

identical(rownames(coldata), colnames(cts))

gene_meta <-
    read_tsv(file.path(ref, "gencode.vM39.gene_metadata.tsv"))

saveRDS(list(counts = cts, coldata = coldata), file.path(outd, "counts_collapsed.rds"))

# Small DESeq2-standard pre-filter: keep genes with >= 10 counts in at least as
# many samples as the smallest group. Independent filtering in results() handles
# the rest. 
prefilter <- function(dds, smallest_group) {
    keep <- rowSums(counts(dds) >= 10) >= smallest_group
    dds[keep, ]
}


# results() by the LAST coefficient in resultsNames() - the group term is written
# last in every design, exactly like coef = ncol(design) in the edgeR script, so
# it is robust to DESeq2 sanitising level names like "IL-1b" -> "IL.1b".
tidy_res <- function(dds) {
    results(dds, name = tail(resultsNames(dds), 1), alpha = 0.05) |>
        as.data.frame() |>
        rownames_to_column("gene_id") |>
        as_tibble() |>
        left_join(gene_meta, by = "gene_id") |>
        arrange(pvalue)
}

# ==============================================================================
# Differential expression
# ==============================================================================

# Exp3: IL-1b vs IL-6  (paired within mouse) -----------------------------------
cd3 <- coldata[coldata$experiment == "exp3", ]
cd3$mouse         <- factor(cd3$mouse)
cd3$condition_grp <- factor(cd3$condition, levels = c("IL-6", "IL-1b"))
m3  <- cts[, rownames(cd3)]

dds3 <- DESeqDataSetFromMatrix(m3, cd3, design = ~ mouse + condition_grp)
dds3 <- prefilter(dds3, min(table(cd3$condition_grp)))
dds3 <- DESeq(dds3)
res3 <- tidy_res(dds3)

write_tsv(res3, file.path(outd, "Exp3_IL1b_vs_IL6.tsv.gz"))

# Exp4: RANKL+ vs RANKL-  (paired within mouse) --------------------------------
cd4 <- coldata[coldata$experiment == "exp4", ]
cd4$mouse         <- factor(cd4$mouse)
cd4$condition_grp <- factor(cd4$condition, levels = c("RANKLneg", "RANKL"))
m4  <- cts[, rownames(cd4)]

dds4 <- DESeqDataSetFromMatrix(m4, cd4, design = ~ mouse + condition_grp)
dds4 <- prefilter(dds4, min(table(cd4$condition_grp)))
dds4 <- DESeq(dds4)
res4 <- tidy_res(dds4)

write_tsv(res4, file.path(outd, "Exp4_RANKLpos_vs_RANKLneg.tsv.gz"))

# Exp5: Thymus dRANKL vs WT  (unpaired; adjust for sex if estimable) -----------
cd5 <- coldata[coldata$experiment == "exp5", ]
cd5$sex          <- factor(cd5$sex)
cd5$genotype_grp <- factor(cd5$genotype, levels = c("WT", "Mut"))
m5  <- cts[, rownames(cd5)]

design5 <- ~ sex + genotype_grp

dds5 <- DESeqDataSetFromMatrix(m5, cd5, design = design5)
dds5 <- prefilter(dds5, min(table(cd5$genotype_grp)))
dds5 <- DESeq(dds5)
res5 <- tidy_res(dds5)  

write_tsv(res5, file.path(outd, "Exp5_Thymus_dRANKL_vs_WT.tsv.gz"))

# ---- quick summary ---------------------------------------------------------
bind_rows(
    tibble(contrast = "Exp3_IL1b_vs_IL6",          res = list(res3)),
    tibble(contrast = "Exp4_RANKLpos_vs_RANKLneg", res = list(res4)),
    tibble(contrast = "Exp5_Thymus_dRANKL_vs_WT",  res = list(res5))) |>
    mutate(tested = map_int(res, nrow),
           up   = map_int(res, ~ sum(.x$padj < 0.05 & .x$log2FoldChange > 0, na.rm = TRUE)),
           down = map_int(res, ~ sum(.x$padj < 0.05 & .x$log2FoldChange < 0, na.rm = TRUE)),
           res = NULL) |>
    print()
