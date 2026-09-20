library(tidyverse)
library(tximport)
library(DESeq2)
library(ggrepel)
library(glue)
library(readxl)

mut <- 
    "../../00_fastq_qc/data/Broad_SK-691I_RLB for Vitor.xls" |>
    read_excel(sheet = "Sample Info", skip = 2, col_names = FALSE) |>
    select(sample_name = 3, genot = 5, rin = 9) |>
    filter(!is.na(genot))

metrics <- 
    "../1_quantification/output/mapping_metrics.tsv" |>
    read_tsv() |>
    select(experiment, run, num_reads, percent_mapped)

meta <- 
    read_tsv("../../00_fastq_qc/data/metadata_fastp.tsv") |>
    select(-read1, -read2) |>
    mutate(
	   run = str_c(well, sample_name, sep = "_"),
	   condition = case_when(
				 str_detect(sample_name, "^Tregs-")    ~ "Treg",
				 str_detect(sample_name, "^Teff-")     ~ "Teff",
				 str_detect(sample_name, "_IL-6$")     ~ "iTregs_IL-6",
				 str_detect(sample_name, "_IL-1b$")    ~ "iTregs_IL-1b",
				 str_detect(sample_name, "_RANKLneg$") ~ "iTregs_RANKLneg",
				 str_detect(sample_name, "_RANKL$")    ~ "iTregs_RANKL",
				 str_detect(sample_name, "^iYFP-")     ~ "iTregs",
				 str_detect(sample_name, "^Thymus-")   ~ "Thymus",
				 TRUE ~ NA_character_),
	   rep_half = if_else(str_detect(well, "^[A-D]"), "rep1", "rep2")) |>
    left_join(mut, join_by(sample_name)) |>
    left_join(metrics, join_by(run))

# Salmon quant  
txi <- read_rds("../1_quantification/output/salmon_txi.rds") 

genes_detected <- 
    tibble(run = colnames(txi$counts),
	   genes_detected = colSums(txi$counts > 0))

# align colData to txi columns
sample_table <- 
    meta |>
    left_join(genes_detected, join_by(run)) |>
    mutate(run = factor(run, levels = colnames(txi$counts))) |>
    arrange(run) |>
    column_to_rownames("run")

# ---- DESeq2 + VST (blind, design ~1 for QC) ----
dds <- DESeqDataSetFromTximport(txi, colData = sample_table, design = ~1)
dds <- dds[rowSums(counts(dds)) > 0, ]

# ---- global PCA: sanity check that samples group by cell type / experiment ----
vsd <- vst(dds, blind = TRUE)

pca_data_all <- 
    plotPCA(vsd, 
	    intgroup = c("experiment"),
	    ntop = 1000, returnData = TRUE) |>
    as_tibble(rownames = "run") |>
    select(run, PC1, PC2) |>
    left_join(as_tibble(sample_table, rownames="run"), join_by(run))

pv_all <- round(100 * attr(pca_data_all, "percentVar"))

pca_global <- 
    ggplot(pca_data_all, 
	   aes(PC1, PC2, shape = experiment, fill = rin)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(data = pca_data_all |> filter(PC1 < 0 | num_reads < 1e6),
		    aes(label = run), 
		    size = 2, 
		    max.overlaps = Inf,
		    segment.size = 0.2, min.segment.length = 0, 
		    show.legend = FALSE) +
    scale_shape_manual(values = 21:25) +
    scale_fill_viridis_c(option = "viridis") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank()) +
    labs(x = glue("PC1 ({pv_all[1]}%)"), y = glue("PC2 ({pv_all[2]}%)"),
	 shape = "Experiment", fill = "RIN", 
	 title = "All samples -- by RIN")

ggsave("plots/pca_all_samples_rin.pdf", pca_global, width = 8, height = 6)

# ---- per-experiment PCA: outlier detection (re-VST within each) ----

plot_experiment <- function(dds, e) {

    dds_e <- dds[, dds$experiment == e]
    dds_e <- dds_e[rowSums(counts(dds_e)) > 0, ]
    vsd_e <- vst(dds_e, blind = TRUE)

    pca_data_e <- 
	plotPCA(vsd_e, 
		intgroup = c("experiment"),
		ntop = 500, returnData = TRUE) |>
	as_tibble(rownames = "run") |>
	select(run, PC1, PC2) |>
	left_join(meta, join_by(run)) |>
	left_join(genes_detected, join_by(run))

    pv_e <- round(100 * attr(pca_data_e, "percentVar"))

    pca_e <- 
	ggplot(pca_data_e, 
	       aes(PC1, PC2, shape = rep_half)) +
	geom_point(size = 3, alpha = 0.9) +
	geom_text_repel(
			aes(label = run), 
			size = 2, 
			max.overlaps = Inf,
			segment.size = 0.2, min.segment.length = 0, 
			show.legend = FALSE) +
	theme_bw(base_size = 11) +
	theme(panel.grid.minor = element_blank()) +
	labs(x = glue("PC1 ({pv_e[1]}%)"), y = glue("PC2 ({pv_e[2]}%)"),
	     shape = "Replicate", 
	     title = e)
}

pca_e1 <- plot_experiment(dds, "Experiment 1")
pca_e2 <- plot_experiment(dds, "Experiment 2")
pca_e3 <- plot_experiment(dds, "Experiment 3")
pca_e4 <- plot_experiment(dds, "Experiment 4")
pca_e5 <- plot_experiment(dds, "Experiment 5")

ggsave("plots/pca_e1.pdf", pca_e1, width = 8, height = 6)
ggsave("plots/pca_e2.pdf", pca_e2, width = 8, height = 6)
ggsave("plots/pca_e3.pdf", pca_e3, width = 8, height = 6)
ggsave("plots/pca_e4.pdf", pca_e4, width = 8, height = 6)
ggsave("plots/pca_e5.pdf", pca_e5, width = 8, height = 6)

library(broom)
library(matrixStats)

# --- N PCs, DESeq2-style (top-ntop variable genes, center, no scale) ---
get_pcs <- function(vsd, ntop = 500, n_pcs = 15) {
    mat <- assay(vsd)
    rv  <- rowVars(mat)
    sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, nrow(mat)))]
    pca <- prcomp(t(mat[sel, ]), center = TRUE, scale. = FALSE)
    k   <- min(n_pcs, ncol(pca$x))
    list(scores     = as_tibble(pca$x[, seq_len(k), drop = FALSE], rownames = "run"),
         percentVar = (pca$sdev^2 / sum(pca$sdev^2))[seq_len(k)])
}

covar_tbl <- as_tibble(sample_table, rownames = "run")

covars <- c("num_reads", "genes_detected", "percent_mapped", "rin",
            "experiment", "condition", "sex", "genot")

# --- PC ~ covariate R^2 + p (lm/glance); drops covariates that don't vary ---
pc_covariate_rsq <- function(scores) {
    
    dat <- left_join(scores, covar_tbl, by = "run")
    
    pc_cols <- setdiff(names(scores), "run")
    
    use <- covars[map_lgl(covars, ~ n_distinct(dat[[.x]][!is.na(dat[[.x]])]) > 1)]
    
    expand_grid(PC = pc_cols, Covariate = use) |>
	mutate(stats = map2(PC, Covariate,
			    ~ glance(lm(reformulate(.y, response = .x), data = dat)))) |>
	unnest(stats) |>
        transmute(PC, Covariate, r.squared, p.value)
}

# --- dotplot, generalized to whatever PCs are present ---
plot_rsq <- function(rsq, title = NULL) {
    
    lev <- str_sort(unique(rsq$PC), numeric = TRUE)
    
    rsq |>
        mutate(fdr = p.adjust(p.value, "fdr"),
               PC = factor(PC, levels = lev),
               neg_log_fdr = -log10(pmax(fdr, 1e-300))) |>
        ggplot(aes(PC, Covariate)) +
        geom_point(aes(size = neg_log_fdr, fill = r.squared), shape = 21, stroke = .2) +
        scale_fill_gradient(limits = c(0, 1), low = "white", high = "firebrick",
                            name = bquote(R^2)) +
        scale_size_continuous(name = expression(-log[10](FDR)), range = c(1, 8)) +
        labs(title = title, x = "Expression principal components", y = NULL) +
        theme_minimal(base_size = 11) +
        theme(panel.grid.major = element_line(color = "gray90"),
              panel.grid.minor = element_blank())
}

# --- global ---
g <- get_pcs(vsd, n_pcs = 15)
plot_rsq(pc_covariate_rsq(g$scores), "All samples")
ggsave("plots/pca_covariates_global.pdf", width = 8, height = 5)

# --- per experiment (PCs capped at n_samples - 1) ---
pdf("plots/pca_covariates_per_experiment.pdf", width = 8, height = 5)
walk(sort(unique(dds$experiment)), function(e) {
    sub <- dds[, dds$experiment == e]; sub <- sub[rowSums(counts(sub)) > 0, ]
    v   <- tryCatch(vst(sub, blind = TRUE),
                    error = function(z) varianceStabilizingTransformation(sub, blind = TRUE))
    pcs <- get_pcs(v, n_pcs = min(10, ncol(sub) - 1))
    print(plot_rsq(pc_covariate_rsq(pcs$scores), e))
})
dev.off()
