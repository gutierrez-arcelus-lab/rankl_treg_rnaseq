library(tidyverse)
library(readxl)
library(tximport)
library(edgeR)

ref  <- "../01_expression/0_index/data/reference"
outd <- "output/edger"

dir.create(outd, showWarnings = FALSE, recursive = TRUE)
dir.create("data", showWarnings = FALSE)

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

write_tsv(meta_filt, "./data/metadata_filtered.tsv")

# import + collapse tech reps -------------------------------------------------
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

y_run <- DGEList(counts = txi$counts, samples = samples)
y     <- sumTechReps(y_run, ID = y_run$samples$sample_name)

gene_meta <- 
    read_tsv(file.path(ref, "gencode.vM39.gene_metadata.tsv"),
	     show_col_types = FALSE) |>
    select(gene_id, gene_name, gene_type)

y$genes <- 
    tibble(gene_id = rownames(y)) |>
    left_join(gene_meta, by = "gene_id") |>
    as.data.frame()

write_rds(y, file.path(outd, "dgelist_collapsed.rds"))
#y <- read_rds(file.path(outd, "dgelist_collapsed.rds"))

# ==============================================================================
# Differential expression 
# ==============================================================================

# Exp2: iYFP Foxp3^WT vs Foxp3^dRANKL  (unpaired; adjust for sex) -------------
s2    <- y[, y$samples$experiment == "exp2"]
sex2  <- factor(s2$samples$sex)
geno2 <- factor(s2$samples$genotype, levels = c("WT", "Mut"))  
design2 <- model.matrix(~ sex2 + geno2) 

s2   <- s2[filterByExpr(s2, group = geno2), , keep.lib.sizes = FALSE]
s2   <- normLibSizes(s2)
s2   <- estimateDisp(s2, design2, robust = TRUE)
fit2 <- glmQLFit(s2, design2, robust = TRUE)

res2 <-
    glmQLFTest(fit2, coef = ncol(design2)) |>         
    topTags(n = Inf) |>
    (\(x) x$table)() |>
    as_tibble()

write_tsv(res2, file.path(outd, "Exp2_iYFP_Mut_vs_WT.tsv.gz"))


# Exp3: IL-1b vs IL-6  (paired within mouse) -----------------------------------
s3     <- y[, y$samples$experiment == "exp3"]
mouse3 <- factor(s3$samples$mouse)
cond3  <- factor(s3$samples$condition, levels = c("IL-6", "IL-1b"))       
design3 <- model.matrix(~ mouse3 + cond3)

s3   <- s3[filterByExpr(s3, group = cond3), , keep.lib.sizes = FALSE]
s3   <- normLibSizes(s3)
s3   <- estimateDisp(s3, design3, robust = TRUE)
fit3 <- glmQLFit(s3, design3, robust = TRUE)

res3 <- 
    glmQLFTest(fit3, coef = ncol(design3)) |>
    topTags(n = Inf) |>
    (\(x) x$table)() |>
    as_tibble()

write_tsv(res3, file.path(outd, "Exp3_IL1b_vs_IL6.tsv.gz"))

# Exp4: RANKL+ vs RANKL-  (paired within mouse) --------------------------------
s4     <- y[, y$samples$experiment == "exp4"]
mouse4 <- factor(s4$samples$mouse)
cond4  <- factor(s4$samples$condition, levels = c("RANKLneg", "RANKL"))
design4 <- model.matrix(~ mouse4 + cond4)                             

s4   <- s4[filterByExpr(s4, group = cond4), , keep.lib.sizes = FALSE]
s4   <- normLibSizes(s4)
s4   <- estimateDisp(s4, design4, robust = TRUE)
fit4 <- glmQLFit(s4, design4, robust = TRUE)

res4 <- 
    glmQLFTest(fit4, coef = ncol(design4)) |>
    topTags(n = Inf) |>
    (\(x) x$table)() |>
    as_tibble()

write_tsv(res4, file.path(outd, "Exp4_RANKLpos_vs_RANKLneg.tsv.gz"))

# Exp5: Thymus dRANKL vs WT  (unpaired; adjust for sex) -----------------------
s5    <- y[, y$samples$experiment == "exp5"]
sex5  <- factor(s5$samples$sex)
geno5 <- factor(s5$samples$genotype, levels = c("WT", "Mut")) 
design5 <- model.matrix(~ sex5 + geno5)

s5   <- s5[filterByExpr(s5, group = geno5), , keep.lib.sizes = FALSE]
s5   <- normLibSizes(s5)
s5   <- estimateDisp(s5, design5, robust = TRUE)
fit5 <- glmQLFit(s5, design5, robust = TRUE)

res5 <- 
    glmQLFTest(fit5, coef = ncol(design5)) |> 
    topTags(n = Inf) |>
    (\(x) x$table)() |>
    as_tibble()

write_tsv(res5, file.path(outd, "Exp5_Thymus_dRANKL_vs_WT.tsv.gz"))

# ---- quick summary ---------------------------------------------------------
bind_rows(
    tibble(contrast = "Exp3_IL1b_vs_IL6",          res = list(res3)),
    tibble(contrast = "Exp1_RANKLpos_vs_RANKLneg", res = list(res4)),
    tibble(contrast = "Exp5_Thymus_dRANKL_vs_WT",  res = list(res5))) |>
    mutate(tested = map_int(res, nrow),
           up     = map_int(res, ~ sum(.x$FDR < 0.05 & .x$logFC > 0)),
           down   = map_int(res, ~ sum(.x$FDR < 0.05 & .x$logFC < 0)),
           res = NULL) |>
    print()
