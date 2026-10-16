library(tidyverse)

mqc <- read_tsv("../../00_fastq_qc/output/multiqc/multiqc_data/multiqc_general_stats.txt")

reads <-
    mqc |>
    select(Sample, passed = matches("passed_filter_reads")) |>
    mutate(run = str_remove(Sample, "_R1$"))

trim_dir <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/colabs/raquel/rankl_treg_rnaseq/00_fastq_qc/output/fastq_trimmed"

meta <- 
    read_tsv("../../00_fastq_qc/data/metadata_fastp.tsv") |>
    mutate(run = str_c(well, sample_name, sep = "_"))

samples_merge <-
    meta |>
    inner_join(reads, by = "run") |>
    filter(passed >= 2e6) |>
    summarise(read1  = str_c(file.path(trim_dir, str_c(run, "_R1.fastq.gz")), collapse = ","),
              read2  = str_c(file.path(trim_dir, str_c(run, "_R2.fastq.gz")), collapse = ","),
	      .by = sample_name)

write_tsv(samples_merge, "data/metadata_salmon_merged.tsv")
