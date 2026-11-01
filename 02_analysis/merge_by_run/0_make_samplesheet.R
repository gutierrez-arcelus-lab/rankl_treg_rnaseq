library(tidyverse)

mqc <- 
    "../../00_fastq_qc/output/multiqc/multiqc_data/multiqc_general_stats.txt" |>
    read_tsv()

reads <-
    mqc |>
    select(Sample, passed = matches("passed_filter_reads")) |>
    mutate(run = str_remove(Sample, "_R1$"))

meta <-
    read_tsv("../../00_fastq_qc/data/metadata_fastp.tsv") |>
    mutate(run = str_c(well, sample_name, sep = "_"))

# one row per library (run) passing QC
samples_run <-
    meta |>
    inner_join(reads, by = "run") |>
    filter(passed >= 2e6) |>
    select(run, sample_name)

write_tsv(samples_run, "data/metadata_salmon_runs.tsv")
