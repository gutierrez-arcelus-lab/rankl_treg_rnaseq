library(tidyverse)
library(readxl)

meta <- 
    read_excel("./data/Broad_SK-691I_RLB.xls", sheet = "Sample Info", skip = 2, col_names = FALSE) |>
    select(well = 1, sample_name = 4, sex = 6) |>
    mutate(well = sub("^([A-Z])0*([0-9]+)$", "\\1\\2", well))  

fastq_list <- 
    read_csv("./data/fastq/Reports/fastq_list.csv")

fastp_meta <- 
    fastq_list |>
    extract(RGSM, "well", "_([A-Z][0-9]+)$") |>
    select(well, lane = Lane, read1 = Read1File, read2 = Read2File) |>
    mutate_at(vars(read1, read2), ~file.path("./data/fastq", basename(.))) |>
    arrange(well, lane) |>
    group_by(well) |>
    summarise_at(vars(read1, read2),
		 ~paste(., collapse = ",")) |>
    ungroup() |>
    left_join(meta, join_by(well)) |>
    select(well, sample_name, sex, read1, read2) |>
    mutate(well = factor(well, levels = unique(str_sort(well, numeric = TRUE)))) |>
    arrange(well)

write_tsv(fastp_meta, "./data/metadata_fastp.tsv")
