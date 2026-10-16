library(tidyverse)
library(jsonlite)
library(tximport)

read_map_rate <- function(path) {
    f <- file.path(path, "aux_info/meta_info.json")
    if (!file.exists(f)) return(tibble(status = "missing"))
    m <- fromJSON(f)
    tibble(num_processed = m$num_processed,
           num_mapped    = m$num_mapped,
           num_decoy     = m$num_decoy_fragments,
           percent_mapped = m$percent_mapped,
           status = "ok")
}

meta <- 
    "../../../00_fastq_qc/data/metadata_fastp.tsv" |>
    read_tsv() |>
    distinct(sample_name, sex)

mapping <- 
    list.dirs("output/salmon", recursive = FALSE) |>
    set_names(basename) |>
    map_dfr(read_map_rate, .id = "sample_name") |>
    mutate(percent_decoy    = (num_decoy  / num_processed) * 100,
           percent_unmapped = (1 - (num_mapped + num_decoy) / num_processed) * 100)

merge_data <- 
    mapping |>
    left_join(meta, join_by(sample_name)) |>
    mutate(experiment = case_when(
				  str_detect(sample_name, "^(Tregs|Teff)-")     ~ "Experiment 1",
				  str_detect(sample_name, "^iYFP-")             ~ "Experiment 2",
				  str_detect(sample_name, "_IL-6$|_IL-1b$")     ~ "Experiment 3",
				  str_detect(sample_name, "_RANKL$|_RANKLneg$") ~ "Experiment 4",
				  str_detect(sample_name, "^Thymus-")           ~ "Experiment 5",
				  TRUE ~ NA_character_)) |>
    select(experiment, sample_name, sex,
	   num_reads = num_processed,
	   percent_mapped, percent_decoy, percent_unmapped) |>
    arrange(experiment)

write_tsv(merge_data, "./output/mapping_metrics.tsv")

# Salmon quant files
files <-
    file.path("./output/salmon", mapping$sample_name, "quant.sf") |>
    set_names(mapping$sample_name)

tx2gene <- read_tsv("../../../01_expression/0_index/data/reference/tx2gene.tsv")

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

write_rds(txi, "./output/salmon_txi.rds")
