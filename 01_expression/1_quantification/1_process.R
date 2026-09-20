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
    "../../00_fastq_qc/data/metadata_fastp.tsv" |>
    read_tsv() |>
    select(1:3) |>
    unite("run", c(well, sample_name), sep = "_", remove = FALSE)

mapping <- 
    list.dirs("output/salmon", recursive = FALSE) |>
    set_names(basename) |>
    map_dfr(read_map_rate, .id = "run") |>
    mutate(percent_decoy    = (num_decoy  / num_processed) * 100,
           percent_unmapped = (1 - (num_mapped + num_decoy) / num_processed) * 100)

mqc <- 
    read_tsv("../../00_fastq_qc/output/multiqc/multiqc_data/multiqc_general_stats.txt") |>
    select(run = Sample, 
	   pct_passed_filter_reads = "fastp_mqc-generalstats-fastp-pct_surviving") |>
    mutate(run = str_remove(run, "_R1$"))
	   
merge_data <- 
    left_join(mapping, mqc, join_by(run)) |>
    left_join(meta, join_by(run)) |>
    mutate(experiment = case_when(
				  str_detect(sample_name, "^(Tregs|Teff)-")     ~ "Experiment 1",
				  str_detect(sample_name, "^iYFP-")             ~ "Experiment 2",
				  str_detect(sample_name, "_IL-6$|_IL-1b$")     ~ "Experiment 3",
				  str_detect(sample_name, "_RANKL$|_RANKLneg$") ~ "Experiment 4",
				  str_detect(sample_name, "^Thymus-")           ~ "Experiment 5",
				  TRUE ~ NA_character_)) |>
    select(experiment, run, well, sample_name, sex,
	   num_reads = num_processed,
	   percent_passed_filter = pct_passed_filter_reads, 
	   percent_mapped, percent_decoy, percent_unmapped) |>
    mutate(well = factor(well, levels = unique(str_sort(well, numeric = TRUE)))) |>
    arrange(well)

write_tsv(merge_data, "./output/mapping_metrics.tsv")


# Salmon quant files
files <-
    file.path("./output/salmon", meta$run, "quant.sf") |>
    set_names(meta$run)

tx2gene <- read_tsv("../0_index/data/reference/tx2gene.tsv")

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

write_rds(txi, "./output/salmon_txi.rds")


# Plot
library(ggh4x)

tmp_data_1 <- 
    merge_data |>
    select(experiment, run, well, sample_name,
	   `Reads after filter (millions)` = num_reads, `% passed filter` = percent_passed_filter) |>
    mutate(`Reads after filter (millions)` = `Reads after filter (millions)` / 1e6) |>
    pivot_longer(c(`Reads after filter (millions)`, `% passed filter`), names_to = "statistic") |>
    mutate(category = "NA_character_")

tmp_data_2 <- 
    merge_data |>
    select(experiment, run, well, sample_name,
	   mapped = percent_mapped, decoy = percent_decoy, unmapped = percent_unmapped) |>
    pivot_longer(c(mapped, decoy, unmapped), names_to = "category") |>
    mutate(statistic = "Mapping (%)")

plot_df <- 
    bind_rows(tmp_data_1, tmp_data_2) |>
    mutate(statistic = factor(statistic,
			      levels = c("Reads after filter (millions)", "% passed filter", "Mapping (%)")),
	   category  = factor(category, levels = c("mapped", "decoy", "unmapped"))) |>
    arrange(experiment, well) |>
    mutate(sample_name = fct_inorder(sample_name))

# pin the percentage columns to 0–100 (reads stays free), via invisible anchors
blank <- 
    plot_df |>
    filter(statistic %in% c("% passed filter", "Mapping (%)")) |>
    distinct(experiment, sample_name, well, statistic) |>
    crossing(value = c(0, 100))

thr <- 
    tibble(statistic = factor("Reads after filter (millions)", levels = levels(plot_df$statistic)),
	   x = 1)   

p <- 
    ggplot(plot_df, aes(x = value, y = well, fill = category)) +
    geom_col(position = position_stack(reverse = TRUE)) +
    geom_vline(data = thr, aes(xintercept = x),
               linetype = "dashed", linewidth = 0.3, colour = "grey30") +
    geom_blank(data = blank, aes(x = value, y = well), inherit.aes = FALSE) +
    facet_nested(experiment + sample_name ~ statistic, scales = "free") +
    scale_fill_manual(values = c(mapped = "#4daf4a", decoy = "#ff7f00", unmapped = "grey75"),
                      na.value = "grey50", breaks = c("mapped", "decoy", "unmapped"),
                      name = NULL) +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = "Well") +
    theme_bw(base_size = 7) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(angle = 0),
          strip.background = element_rect(fill = "grey95", colour = NA),
          legend.position = "top",
          panel.spacing.x = unit(4, "pt"))

ggsave("./plots/qc_dashboard.pdf", p, width = 8.5, height = 11, units = "in")
