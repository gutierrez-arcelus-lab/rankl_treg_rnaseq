library(tidyverse)
library(edgeR)
library(ggsignif)   # for the comparison bracket
library(ggrepel)
library(patchwork)
library(gridExtra)
library(fgsea)
library(msigdbr)

# Expression data
y <- read_rds("../02_dge/output/edger/dgelist_collapsed.rds")

gene_meta <- 
    "../01_expression/0_index/data/reference/gencode.vM39.gene_metadata.tsv" |>
    read_tsv()

# Figure 3-h,i: RANKL+ vs RANKL- -----------------------------------------------
# Exp4: RANKL+ vs RANKL- (paired within mouse) ---------------------------------
s4_0     <- y[, y$samples$experiment == "exp4"]
mouse4 <- factor(s4_0$samples$mouse)
cond4  <- factor(s4_0$samples$condition, levels = c("RANKLneg", "RANKL"))
design4 <- model.matrix(~ mouse4 + cond4)                             

s4   <- s4_0[filterByExpr(s4_0, group = cond4), , keep.lib.sizes = FALSE]
s4   <- normLibSizes(s4)
s4   <- estimateDisp(s4, design4, robust = TRUE)

res4 <- read_tsv("../02_dge/output/edger/Exp4_RANKLpos_vs_RANKLneg.tsv.gz")

# 1. top 50 genes by p-value
top_genes <- 
    res4 |>
    slice_min(PValue, n = 50) |>
    arrange(sign(logFC), desc(PValue))

# 2. expression matrix (logCPM)
mat <- cpm(s4, log = TRUE, prior.count = 2)[top_genes$gene_id, ] 

# 3. row-scale: z-score per gene (scale() works on columns, hence the transposes)
mat_z <- t(scale(t(mat)))

# 4. long format + attach sample metadata
long <- 
    mat_z |>
    as.data.frame() |>
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "sample_name", values_to = "z") |>
    left_join(s4$sample, join_by(sample_name)) |>
    left_join(s4$genes, join_by(gene_id)) |>
    select(mouse, sample_name, sex, genotype, condition, gene_id, gene_name, gene_type, z) |>
    mutate(mouse = paste0("M", mouse),
	   condition = recode(condition, "RANKL" = "RANKL+", "RANKLneg" = "RANKL-"),
	   condition = factor(condition, levels = c("RANKL-", "RANKL+")),
	   gene_name = factor(gene_name, levels = top_genes$gene_name),
	   mouse = str_remove(mouse, "exp4_"))

fig3_h <- 
    ggplot(long, aes(x = mouse, y = gene_name, fill = z)) +
    geom_tile() +
    facet_grid(~ condition, scales = "free_x") +
    scico::scale_fill_scico(palette = "vik", 
			    midpoint = 0,
			    limit = c(-2, 2),
			    breaks = c(-2, -1, 0, 1, 2)
			    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(
	  axis.text.x.top = element_text(size = 10, margin = margin(t = 0, b = 0)),
          axis.text.y = element_text(size = 9, face = "italic"), 
	  panel.grid = element_blank(),
	  #plot.title = element_text(size = 10, face = "bold"),
	  strip.text = element_text(size = 10, margin = margin(t = 0, b = 0)),
	  strip.placement = "outside"
	  ) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colorbar(barwidth = .5))

ggsave("pdf/fig3_h.pdf", fig3_h, 
       width = 120, height = 160, 
       dpi = 600,
       units = "mm", 
       device = cairo_pdf)

# Enrichment analysis
ranks <- 
    res4 |>
    mutate(tstat = sign(logFC) * sqrt(F)) |> 
    group_by(gene_name) |> 
    slice_max(abs(tstat), n = 1) |> 
    ungroup() |>  
    select(gene_name, tstat) |> 
    deframe()

pathways <- 
    msigdbr(species = "Mus musculus", category = "H") |>  
    (\(x) split(x = x$gene_symbol, f = x$gs_name))()

set.seed(1)
gsea <- 
    fgsea(pathways, ranks, minSize = 15, maxSize = 500) |>
    as_tibble() |> 
    arrange(padj)

plot_data_gsea <- 
    gsea |>
    filter(padj < 0.05) |>
    top_n(n = 16, wt = -log10(pval)) |>
    mutate(pathway = str_remove(pathway, "^HALLMARK_") |> str_replace_all("_", " "),
           pathway = fct_reorder(pathway, NES),
           le  = lengths(leadingEdge),
           dir = if_else(NES > 0, "Higher in RANKL+", "Higher in RANKL−"))


fig3_i <- 
    ggplot(plot_data_gsea, aes(NES, pathway)) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_segment(aes(x = 0, xend = NES, yend = pathway), color = "grey75") +
    geom_point(aes(size = le, color = dir)) +
    scale_x_continuous(limit = c(-3, 3.5), 
		       breaks = c(-3, -2, -1, 0, 1, 2, 3)
		       ) +
    scale_color_manual(values = c("Higher in RANKL+" = "#B2182B",
				  "Higher in RANKL−" = "#2166AC")) +
    scale_size_area(max_size = 4, name = "Leading-edge\ngenes") +
    theme_minimal(base_size = 10) +
    theme(
	  panel.grid.minor.x = element_blank(),
	  panel.grid.major = element_line(color = "gray96"),
	  legend.text = element_text(size = 8),
	  plot.margin = margin(t = 2, b = 2, unit = "pt")
	  ) +   
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(x = "Normalized enrichment score", 
	 y = NULL, color = NULL) 


ggsave("pdf/fig3_i.pdf", fig3_i, 
       width = 120, height = 66, 
       dpi = 600,
       units = "mm", 
       device = cairo_pdf)

# Figure S3-g: Experiment 4 selected genes--------------------------------------
selected_genes <-
    c("Tnfsf11", "Runx2", "Thop1", "Mlf1", "9130401M01Rik", "Idi2", "Dgat2",
      "Dctpp1", "Ccr6", "Gm57322", "Ckb", "Ccl20", "Scin", "Penk", "Nfkb2",
      "Tnfsf9", "Il18", "Ifng", "Rora", "Pdcd1", "Ccl5", "Tgfb2", "Gata3",
      "Il1rl1", "Sell", "Tlr1", "Il6ra", "Usp18", "Foxp3", "Socs3", "Ifi27l2a",
      "Usp3", "Ifi203", "Acss1", "Rapgef6", "Irf4", "Il1rn", "Il2ra", "Areg",
      "Notch2", "Foxp1", "Aire", "Il10", "Icos", "Foxp2", "Cxcr5", "Il2",
      "Cd69", "Cd40lg", "Ctla4", "Itgae", "Cd28", "Foxp4", "Nfkb1", "Tgfb1",
      "Tgfb3")

sel_genes <-
    filter(gene_meta, gene_name %in% selected_genes)

mat_sel <- 
    cpm(s4_0, log = TRUE, prior.count = 2) |>
    {\(x) x[rownames(x) %in% sel_genes$gene_id, ]}()

# 3. row-scale: z-score per gene (scale() works on columns, hence the transposes)
mat_sel_z <- t(scale(t(mat_sel)))

# 4. long format + attach sample metadata
long_sel <- 
    mat_sel_z |>
    as.data.frame() |>
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "sample_name", values_to = "z") |>
    left_join(s4_0$sample, join_by(sample_name)) |>
    left_join(s4_0$genes, join_by(gene_id)) |>
    select(mouse, sample_name, sex, genotype, condition, gene_id, gene_name, gene_type, z) |>
    mutate(mouse = paste0("M", mouse),
	   condition = recode(condition, "RANKL" = "RANKL+", "RANKLneg" = "RANKL-"),
	   condition = factor(condition, levels = c("RANKL-", "RANKL+")),
	   gene_name = factor(gene_name, levels = rev(selected_genes)),
	   mouse = str_remove(mouse, "exp4_"))

fig_s3g <- 
    ggplot(long_sel, aes(x = mouse, y = gene_name, fill = z)) +
    geom_tile() +
    facet_grid(~ condition, scales = "free_x") +
    scico::scale_fill_scico(palette = "vik", 
			    midpoint = 0,
			    limit = c(-2, 2), oob = scales::squish,
			    breaks = c(-2, -1, 0, 1, 2)
			    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(
	  axis.text.x.top = element_text(size = 10, margin = margin(t = 0, b = 0)),
          axis.text.y = element_text(size = 9, face = "italic"), 
	  panel.grid = element_blank(),
	  strip.text = element_text(size = 10, margin = margin(t = 0, b = 0)),
	  strip.placement = "outside"
	  ) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colorbar(barwidth = .5))

ggsave("pdf/Sfig3_g.pdf", fig_s3g, 
       width = 120, height = 160, 
       dpi = 600,
       units = "mm", 
       device = cairo_pdf)
