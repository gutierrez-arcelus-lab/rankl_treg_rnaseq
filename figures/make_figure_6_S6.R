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

# Figure 6-f: Tnfsf11 expression -----------------------------------------------

# ==============================================================================
# Exp2: iYFP Foxp3^WT vs Foxp3^dRANKL
# ==============================================================================
s2    <- y[, y$samples$experiment == "exp2"]
sex2  <- factor(s2$samples$sex)
geno2 <- factor(s2$samples$genotype, levels = c("WT", "Mut"))

design2 <- model.matrix(~ sex2 + geno2)

s2   <- s2[filterByExpr(s2, group = geno2), , keep.lib.sizes = FALSE]
s2   <- normLibSizes(s2)
s2   <- estimateDisp(s2, design2, robust = TRUE)

res2 <- read_tsv("../02_dge/output/edger/Exp2_iYFP_Mut_vs_WT.tsv.gz")

# ==============================================================================
# Exp3: IL-1b vs IL-6
# ==============================================================================
s3_0   <- y[, y$samples$experiment == "exp3"]
mouse3 <- factor(s3_0$samples$mouse)
cond3  <- factor(s3_0$samples$condition, levels = c("IL-6", "IL-1b"))
design3 <- model.matrix(~ mouse3 + cond3)

s3   <- s3_0[filterByExpr(s3_0, group = cond3), , keep.lib.sizes = FALSE]
s3   <- normLibSizes(s3)
s3   <- estimateDisp(s3, design3, robust = TRUE)

res3 <- read_tsv("../02_dge/output/edger/Exp3_IL1b_vs_IL6.tsv.gz")

# ==============================================================================
# Figure
# ==============================================================================

# CPM on each experiment's own TMM normalization
gene_cpm <- function(s, gene) {
    i <- which(s$genes$gene_name == gene)
    tibble(sample = colnames(s), CPM = cpm(s)[i, ]) |>
        left_join(as_tibble(s$samples, rownames = "sample"), by = "sample")
}

d6f <-
    bind_rows(gene_cpm(s2, "Tnfsf11") |> mutate(group = as.character(genotype)),
              gene_cpm(s3, "Tnfsf11") |> mutate(group = condition)) |>
    mutate(group = factor(group, levels = c("WT", "Mut", "IL-6", "IL-1b")))

# p-values from the genome-wide fits 
pv <- 
    c(filter(res2, gene_name == "Tnfsf11")$PValue,
      filter(res3, gene_name == "Tnfsf11")$PValue) |>
    {\(p) case_when(p < 0.001 ~ as.character(signif(p, 2)), TRUE ~ sprintf("%.2f", p))}()

# one bracket per facet
ann <-
    tibble(experiment = c("exp2", "exp3"),
           group1     = c("WT",   "IL-6"),
           group2     = c("Mut",  "IL-1b"),
           label      = pv) |>
    left_join(d6f |>
                  group_by(experiment) |>
                  summarise(y_position = max(CPM) * 1.08),
              join_by(experiment))


xlab_6f <- \(x) parse(text = c("WT"    = "WT",
                               "Mut"   = "Delta*RANKL",
                               "IL-6"  = "'IL-6'",
                               "IL-1b" = "'IL-1'*beta")[x])

fig_6f <-
    ggplot(d6f, aes(group, CPM)) +
    geom_point(aes(fill = group), 
	       position = position_jitter(width = 0.08, height = 0, seed = 1),
               size = 2.5, shape = 21, stroke = .25) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.4) +
    geom_signif(data = ann,
                aes(xmin = group1, xmax = group2,
                    annotations = label, y_position = max(ann$y_position)),
                manual = TRUE, tip_length = 0.02, textsize = 3, vjust = -0.3) +
    facet_wrap(~ experiment, scales = "free",
               labeller = as_labeller(c(exp2 = "Foxp3 genotype",
                                        exp3 = "Cytokine"))) +
    scale_x_discrete(labels = xlab_6f) +
    scale_y_continuous(limits = c(0, 25), expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = c("WT" = "#BEBEBE", "Mut" = "#7D2C4F", "IL-6" = "#F2D791", "IL-1b" = "#548D91"),
		      guide = "none") +
    labs(x = NULL, y = "Expression (CPM)",
	 title = "Tnfsf11"
	 ) +
    theme_classic(base_size = 10) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
	  axis.text = element_text(face = "bold"),
	  axis.title = element_text(face = "bold"),
	  plot.title = element_text(size = 10, hjust = 0.5, face = "italic")
	  )

ggsave("pdf/fig6_f.pdf", fig_6f, 
       width = 60, height = 60, 
       dpi = 600,
       units = "mm", 
       device = cairo_pdf)





# Figure S6-c: Experiment 3 ----------------------------------------------------

# top 20 genes by p-value
top_genes_exp3 <- 
    res3 |>
    slice_min(PValue, n = 30) |>
    arrange(sign(logFC), desc(PValue))

# expression matrix (logCPM)
mat_exp3 <- cpm(s3, log = TRUE, prior.count = 2)[top_genes_exp3$gene_id, ] 

# row-scale: z-score per gene (scale() works on columns, hence the transposes)
mat_exp3_z <- t(scale(t(mat_exp3)))

# long format + attach sample metadata
long_exp3 <- 
    mat_exp3_z |>
    as.data.frame() |>
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "sample_name", values_to = "z") |>
    left_join(s3$sample, join_by(sample_name)) |>
    left_join(s3$genes, join_by(gene_id)) |>
    select(mouse, sample_name, sex, genotype, condition, gene_id, gene_name, gene_type, z) |>
    mutate(mouse = paste0("M", mouse),
	   condition = factor(condition, levels = c("IL-1b", "IL-6")),
	   gene_name = factor(gene_name, levels = top_genes_exp3$gene_name),
	   mouse = str_remove(mouse, "exp3_"))

fig_s6c <- 
    ggplot(long_exp3, aes(x = mouse, y = gene_name, fill = z)) +
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
	  strip.text = element_text(size = 10, margin = margin(t = 0, b = 0)),
	  strip.placement = "outside"
	  ) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colorbar(barwidth = .5))

ggsave("pdf/Sfig6_c.pdf", fig_s6c, 
       width = 120, height = 100, 
       dpi = 600,
       units = "mm", 
       device = cairo_pdf)

# Figure S6-d: Experiment 3 ----------------------------------------------------

# selected genes
selected_genes_exp3 <- 
    c("Tgfb3", "Il18", "Il2ra", "Cxcr5", "Cd40lg", "Ctla4", "Cd28", "Il10",
      "Il1rn", "Icos", "Il2", "Nfkb1", "Nfkb2", "Ifng", "Tnfsf11", "Foxp4",
      "Notch2", "Aire", "Tgfb1", "Ccr6", "Foxp1", "Foxp3", "Tgfb2", "Foxp2",
      "Runx2")

sel_genes_exp3 <- 
    gene_meta |>
    filter(gene_name %in% selected_genes_exp3)

# expression matrix (logCPM)
mat_exp3_sel <- cpm(s3_0, log = TRUE, prior.count = 2)[sel_genes_exp3$gene_id, ] 

# row-scale: z-score per gene (scale() works on columns, hence the transposes)
mat_exp3_sel_z <- t(scale(t(mat_exp3_sel)))

# long format + attach sample metadata
long_exp3_sel <- 
    mat_exp3_sel_z |>
    as.data.frame() |>
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "sample_name", values_to = "z") |>
    left_join(s3_0$sample, join_by(sample_name)) |>
    left_join(s3_0$genes, join_by(gene_id)) |>
    select(mouse, sample_name, sex, genotype, condition, gene_id, gene_name, gene_type, z) |>
    mutate(mouse = paste0("M", mouse),
	   condition = factor(condition, levels = c("IL-1b", "IL-6")),
	   gene_name = factor(gene_name, levels = rev(selected_genes_exp3)),
	   mouse = str_remove(mouse, "exp3_"))

fig_s6d <- 
    ggplot(long_exp3_sel, aes(x = mouse, y = gene_name, fill = z)) +
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
	  strip.text = element_text(size = 10, margin = margin(t = 0, b = 0)),
	  strip.placement = "outside"
	  ) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colorbar(barwidth = .5))

ggsave("pdf/Sfig6_d.pdf", fig_s6d, 
       width = 120, height = 100, 
       dpi = 600,
       units = "mm", 
       device = cairo_pdf)



# Bonus: Pathway enrichment for IL1-b vs IL-6

ranks_exp3 <- 
    res3 |>
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
gsea_exp3 <- 
    fgsea(pathways, ranks_exp3, minSize = 15, maxSize = 500) |>
    as_tibble() |> 
    arrange(padj)

plot_data_gsea_exp3 <- 
    gsea_exp3 |>
    filter(padj < 0.05) |>
    top_n(n = 25, wt = -log10(pval)) |>
    mutate(pathway = str_remove(pathway, "^HALLMARK_") |> str_replace_all("_", " "),
           pathway = fct_reorder(pathway, NES),
           le  = lengths(leadingEdge),
           dir = if_else(NES > 0, "Higher in IL-1β", "Higher in IL-6"))


gsea_exp3_plot <- 
    ggplot(plot_data_gsea_exp3, aes(NES, pathway)) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_segment(aes(x = 0, xend = NES, yend = pathway), color = "grey75") +
    geom_point(aes(size = le, color = dir)) +
    scale_x_continuous(limit = c(-3, 3.5), 
		       breaks = c(-3, -2, -1, 0, 1, 2, 2)
		       ) +
    scale_color_manual(values = c("Higher in IL-1β" = "#B2182B",
				  "Higher in IL-6" = "#2166AC")) +
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


ggsave("pdf/gsea_IL1BvsIL6.pdf", gsea_exp3_plot, 
       width = 120, height = 90, 
       dpi = 600,
       units = "mm", 
       device = cairo_pdf)


