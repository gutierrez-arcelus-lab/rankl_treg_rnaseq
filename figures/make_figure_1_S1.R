library(tidyverse)
library(edgeR)
library(ggsignif)   # for the comparison bracket
library(ggrepel)
library(patchwork)
library(gridExtra)

# Expression data
y <- read_rds("../02_dge/output/edger/dgelist_collapsed.rds")

# ---- Figure 1-g: Aire, dRANKL vs WT ----------------------------------------
# Exp5: Thymus dRANKL vs WT  (unpaired; adjust for sex) -----------------------
s5    <- y[, y$samples$experiment == "exp5"]
sex5  <- factor(s5$samples$sex)
geno5 <- factor(s5$samples$genotype, levels = c("WT", "Mut")) 
design5 <- model.matrix(~ sex5 + geno5)

s5   <- s5[filterByExpr(s5, group = geno5), , keep.lib.sizes = FALSE]
s5   <- normLibSizes(s5)
s5   <- estimateDisp(s5, design5, robust = TRUE)

res5 <- read_tsv("../02_dge/output/edger/Exp5_Thymus_dRANKL_vs_WT.tsv.gz")

# dots = CPM on Exp5's own normalisation; 
# the p-value/FDR come from the same fit (res5)
aire <- which(s5$genes$gene_name == "Aire")

aire_expr <- 
    tibble(sample_name = colnames(s5),
	   CPM = cpm(s5, log = FALSE)[aire, ]) |>
    left_join(s5$samples, join_by(sample_name)) |>
    select(sample_name, sex, genotype, CPM) |>
    mutate(genotype = factor(genotype, levels = c("WT", "Mut")))

aire_stats <- 
    filter(res5, gene_name == "Aire") |>
    select(gene_name, logFC, logCPM, PValue, FDR)


fig1_g <- 
    ggplot(aire_expr, aes(x = genotype, y = CPM)) +
    stat_summary(fun = mean, geom = "col",
                 width = .6, fill = "white", color = "black", linewidth = .4) +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 width = .2, linewidth = .4) +
    geom_jitter(aes(fill = genotype), 
		shape = 21, stroke = .25,
		width = .15, size = 2) +
    geom_signif(
		comparisons = list(c("WT", "Mut")),
		annotations = round(aire_stats$PValue, 2),
		y_position = max(aire_expr$CPM) * 1.05,
		tip_length = 0, textsize = 8/.pt
    ) +
    scale_fill_manual(values = c("WT" = "#999999", "Mut" = "#75143F")) +
    scale_y_continuous(limits = c(0, NA),
		       expand = expansion(mult = c(0, .1))) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_blank(),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold"),
	  axis.ticks.x = element_blank(),
	  legend.position = "none",
	  plot.title = element_text(size = 10, hjust = .5, face = "plain")
	  ) +
    labs(x = NULL, y = "Expression (CPM)",
	 title = expression(italic("Aire")))
    
ggsave("pdf/fig1_g.pdf", fig1_g, 
       width = 30, height = 40, 
       dpi = 600,
       units = "mm", device = cairo_pdf)

# Figure S1-c: Experiment 5 volcano plot ---------------------------------------

## moderate LFCs
lfc_mod5 <- predFC(s5, design5, prior.count = 2)[, ncol(design5)]
res5 <- res5 |> mutate(logFC_mod = lfc_mod5[gene_id])


volcano_s1c <-
    ggplot(res5, aes(logFC_mod, -log10(PValue))) +
    geom_vline(xintercept = 0, color = "grey85", linetype = 2) +
    geom_point(aes(color = FDR < 0.05), 
	       size = 0.6, alpha = 0.6) +
    geom_text_repel(data =  top_n(res5, n = 20, wt = -log10(PValue)), 
		    aes(label = gene_name),
                    size = 8/.pt, fontface = "italic",
                    min.segment.length = 0, segment.size = 0.2,
                    max.overlaps = Inf) +
    scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "#B2182B"),
		       guide = "none") +
    theme_classic(base_size = 12) +
    theme(axis.title = element_text(face = "bold"),
	  axis.text = element_text(face = "bold")
	  ) +
    labs(x = expression(log[2]~"fold change (" * Delta * "RANKL vs WT)"),
         y = expression(-log[10]~"p-value"))


# Summary table
de_counts <- 
    res5 |>
    summarise(
	      across(c(PValue, FDR),
		     list(Total = ~ sum(.x < 0.05),
			  Up    = ~ sum(.x < 0.05 & logFC > 0),
			  Down  = ~ sum(.x < 0.05 & logFC < 0)))) |>
    pivot_longer(everything(),
                 names_to = c("Threshold", ".value"), names_sep = "_") |>
    mutate(Threshold = recode(Threshold,
                              PValue = "P < 0.05",
                              FDR    = "FDR < 5%"))

tbl_plot <- 
    tableGrob(de_counts, rows = NULL,
	      theme = ttheme_minimal(base_size = 12, 
				     colhead = list(
						    bg_params = list(fill = "#332288", col = NA),        
						    fg_params = list(col = "white", fontface = "bold"))
				     ))

fig_s1c <- 
    volcano_s1c / wrap_elements(tbl_plot) +
    plot_layout(heights = c(4, 1))

ggsave("pdf/Sfig1_c.pdf", fig_s1c, 
       width = 130, height = 130, 
       dpi = 600,
       units = "mm", 
       device = cairo_pdf)
