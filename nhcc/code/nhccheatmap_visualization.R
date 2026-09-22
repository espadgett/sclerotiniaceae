#last updated: Aug 17, 2026
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(ggridges)
library(forcats)
library(readxl)
library(tidyverse)
library(viridis)

#aim: graphically represent interchromosomal (trans-chromosome_ interactions from Streptotinia caulophylli/Sclerotinia minor HiC interaction data

#load and arrange Sminor interaction data
#all listed interactions have q-values < 0.05
#Smin_zscore <- read_delim("/Users/user/postgrad/lab_docx/fungal_asm_data/genome/Smin.Smin.Smin.trans1vsAll.10KB.trans_all.Smin.10000.zscores.zscore.mergeR.txt",
Smin_zscore <- read_delim("/Users/user/postgrad/lab_docx/fungal_asm_data/genome/Smin.Smin.Smin.trans1vsAll.10KB.trans_all.Smin.10000.zscores.zscore.mergeR.txt",
#Smin_zscore <- read_delim("/Users/user/postgrad/lab_docx/fungal_asm_data/genome/Smin.Smin.Smin.trans1vsAll.10KB.trans_all.Smin.10000.zscores.zscore.mergeR_qvalue.txt",
                          delim = "\t", escape_double = FALSE, 
                          col_names = TRUE, trim_ws = TRUE)
#View(Smin_zscore)
Smin_zscore$Sminor <- as.numeric(Smin_zscore$Smin)
Smin_zscore <- Smin_zscore %>%
  separate_wider_regex(
    ID,
    patterns = c(
      "Achr", CHR_A_NUMBER = "\\d+", "\\.",
      a_start = "\\d+", "\\.", a_end = "\\d+",
      "Bchr", CHR_B_NUMBER = "\\d+", "\\.",
      b_start = "\\d+", "\\.", b_end = "\\d+"
    )
  ) %>%
  mutate(
    CHR_A_REG = paste0(a_start, "-", a_end),
    CHR_B_REG = paste0(b_start, "-", b_end)
  ) %>%
  select(-a_start, -a_end, -b_start, -b_end)

Smin_zscore <- Smin_zscore %>%
  mutate(CHR_A_NUMBER = paste0("chr", CHR_A_NUMBER),
         CHR_B_NUMBER = paste0("chr", CHR_B_NUMBER))

#load and arrange Sminor chromosome size data
Smin.chrom.sizes <- read.table("/Users/user/postgrad/lab_docx/fungal_asm_data/genome/Smin.chrsizes.txt",
                               quote="\"", comment.char="",
                               col.names = c("chrom", "size"))

#View(Smin.chrom.sizes)
# ── Bin sizes ──────────────────────────────────────────────────────────────────
binsize     <- 10000   # resolution of the raw data (10 kb) — do not change
plotbinsize <- 100000  # display bin size — increase to merge more tiles (e.g. 50000, 100000, 200000)
# ──────────────────────────────────────────────────────────────────────────────

#build genome matrix using the *plot* bin size for offsets/boundaries
Smin.chrom.sizes <- Smin.chrom.sizes %>%
  mutate(
    n_bins = ceiling(size / plotbinsize),
    offset = cumsum(lag(n_bins, default = 0))
  )
Smin_total_bins <- sum(Smin.chrom.sizes$n_bins)
Smin_offsets    <- setNames(Smin.chrom.sizes$offset, Smin.chrom.sizes$chrom)
Smin_n_bins     <- setNames(Smin.chrom.sizes$n_bins,  Smin.chrom.sizes$chrom)

# parse_bin now maps a raw region string → plot-resolution bin index
parse_bin <- function(region_str) {
  start <- as.numeric(sub("-.*", "", region_str))
  floor(start / plotbinsize)          # <-- uses plotbinsize, not binsize
}

#convert df to genome-wide plot-bin coordinates and aggregate z-scores
print(names(Smin_offsets))
print(head(Smin_zscore$CHR_A_NUMBER))

Smin_plot_df <- Smin_zscore %>%
  mutate(
    bin_a = parse_bin(CHR_A_REG) + Smin_offsets[CHR_A_NUMBER],
    bin_b = parse_bin(CHR_B_REG) + Smin_offsets[CHR_B_NUMBER]
  ) %>%
  group_by(bin_a, bin_b) %>%                        # aggregate within each coarser tile
  summarise(Sminor = mean(Sminor, na.rm = TRUE), .groups = "drop")
View(Smin_plot_df)
#lay out plot
Smin_chrom_breaks <- Smin.chrom.sizes %>%
  mutate(
    mid   = Smin_offsets + n_bins / 2,
    start = Smin_offsets,
    end   = Smin_offsets + n_bins
  )


Smin_plot <- ggplot(Smin_plot_df, aes(x = bin_a, y = bin_b, fill = Sminor)) +
  geom_tile() +
  #chr boundaries
  geom_vline(
    xintercept = Smin_chrom_breaks$start[-1],
    colour = "grey40", linewidth = 0.3, linetype = "solid"
  ) +
  geom_hline(
    yintercept = Smin_chrom_breaks$start[-1],
    colour = "grey40", linewidth = 0.3, linetype = "solid"
  ) +
  scale_x_continuous(
    breaks = Smin_chrom_breaks$mid,
    labels = Smin_chrom_breaks$chrom,
    expand = c(0, 0),
  ) +
  scale_y_continuous(
    breaks = Smin_chrom_breaks$mid,
    labels = Smin_chrom_breaks$chrom,
    expand = c(0, 0),
  ) +
  scale_fill_viridis(
    name     = "z-score",
    option   = "G",
    na.value = "white",
    direction = -1,
    limits = c(-5, 5)
  ) +
  coord_fixed() +
  labs(
    title = "S. minor NHCC",
    subtitle = paste0("Display resolution: ", plotbinsize / 1000, " kb bins"),
    x = "Chromosome (A)",
    y = "Chromosome (B)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 8),
    panel.border    = element_rect(colour = "grey60", fill = NA, linewidth = 0.5),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5, size = 9, colour = "grey40"),
    legend.position = "right"
  )
Smin_plot

#ggsave("/Users/user/postgrad/lab_docx/fungal_asm_data/genome/hic_020426/interchromosomal/Scau_NHCC.pdf",
#       Scau_plot,
#       width = 9, height = 7, units = "in")

