library(rtracklayer)
library(ggplot2)
library(dplyr)

# ── CONFIG ───────────────────────────────────────────────────
GFF3_FILE   <- "/Users/user/postgrad/lab_docx/fungal_asm_data/genome/hic_020426/Scaulophylli_Scaulophylli.gff3"  # ← your GFF3 path
FEATURE     <- "gene"                  # gene / mRNA / CDS …
BIN_SIZE    <- 0.25e5                     # bin width in bp (1 Mb)
FILL_LOW    <- "#ece6ac"              # low-density colour
FILL_HIGH   <- "#2C715FFF"              # high-density colour
CHROM_ORDER <- c("St_cau_sc_1","St_cau_sc_2","St_cau_sc_3","St_cau_sc_4",
                 "St_cau_sc_5","St_cau_sc_6","St_cau_sc_7","St_cau_sc_8",
                 "St_cau_sc_9","St_cau_sc_10","St_cau_sc_11","St_cau_sc_12",
                 "St_cau_sc_13","St_cau_sc_14","St_cau_sc_15","St_cau_sc_16")
                                      # NULL = auto sort by length (largest first)
                                      # or e.g.: c("chr1","chr2", ... "chr16")
# ─────────────────────────────────────────────────────────────


# ── 1. Import GFF3 ───────────────────────────────────────────
message("Reading GFF3 …")
gff   <- import(GFF3_FILE)
genes <- gff[gff$type == FEATURE]

if (length(genes) == 0)
  stop("No '", FEATURE, "' features found — check your FEATURE variable.")

genes_df <- as.data.frame(genes) |>
  select(seqnames, start, end) |>
  rename(chrom = seqnames)

genes_df <- genes_df |>
  filter(chrom %in% CHROM_ORDER)


# ── 2. Chromosome lengths (from max GFF3 coordinate) ─────────
# Uses observed max — avoids the karyoploteR length-mismatch bug
chrom_lengths <- genes_df |>
  group_by(chrom) |>
  summarise(chrom_len = max(end), .groups = "drop")

# ── 3. Build fixed-width bins per chromosome ─────────────────
bins <- chrom_lengths |>
  rowwise() |>
  reframe(
    chrom     = chrom,
    bin_start = seq(1, chrom_len, by = BIN_SIZE),
    bin_end   = pmin(seq(1, chrom_len, by = BIN_SIZE) + BIN_SIZE - 1,
                     chrom_len)
  )


# ── 4. Count overlapping features per bin ────────────────────
message("Counting features per bin …")

bins_counted <- bins |>
  rowwise() |>
  mutate(
    gene_count = sum(
      genes_df$chrom  == chrom &
        genes_df$start  <= bin_end &
        genes_df$end    >= bin_start
    )
  ) |>
  ungroup() |>
  mutate(
    chrom   = as.character(chrom),          # ← force character before rename
    density = gene_count / (BIN_SIZE / 1e6),
    bin_mid = (bin_start + bin_end) / 2
  )

# ── 5. Forcibly chromosomes from Funannotate reorder (returning to B.cinerea-aligned naming convention) ───────────────────────────────────────────────────
chrom_rename <- c("St_cau_sc_1" = "1",
                  "St_cau_sc_2" = "2",
                  "St_cau_sc_3" = "5",
                  "St_cau_sc_4" = "3",
                  "St_cau_sc_5" = "7",
                  "St_cau_sc_6" = "6",
                  "St_cau_sc_7" = "9",
                  "St_cau_sc_8" = "10",
                  "St_cau_sc_9" = "8",
                  "St_cau_sc_10" = "12",
                  "St_cau_sc_11" = "11",
                  "St_cau_sc_12" = "4",
                  "St_cau_sc_13" = "15",
                  "St_cau_sc_14" = "13",
                  "St_cau_sc_15" = "16",
                  "St_cau_sc_16" = "14")

bins_counted <- bins_counted |>
  mutate(chrom = chrom_rename[as.character(chrom)])  

# ── 6. Set numeric chromosome order ──────────────────────────
new_chrom_order <- bins_counted |>
  distinct(chrom) |>
  mutate(num = as.integer(chrom)) |>
  arrange(num) |>
  pull(chrom)

bins_counted <- bins_counted |>
  mutate(chrom = factor(chrom, levels = rev(new_chrom_order)))

# ── 7. Plot ───────────────────────────────────────────────────
p <- ggplot(bins_counted,
            aes(x = bin_mid, y = chrom, fill = density)) +  # ← fill = density
  
  geom_tile(aes(width = BIN_SIZE, height = 0.75), colour = NA) +
  
  scale_fill_gradient(
    low  = FILL_LOW,
    high = FILL_HIGH,
    name = "Gene density\n(Genes / Mb)"  # ← legend title
  ) +
  
  scale_x_continuous(
    labels = \(x) paste0(round(x / 1e6), " Mb"),
    expand = c(0.01, 0)
  ) +
  
  labs(
    title    = "Gene Density — Streptotinia caulophylli",   # ← your title
    #subtitle = paste0("Bin size: ", round(BIN_SIZE / 1e6), " Mb"),
    x        = "Genomic position (Mb)",          # ← x-axis label
    
    y        = "Chromosome"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(face = "bold"),
    legend.position    = "right",
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(colour = "grey50")
  )

print(p)


# ── 7. Save ───────────────────────────────────────────────────
height_in <- max(4, 0.45 * length(chrom_order) + 2)
ggsave("/Users/user/postgrad/lab_docx/fungal_asm_data/genome/hic_020426/gene_density_map.pdf", p,
       width = 12, height = height_in, units = "in", limitsize = FALSE)
message("Saved → gene_density_map.pdf")