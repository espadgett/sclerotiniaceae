library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(KEGGREST)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(enrichplot)
library(ggupset)

#process BLASTkoala (web server) output
#ko terms ONLY file (2 column with protein name/ko term)
  #currently, KO is oriented to hsa (default of BLASTkoala)

#reference used here for SM is SS (ssl) SC is BC (bfu)
sm_ko_raw <- read_delim("postgrad/lab_docx/fungal_asm_data/secretome/final/koala/sample_koala/sample_ko_definition.txt", 
                                                delim = "\t", escape_double = FALSE, 
                                                col_names = FALSE, trim_ws = TRUE)
sm_ko_raw <- sm_ko_raw[,c(1,2)]
colnames(sm_ko_raw) <- c("gene_id", "KO")
sm_ko_raw <- sm_ko_raw %>%
  filter(!is.na(sm_ko_raw$KO)) #remove genes without KO term
head(sm_ko_raw)

#pull bfu-specific KO terms
ko_hsa_ids <- unique(sm_ko_raw$KO)
ko_to_ssl <- keggLink("ssl", paste0("ko:", ko_hsa_ids))
ko_ssl <- tibble(
  KO  = gsub("ko:", "", names(ko_to_ssl)),
  ssl = gsub("ssl:", "", ko_to_ssl)
)
View(ko_ssl)

annotated_df <- sm_ko_raw %>%
  left_join(ko_ssl, by = "KO")
View(annotated_df)

ssl_to_pathway <- keggLink("pathway", paste0("ssl:", ko_ssl$ssl))
pathway_df <- tibble(
  bfu     = gsub("ssl:", "", names(ssl_to_pathway)),
  pathway = gsub("path:", "", ssl_to_pathway))

#ANALYSIS: enrichment
enrich_res <- enrichKEGG(
  gene         = annotated_df$ssl, #gene/orgnaism terms must match (e.g. annotated_df$KO -> "ko")
  organism     = "ssl",       # use "ko" for KO-level, or "bfu" for organism-level
  keyType      = "kegg",
  pvalueCutoff = 0.1
)

#barplot(enrich_res, showCategory = 20) +
#  ggtitle("S. minor KEGG Enrichment Among Effector Candidates") +
#  xlab("Protein Count")

enrich_df <- as.data.frame(enrich_res) %>%
  mutate(
    GeneRatio_numeric = sapply(GeneRatio, function(x) eval(parse(text = x))),
    Description = fct_reorder(Description, Count)  # order bars by count
  ) %>%
  slice_min(p.adjust, n = 30)   # top 30 by adjusted p-value
View(enrich_df)

enrichmentbar <- ggplot(enrich_df, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.7) +   # <-- bar width
    scale_fill_gradient(
    high = "#085041",
    low = "#bbd491",
    name = "Adjusted\np-value"
  ) +
  labs(
    title    = "KEGG Enrichment Among Effector Candidates",
    subtitle = "Top 30 terms by adjusted p-value; KEGG Reference: Sclerotinia sclerotiorum (ssl)",
    x        = "Protein Count",
    y        = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.y      = element_text(size = 10),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 8)
  )
  
print(enrichmentbar)

#saving settings
ggsave("/Users/user/postgrad/lab_docx/fungal_asm_data/secretome/final/koala/sminor_koala/keggenrichmentbar_sminor.pdf",
           enrichmentbar,
           width = 8, height = 4, unit = "in", limitsize = F)

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(KEGGREST)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(enrichplot)
library(ggupset)

#process BLASTkoala output
#ko terms ONLY file (2 column with protein name/ko term)
#currently, KO is oriented to hsa (default of BLASTkoala)

#reference used here for SM is SS (ssl) SC is BC (bfu)
sc_ko_raw <- read_delim("~/postgrad/lab_docx/fungal_asm_data/secretome/final/koala/scaulophylli_koala/scaulophylli_ko_definition.txt", 
                        delim = "\t", escape_double = FALSE, 
                        col_names = FALSE, trim_ws = TRUE)
sc_ko_raw <- sc_ko_raw[,c(1,2)]
colnames(sc_ko_raw) <- c("gene_id", "KO")
sc_ko_raw <- sc_ko_raw %>%
  filter(!is.na(sc_ko_raw$KO)) #remove genes without KO term
head(sc_ko_raw)

#pull bfu-specific KO terms
ko_hsa_ids <- unique(sc_ko_raw$KO)
ko_to_bfu <- keggLink("bfu", paste0("ko:", ko_hsa_ids))
ko_bfu <- tibble(
  KO  = gsub("ko:", "", names(ko_to_bfu)),
  bfu = gsub("bfu:", "", ko_to_bfu)
)
View(ko_bfu)

annotated_df <- sc_ko_raw %>%
  left_join(ko_bfu, by = "KO")
View(annotated_df)

bfu_to_pathway <- keggLink("pathway", paste0("bfu:", ko_bfu$bfu))
pathway_df <- tibble(
  bfu     = gsub("bfu:", "", names(bfu_to_pathway)),
  pathway = gsub("path:", "", bfu_to_pathway))

#ANALYSIS: enrichment
enrich_res <- enrichKEGG(
  gene         = annotated_df$bfu, #gene/orgnaism terms must match (e.g. annotated_df$KO -> "ko")
  organism     = "bfu",       # use "ko" for KO-level, or "bfu" for organism-level
  keyType      = "kegg",
  pvalueCutoff = 0.1
)

#barplot(enrich_res, showCategory = 20) +
#  ggtitle("S. minor KEGG Enrichment Among Effector Candidates") +
#  xlab("Protein Count")

enrich_df <- as.data.frame(enrich_res) %>%
  mutate(
    GeneRatio_numeric = sapply(GeneRatio, function(x) eval(parse(text = x))),
    Description = fct_reorder(Description, Count)  # order bars by count
  ) %>%
  slice_min(p.adjust, n = 30)   # top 30 by adjusted p-value
View(enrich_df)

enrichmentbar <- ggplot(enrich_df, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.7) +   # <-- bar width
  scale_fill_gradient(
    high = "#085041",
    low = "#bbd491",
    name = "Adjusted\np-value"
  ) +
  labs(
    title    = "KEGG Enrichment Among Effector Candidates",
    subtitle = "Top 30 terms by adjusted p-value; KEGG Reference: Botrytis cinerea (bfu)",
    x        = "Protein Count",
    y        = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.y      = element_text(size = 10),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 8)
  )

print(enrichmentbar)

#saving settings
ggsave("/Users/user/postgrad/lab_docx/fungal_asm_data/secretome/final/koala/scaulophylli_koala/keggenrichmentbar_scaulophylli.pdf",
       enrichmentbar,
       width = 8, height = 4, unit = "in", limitsize = F)

