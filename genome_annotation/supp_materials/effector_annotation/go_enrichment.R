library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(readr)
library(enrichplot)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationHub")
library(AnnotationHub)
library(AnnotationDbi)
library(GOSemSim)
library(patchwork)


pos_sample_go <- read_delim("postgrad/lab_docx/fungal_asm_data/secretome/phib-blast/go/phi_pos_go.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         col_names = FALSE, trim_ws = TRUE)
#View(pos_sample_go)
colnames(pos_sample_go) <- c("sample_id", "go_id")

term_list <- go2term(pos_sample_go$go_id)
ont_list <- go2ont(pos_sample_go$go_id)
go_df <- left_join(pos_sample_go, term_list, by = "go_id")
go_df <- left_join(go_df, ont_list, by = "go_id")
head(go_df)

#enrichment
#creating OrgDb object
hub <- AnnotationHub() #create AnnotationHub directory at /Users/user/Library/Caches/org.R-project.R/R/AnnotationHub
q <- query(hub, "Sclerotinia") #selected reference
id <- q$ah_id[length(q)] #id = AH119190
hub[["AH119190"]] #indeed S. sclerotiorum
Sclerotinia_orgdb <- hub[[id]]
class(Sclerotinia_orgdb) #OrgDB

#mapping GO ids to ENTREZ ids
entrez_map <- select(
  Sclerotinia_orgdb,
  keys = as.character(go_df$go_id),
  keytype = "GO",
  columns = c("ENTREZID", "SYMBOL")
)
View(entrez_map)

mapped_go_df <- merge(go_df, entrez_map, by.x="go_id", by.y="GO", all.x = TRUE)
View(mapped_go_df)

#assessing semantic similarity of GO terms - MF ontology
SS_MFref <- godata(Sclerotinia_orgdb, ont = "MF")
SS_BPref <- godata(Sclerotinia_orgdb, ont = "BP")
SS_CCref <- godata(Sclerotinia_orgdb, ont = "CC")
class(SS_MFref) #GOSemSim

go_mf_similarity <- mgoSim(go_df$go_id, go_df$go_id, semData = SS_MFref, measure="Wang", combine=NULL)
go_mf_similarity %>%
  as_tibble() %>%
  rowid_to_column((var = "X")) %>%
  gather(key = "Y", value = "Z", -1) %>%
  mutate(Y=as.numeric(gsub("V","",Y))) %>%

  ggplot(aes(X, Y, fill= Z)) + 
    geom_tile() +
    theme_minimal()+
    scale_fill_gradient(low="white", high="#ff3d03") +
    labs(x = NULL, y = NULL, title = "Molecular Function GO Similarity", fill = "Term Similarity")

go_bp_similarity <- mgoSim(go_df$go_id, go_df$go_id, semData = SS_BPref, measure="Wang", combine=NULL)
go_bp_similarity %>%
  as_tibble() %>%
  rowid_to_column((var = "X")) %>%
  gather(key = "Y", value = "Z", -1) %>%
  mutate(Y=as.numeric(gsub("V","",Y))) %>%
  
  ggplot(aes(X, Y, fill= Z)) + 
  geom_tile() +
  theme_minimal()+
  scale_fill_gradient(low="white", high="#ff3d03") +
  labs(x = NULL, y = NULL, title = "Biological Process GO Similarity", fill = "Term Similarity")

go_cc_similarity <- mgoSim(go_df$go_id, go_df$go_id, semData = SS_CCref, measure="Wang", combine=NULL)
go_cc_similarity %>%
  as_tibble() %>%
  rowid_to_column((var = "X")) %>%
  gather(key = "Y", value = "Z", -1) %>%
  mutate(Y=as.numeric(gsub("V","",Y))) %>%
  
  ggplot(aes(X, Y, fill= Z)) + 
  geom_tile() +
  theme_minimal()+
  scale_fill_gradient(low="white", high="#ff3d03") +
  labs(x = NULL, y = NULL, title = "GO Similarity", fill = "Term Similarity")

#GO/gene enrichment
#need to pull genes from ENTREZ terms
gene_list <- unique(mapped_go_df$ENTREZID)

term2gene <- go_df[,c("go_id", "sample_id")]
term2name <- go_df[,c("go_id", "Term")

enrich_go_res <- enrichGO(
  gene          = gene_list,
  OrgDb         = Sclerotinia_orgdb,
  keyType       = "ENTREZID",   # or whatever matches your OrgDb
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
head(enrich_go_res)

#enrichment analyses
  #low-redundancy enrichment network
all_gene_list <- mapped_go_df$ENTREZID
lr_res <- enrichGO(all_gene_list, ont="BP", OrgDb = Sclerotinia_orgdb)
View(lr_res)

lr_d <- pairwise_termsim(lr_res)
lr_ds <- simplify(lr_res, cutoff=0.7, by = "p.adjust", select_fun=min)
lr_ds <- pairwise_termsim(lr_ds) #recreate similarity matrix

lr_m1 <- emapplot(lr_d, layout = "fr", node_label = "category")
lr_m1
lr_m2 <- emapplot(lr_ds, layout = "fr", node_label = "category")
lr_m2

ggsave("/Users/user/postgrad/lab_docx/fungal_asm_data/secretome/figures/goenrichmentnetwork_sminor.pdf",
       lr_m1,
       width = 8, height = 9, unit = "in", limitsize = F)
ggsave("/Users/user/postgrad/lab_docx/fungal_asm_data/secretome/figures/goenrichmentnetwork_sminor_reducedsimilarity.pdf",
       lr_m2,
       width = 8, height = 9, unit = "in", limitsize = F)

#plotting bar plot top GO terms (from reduced redundancy df)
plot_df <- as.data.frame(lr_res)
#View(plot_df)
plot_df <- plot_df %>%
  arrange(p.adjust) %>%
  slice_head(n = 30) %>%
  mutate(
    Description = factor(Description, levels = rev(Description)),
    log_padj    = -log10(p.adjust),
    GeneRatio   = sapply(GeneRatio, function(x) eval(parse(text = x)))
  )

names(plot_df)[1] <- "go_id"
plot_df <- plot_df %>%
  left_join(mapped_go_df %>% dplyr::select(go_id, Ontology), by = "go_id")
plot_df <- plot_df %>%
  arrange(GeneRatio)

enrichment_dot <- ggplot(plot_df, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = log_padj)) +
  geom_point() +
  scale_color_gradient(high = "#085041",
                       low = "#bbd491") +
  labs(
    x     = "GeneRatio",
    y     = NULL,
    title = "GO Enrichment Among Effector Candidates",
    subtitle = "Top 30 terms by adjusted p-value; GO reference: Sclerotinia sclerotiorum (AH119190)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y      = element_text(size = 10),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 8)
  )
enrichment_dot

ggsave("/Users/user/postgrad/lab_docx/fungal_asm_data/secretome/figures/goenrichmentdot.pdf",
       enrichment_dot,
       width = 8, height = 10, unit = "in", limitsize = F)