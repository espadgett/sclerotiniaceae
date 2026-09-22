#**NOTE** the genespace parameter object is returned or can be loaded
#        into R via
#        `load('/home/user/scratch/streptotinia_caulophylli/genespace/run280726/output/results/gsParams.rda',
#        verbose = TRUE)`. Then you can customize your riparian plots by
#        calling `plot_riparian(gsParam = gsParam, ...)`. The source
#        data and ggplot2 objects are also stored in the /riparian
#        directory and can also be accessed by `load(...)`.
#**NOTE** To query genespace results by position or gene, use
#        `query_genespace(...)`. See specifications in ?query_genespace
#        for details.

# replicate environment from genespace_initalize.sh
module load StdEnv/2023 gcc/14.3 python diamond/2.1.11 r/4.5.0
source ~/envs/genespace_env/bin/activate

export PATH=$PATH:/home/user/scratch/streptotinia_caulophylli/genespace/MCScanX
export PATH=$PATH:home/user/scratch/streptotinia_caulophylli/genespace/MCScanX/downstream_analyses
export PATH=$PATH:/home/user/scratch/streptotinia_caulophylli/genespace/OrthoFinder_source
export R_LIBS_USER=/home/user/envs/R/x86_64-pc-linux-gnu-library/4.5

# in interactive session, all modules/paths loaded already
library(devtools)
#install.packages("data.table")
library(data.table)
library(BiocManager)
library(GENESPACE)

load("genespace_run280726.RData")

pdf(file = 'riparian_020826_10.pdf')

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("#742C14FF","#274637FF", "#58A449FF", "#CEC917FF", "#E9D097FF" ))

invchr <- data.frame(
  genome = c("sminor", "sminor", "sminor", "sminor", "sminor", "sminor",
  "bcinerea","bcinerea","bcinerea","bcinerea","bcinerea","bcinerea",
  "scaulophylli","scaulophylli","scaulophylli","scaulophylli","scaulophylli"
  ),
  chr = c("Sc_min_sc_1", "Sc_min_sc_4", "Sc_min_sc_5","Sc_min_sc_6", "Sc_min_sc_11", "Sc_min_sc_13",
  "NC_037316.1","NC_037310.1","NC_037314.1","NC_037320.1","NC_037319.1","NC_037313.1",
  "St_cau_sc_2","St_cau_sc_9","St_cau_sc_7","St_cau_sc_12","St_cau_sc_13"))
  
ripDat <- plot_riparian(
  gsParam = out,
  refGenome = "ssclerotiorum", 
  genomeIDs = c("scaulophylli","bcinerea","sminor","ssclerotiorum"),
  forceRecalcBlocks = FALSE,
  palette = customPal,
  chrFill = "lightgrey",
  minChrLen2plot = 0,
	chrLabFun = function(x) {
	  # strict one-to-one mapping
	  map <- c(
    "CP017814.1" = "1",  "CP009805.1" = "1",
    "CP017815.1" = "2",  "CP009806.1" = "2",
    "CP017816.1" = "3",  "CP009807.1" = "3",
    "CP017817.1" = "4",  "CP009808.1" = "4",
    "CP017818.1" = "5",  "CP009809.1" = "5",
    "CP017819.1" = "6",  "CP009810.1" = "6",
    "CP017820.1" = "7",  "CP009811.1" = "7",
    "CP017821.1" = "8",  "CP009812.1" = "8",
    "CP017822.1" = "9",  "CP009813.1" = "9",
    "CP017823.1" = "10", "CP009814.1" = "10",
    "CP017824.1" = "11", "CP009815.1" = "11",
    "CP017825.1" = "12", "CP009816.1" = "12",
    "CP017826.1" = "13", "CP009817.1" = "13",
    "CP017827.1" = "14", "CP009818.1" = "14",
    "CP017828.1" = "15", "CP009819.1" = "15",
    "CP017829.1" = "16", "CP009820.1" = "16",
    "CP009821.1" = "17",
    "CP009822.1" = "18",

    "NC_037316.1" = "7",
    "NC_037310.1" = "1",
    "NC_037312.1" = "3",
    "NC_037315.1" = "6",
    "NC_037314.1" = "5",
    "NC_037319.1" = "10",
    "NC_037311.1" = "2",
    "NC_037317.1" = "8",
    "NC_037318.1" = "9",
    "NC_037313.1" = "4",
    "NC_037321.1" = "12",
    "NC_037324.1" = "15",
    "NC_037320.1" = "11",
    "NC_037322.1" = "13",
    "NC_037325.1" = "16",
    "NC_037323.1" = "14",
    
    #contig swapping
	  "St_cau_sc_3" = "contig5",
	  "St_cau_sc_4" = "contig3",
	  "St_cau_sc_5" = "contig7",
	  "St_cau_sc_7" = "contig9",
	  "St_cau_sc_8" = "contig10",
	  "St_cau_sc_9" = "contig8",
	  "St_cau_sc_10" = "contig12",
	  "St_cau_sc_12" = "contig4",
	  "St_cau_sc_13" = "contig15",
	  "St_cau_sc_14" = "contig13",
	  "St_cau_sc_15" = "contig16",
	  "St_cau_sc_16" = "contig14"
	  )

    # apply mapping
	  out <- map[x]
	  out[is.na(out)] <- x[is.na(out)]
	
	  # universal cleanup
	  out <- gsub("St_cau_sc_|Sc_min_sc_|_|contig", "", out)
	  return(out)
	},
  invertTheseChrs = invchr,
	addThemes = ggthemes,
  #useOrder = FALSE,
  #useRegions = FALSE,
  backgroundColor = NULL,
  scalePlotHeight = 1,
  scalePlotWidth = 1)
dev.off()