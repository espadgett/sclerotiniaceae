#R

library(devtools)
#install.packages("data.table")
library(data.table)
library(BiocManager)
library(GENESPACE)

system("python --version") #Python 3.11.5
system("diamond --version") #diamond version 2.1.11
system("which orthofinder") #~/scratch/streptotinia_caulophylli/genespace/OrthoFinder_source/orthofinder #OrthoFinder v2.55

###############################################
# -- change paths to those valid on your system
genomeRepo <- "/home/user/scratch/streptotinia_caulophylli/genespace/run280726/rawGenomes"
wd <- "/home/user/scratch/streptotinia_caulophylli/genespace/run280726/output"
path2mcscanx <- "/home/user/scratch/streptotinia_caulophylli/genespace/MCScanX/"
###############################################

# -- download raw data from NCBI for human and chicken genomes
#dir.create(genomeRepo)
#rawFiles <- download_exampleData(filepath = genomeRepo)

# -- parse the annotations to fastas with headers that match a gene bed file
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = c("sscl","bcin"),
  genomeIDs = c("ssclerotiorum","bcinerea"),
  gffString = "gff3",
  faString = "fa",
  headerSep = " ", 
  genespaceWd = wd)
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = c("smin","scau"),
  genomeIDs = c("sminor","scaulophylli"),
  gffString = "gff3",
  faString = "fa",
  headerEntryIndex = 1, 
  gffIdColumn = "ID",
  genespaceWd = wd)
   
# -- initalize the run and QC the inputs
gpar <- init_genespace(
  wd = wd, 
  nCores = 30,
  path2mcscanx = path2mcscanx,
  dotplots = "never"
  )

# -- accomplish the run
out <- run_genespace(gpar, overwrite = T)
save.image("genespace_run280726.RData")