#!/bin/bash

cd /home/user/scratch/streptotinia_caulophylli/genespace/OrthoFinder_source
set -euo pipefail 

module load StdEnv/2023 gcc/14.3 python diamond/2.1.11 r/4.5.0
source ~/envs/genespace_env/bin/activate

export PATH=$PATH:/home/user/scratch/streptotinia_caulophylli/genespace/MCScanX
export PATH=$PATH:home/user/scratch/streptotinia_caulophylli/genespace/MCScanX/downstream_analyses
export PATH=$PATH:/home/user/scratch/streptotinia_caulophylli/genespace/OrthoFinder_source
#echo "Confirm Orthofinder:"
#orthofinder -h
#-----------------------------------------------
# ENVIRONMENT STAGING
cd /home/user/scratch/streptotinia_caulophylli/genespace/
mkdir run280726
cd run280726

#copy over most up to date files (as of Jul 28) from /home/user/scratch/streptotinia_caulophylli/hic_asm/for_weijie/
mkdir scau
cp /home/user/scratch/streptotinia_caulophylli/hic_asm/for_weijie/annotations/scaulophylli/Scaulophylli.proteins.fa scau/
cp /home/user/scratch/streptotinia_caulophylli/hic_asm/for_weijie/annotations/scaulophylli/Scaulophylli.gff3 scau/
mkdir smin
cp /home/user/scratch/streptotinia_caulophylli/hic_asm/for_weijie/annotations/sminor/renamed.Sminor_Sminor.proteins.fa smin/
cp /home/user/scratch/streptotinia_caulophylli/hic_asm/for_weijie/annotations/sminor/renamed.Sminor_Sminor.gff3 smin/
#copy up to date files from 
mkdir bcin
cp /home/user/scratch/streptotinia_caulophylli/ref_genomes/Bcinerea/bcinerea.protein.faa bcin/
cd bcin
mv bcinerea.protein.faa bcinerea.protein.fa
cd ..
cp /home/user/scratch/streptotinia_caulophylli/ref_genomes/Bcinerea/bcinerea_genomic.gff bcin/
mv bcin/bcinerea_genomic.gff bcin/bcinerea.gff3
mkdir sscl
cp /home/user/scratch/streptotinia_caulophylli/ref_genomes/Ssclerotiorum_GCA_001857865.1/ncbi_dataset/data/GCA_001857865.1/ssclerotiorum.protein.fa sscl/
cp /home/user/scratch/streptotinia_caulophylli/ref_genomes/Ssclerotiorum_GCA_001857865.1/GCA_001857865.1_ASM185786v1_genomic.gff sscl/
mv sscl/GCA_001857865.1_ASM185786v1_genomic.gff sscl/ssclerotiorum.gff3


mkdir rawGenomes
mv * rawGenomes
mkdir output
#-----------------------------------------------
Rscript /home/user/scratch/streptotinia_caulophylli/genespace/scripts/genespace_280726.R

sacct -j $SLURM_JOBID --format=JobName,MaxRSS,TotalCPU,Elapsed --units=M
