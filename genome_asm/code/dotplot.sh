#!/bin/bash

#dotplot.sh to visualize genome asm
module load minimap2/2.28

minimap2 -x asm5 -t 36 /home/user/scratch/streptotinia_caulophylli/hic_asm/SC_hifihic_020326.asm.hic.p_ctg.fa /home/user/scratch/streptotinia_caulophylli/ref_genomes/Bcinerea/GCA_000143535.4_ASM14353v4_genomic.fna > SC_BC.minimap2.paf 

#pafplot installation
export PATH="/home/user/.cargo/bin:$PATH"
pafplot --version
#pafplot 0.1.0

pafplot SC_BC.minimap2.paf #written to SC_BC.minimap2.paf.png