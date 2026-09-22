#!/bin/bash
# ---------------------------------------------------------------------
# PREPROCESSING 0.0
# ---------------------------------------------------------------------
# S. caulophylli NHCC run with chromosome-only 200526.fasta file (not shown fully, as code is the same from [PREPROCESSING - BWAMEM] and onwards
cd /home/user/scratch/streptotinia_caulophylli/hic_asm
awk '/^>/ { keep = /^>chr/ } keep' Scaulophylli_200526.fasta > Scaulophylli_200526.chronly.fasta

# S. minor NHCC run with chromosome-only fasta as well
awk '
/^>/{
    if ($0 ~ /^>scaffold_([1-9]|1[0-6])$/) {
        keep = 1
        n = substr($0, 11)          # strip ">scaffold_"
        print ">chr" n
        next
    } else {
        keep = 0
        next
    }
}
keep
' Sminor_sorted.fasta > Sminor_chronly.fasta
# ---------------------------------------------------------------------
# PREPROCESSING 0.1 - BWAMEM
# ---------------------------------------------------------------------
# MODULE LOAD
module load StdEnv/2020 bwa/0.7.17 samtools/1.15.1 #bwa version used by signature pipeline

# FILES/PATHS
ref_genome="/home/user/scratch/streptotinia_caulophylli/ref_genomes/Sminor/Sminor_chronly.fasta"
hic_read1="/home/user/projects/investigator/user/fungal_genomes/raw_data/HiC/Smin_22LYWWLT4_8_1_TGGTCCTT-GATGCTAC_150bp.concat_chastity_passed.fastq.gz" 
hic_read2="/home/user/projects/investigator/user/fungal_genomes/raw_data/HiC/Smin_22LYWWLT4_8_2_TGGTCCTT-GATGCTAC_150bp.concat_chastity_passed.fastq.gz"
sp_prefix="Smin"
# ---------------------------------------------------------------------
cd /home/user/scratch/streptotinia_caulophylli/hic_asm/interchr_interact/sminor_rerun

mkdir signature_analysis
cd /home/user/scratch/streptotinia_caulophylli/hic_asm/interchr_interact/sminor_rerun
###bwa-mem maps raw HiC reads to reference genome
bwa index $ref_genome
bwa mem $ref_genome $hic_read1 $hic_read2 > "$sp_prefix"_hic_aln.sam

###mapping statistics pulled by samtools view/flagstat
samtools view -bhS -o "$sp_prefix"_hic_aln.bam "$sp_prefix"_hic_aln.sam #converting sam to bam
samtools flagstat "$sp_prefix"_hic_aln.bam #summarize alignment flags

# ---------------------------------------------------------------------
# PREPROCESSING 0.2 - PAIRTOOLS
# ---------------------------------------------------------------------
# MODULE LOAD
module load python/3.11 samtools htslib/1.18 

# FILES/PATHS
sp_prefix="Smin"
# ---------------------------------------------------------------------
#make chrom.sizes file
cd /home/user/scratch/streptotinia_caulophylli/hic_asm/interchr_interact/sminor_rerun
samtools view -H "$sp_prefix"_hic_aln.bam | grep "@SQ" | wc #checking to see if there is "chrom" information to extract
samtools view -H "$sp_prefix"_hic_aln.bam | grep "@SQ" | \
    awk '{split($2,a,":"); split($3,b,":"); print a[2] "\t" b[2]}' | \
    sort -t$'\t' -k1,1V > chrom.sizes

cd /home/user/scratch/streptotinia_caulophylli/hic_asm/interchr_interact/sminor_rerun
source ~/envs/pairtools_env/bin/activate

echo "confirming pairtools load at `date`"
which pairtools
pairtools --version
echo ""

echo "running pairtools parse at `date`" #organizing into HiC pairs
samtools view -h "$sp_prefix"_hic_aln.bam | pairtools parse -c chrom.sizes -o parsed.pairs.gz --drop-sam

echo "sorting pairs at `date`" #sort pairs by location
pairtools sort -o parsed.sorted.pairs.gz parsed.pairs.gz

echo "deduplicating pairs at `date`"
pairtools dedup -o parsed.sorted.dedup.pairs.gz --mark-dups parsed.sorted.pairs.gz

echo "selecting UU/UR/RU pair types at `date`" #saving unique pairs (R = rescued, U = unique)
pairtools select \
    '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
    -o parsed.sorted.dedup.selected.pairs.gz \
    parsed.sorted.dedup.pairs.gz

echo "splitting pairs at `date`" #make .pairs output
pairtools split --output-pairs split.pairs parsed.sorted.dedup.selected.pairs.gz

echo "running final selection at `date`"
pairtools select 'True' -o parsed.sorted.dedup.selected.final.pairs.gz \
    parsed.sorted.dedup.selected.pairs.gz

echo ""
echo "pairtools pipe completed at: $(date)"

deactivate

# ---------------------------------------------------------------------
# PREPROCESSING 0.3 - PAIRIX + COOLER
# ---------------------------------------------------------------------
# MODULE LOAD
module load StdEnv/2020 pairix/0.3.7
echo "confirming pairix load at `date`"
pairix --version #0.3.7
echo ""

# FILES/PATHS
sp_prefix="Smin"
# ---------------------------------------------------------------------
cd /home/user/scratch/streptotinia_caulophylli/hic_asm/interchr_interact/sminor_rerun

#pairix indexing
file parsed.sorted.dedup.selected.final.pairs.gz #confirm bgzip

echo "running pairix indexing at `date`" #bin a pairix-indexed contact list file
pairix parsed.sorted.dedup.selected.final.pairs.gz
ls -l parsed.sorted.dedup.selected.final.pairs.gz.px2

#cooler setup
ls -l chrom.sizes #confirm chrom.sizes file - other requisite of cooler

module load python/3.11
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index cooler pypairix

#cooler
echo "running cload pairix at `date`" #bin a pairix-indexed contact list file
cooler cload pairix chrom.sizes:10000 parsed.sorted.dedup.selected.final.pairs.gz cool.output

echo "running balance at `date`" #out-of-core matrix balancing
cooler balance cool.output

echo "running dump at `date`"	#dump cooler's data into a text stream
#cooler dump --balanced -o cool.output300426.txt cool.output
#NOTE that --balanced and --output options were not included on Signature Github but (1) file had incorrect column number for conversion to Signature format wihtout --balanced and (2) I want an output file, I don't want to parse through my .out note file

cooler dump --join --balanced -o cool.output.joined.txt cool.output

# ---------------------------------------------------------------------
# PREPROCESSING 0.4 - COOLER DUMP TO SIGNATURE
# ---------------------------------------------------------------------
# MODULE LOAD
module load python/3.11 scipy-stack StdEnv/2023

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install biopython
pip install Bio
# ---------------------------------------------------------------------
cd /home/user/scratch/streptotinia_caulophylli/hic_asm/interchr_interact/sminor_rerun

echo "running custom coolerdumptosignature.py at `date`"
python ~/projects/investigator/user/bin/Signature/Pipeline/Pre-Processing_HiC/coolerdumpToSignature.py \
/home/user/scratch/streptotinia_caulophylli/hic_asm/interchr_interact/sminor_rerun/cool.output.joined.txt \
/home/user/scratch/streptotinia_caulophylli/hic_asm/interchr_interact/sminor_rerun/ \
10000

#after script has completed (I used a SLURM system), in terminal, do:

#cd /home/user/scratch/streptotinia_caulophylli/hic_asm/interchr_interact/sminor_rerun
#mkdir -p Smin
#mv trans.10000_iced.sorted.txt Smin/