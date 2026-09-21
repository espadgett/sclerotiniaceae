#!/bin/bash

#genome_asm.sh assembles HiFi + HiC sequences, and outputs a few summary QC stats

# ---------------------------------------------------------------------
# ASSEMBLY
# ---------------------------------------------------------------------
# MODULE LOAD
module load fastqc/0.12.1 samtools/1.22.1 StdEnv/2020 hifiasm/0.19.5 minimap2

# FILES/PATHS
hifi_seq="/home/user/scratch/streptotinia_caulophylli/hifi_seq/Scaulophylli_fasta/Scaulophylli_hifireads.fq.gz"
hic_R1="/home/user/scratch/streptotinia_caulophylli/hic_raw/Scaulophylli_hic/23JWGJLT4_6_1_CAGCGATT-TCAGCCTT_150bp_St.concat.fastq.gz"
hic_R2="/home/user/scratch/streptotinia_caulophylli/hic_raw/Scaulophylli_hic/23JWGJLT4_6_2_CAGCGATT-TCAGCCTT_150bp_St.concat.fastq.gz"
# ---------------------------------------------------------------------

echo "Quality check fastq files with FastQC..."
mkdir -p fastqc_output #hic_asm/fastqc_output
fastqc -o fastqc_output \
	-f fastq \
	$hic_R1 \
	$hic_R2


#make SRF argument .txt with absolute paths to fastq files/
echo "Generating de-novo [hifiasm] assembly with HiFi and HiC reads"
echo -e "\n[hifiasm] with light purging/1 haplotype"
hifiasm -o SC_hifihic_l1_136homcov_190526.asm -t32 \
 -l1 --n-hap 1 \
 --hom-cov 136 \
 --h1 $hic_R1 --h2 $hic_R2 \
 $hifi_seq
 
awk '/^S/{print ">"$2;print $3}' SC_hifihic_l1_136homcov_190526.asm.hic.p_ctg.gfa > SC_hifihic_l1_136homcov_190526.asm.hic.p_ctg.fa
 
#optional - purgedups makes no difference after hifiasm options
echo "running [dedups configuration] with hifiasm-l0 option"
duplicate_scan/purge_dups/scripts/pd_config.py -l l0_purge \
-n purge_dups_l0_config.json \
/home/user/scratch/streptotinia_caulophylli/hic_asm/duplicate_scan/hifiasm_reruns/l0/SC_hifihic_l0_180526.asm.hic.p_ctg.fa \
pbfofn.txt
#conti.json created 
echo "running [dedups purging] with hifiasm-l0 option"
	#run purge_dups with -p bash option for operating in SLURM/CC
duplicate_scan/purge_dups/scripts/run_purge_dups.py \
-p bash \
purge_dups_l0_config.json \
duplicate_scan/purge_dups/bin \
S_caulophylli
echo -e "\n[dedups] with hifiasm-l0 option complete."

echo "running [dedups configuration] with hifiasm-l2 option"
duplicate_scan/purge_dups/scripts/pd_config.py -l l2_purge \
-n purge_dups_l2_config.json \
/home/user/scratch/streptotinia_caulophylli/hic_asm/duplicate_scan/hifiasm_reruns/l2/SC_hifihic_l2_180526.asm.hic.p_ctg.fa \
pbfofn.txt
#conti.json created 
echo "running [dedups purging] with hifiasm-l2 option"
	#run purge_dups with -p bash option for operating in SLURM/CC
duplicate_scan/purge_dups/scripts/run_purge_dups.py \
-p bash \
purge_dups_l2_config.json \
duplicate_scan/purge_dups/bin \
S_caulophylli
echo -e "\n[dedups] with hifiasm-l2 option complete."

# ---------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------
# MODULE LOAD
module load StdEnv/2020 gcc/9.3.0 python/3.10 augustus/3.5.0 hmmer/3.3.2 blast+/2.13.0 metaeuk/6 prodigal/2.6.3 r/4.3.1 bbmap/38.86
#virtual environment will be created for BUSCO and QUAST, so no need to load those modules here

# FILES/PATHS
asm="/home/user/scratch/streptotinia_caulophylli/hic_asm/SC_hifihic_020326.asm.hic.p_ctg.fa"
ref_asm="/home/user/scratch/streptotinia_caulophylli/ref_genomes/Bcinerea/GCA_000143535.4_ASM14353v4_genomic.fna"
#bcinerea is ref genome
mkdir -p ~/scratch/streptotinia_caulophylli/busco/

# ---------------------------------------------------------------------
echo -e "\nBegin BUSCO analysis..."
# Generate your virtual environment in $SLURM_TMPDIR.
virtualenv --no-download ${SLURM_TMPDIR}/busco_env
source ${SLURM_TMPDIR}/busco_env/bin/activate
# Install busco and its dependencies.
pip install --no-index --upgrade pip
pip install --no-index --requirement /home/user/scratch/sunflower_genomes/busco/busco-requirements.txt


busco \
--in $asm \
--out ~/scratch/streptotinia_caulophylli/busco/hic_hifi_asm \
--lineage_dataset ascomycota_odb10 \
--mode genome \
--cpu ${SLURM_CPUS_PER_TASK:-1} \
--download_path ~/scratch/sunflower_genomes/busco/busco_downloads/

deactivate

echo -e "\nBeginning QUAST Analysis"
echo ""
module purge
module load StdEnv/2023 python/3.11 scipy-stack 

echo "confirm quast exists:"
avail_wheels quast

virtualenv --no-download $SLURM_TMPDIR/quast_env
source $SLURM_TMPDIR/quast_env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index quast

python $SLURM_TMPDIR/quast_env/bin/quast.py \
  -o /home/user/scratch/streptotinia_caulophylli/quast/sc_hic \
  "$asm" \
  -r "$ref_asm"
	echo "QUAST complete."

#(harmonize naming to reference before proceeding...)
bash renaming_asm.sh #bash or sbatch depending on system

# ---------------------------------------------------------------------
# REPEAT ANNOTATION
# ---------------------------------------------------------------------
# MODULE LOAD
module load StdEnv/2023 gcc/12.3 rmblast hmmer trf

# FILES/PATHS
RM="/scratch/user/streptotinia_caulophylli/repeatmasker_dir/RepeatMasker/RepeatMasker"
asm="/home/user/scratch/streptotinia_caulophylli/hic_asm/Scaulophylli_270426.fa"
# ---------------------------------------------------------------------
echo "Beginning RepeatMasker analysis for repeat annotation..."
$RM -pa $SLURM_CPUS_PER_TASK -species ascomycota -a -inv "$asm"

# ---------------------------------------------------------------------
# TELOMERE ANNOTATION
# ---------------------------------------------------------------------
# MODULE LOAD
echo "Beginning telomere probe scan across genome asm sequences using probe_repeat_scan.py"
module load python/3.14.2

# FILES/PATHS
sc_genome="/home/user/scratch/streptotinia_caulophylli/hic_asm/Scaulophylli_270426.fa"
forward_motif="CCCTAA"
reverse_motif="TTAGGG"
min_rep="5"
output_tsv="/home/user/scratch/streptotinia_caulophylli/hic_asm/annotat_repeat/probe_scan.tsv"
# ---------------------------------------------------------------------
echo "
Usage:
    e.g. python probe_repeat_scan.py -f genome.fasta -m ATCG -r 5 -o results.tsv

Arguments:
    -f / --fasta      Path to the input genome FASTA file
    -m / --motif      Repeat motif to search for (e.g. AT, CGG, TTAGGG)
    -rc / --rc-motif  Reverse complement motif to search
    -r / --repeats    Minimum number of consecutive repeats to report (e.g. 5)
    -o / --output     Output TSV file path (default: repeat_hits.tsv)
    --case-sensitive  Treat motif search as case-sensitive (default: case-insensitive)

Output TSV columns:
    contig          Contig/sequence name (text after '>' up to first whitespace)
    start           1-based start coordinate of the repeat run
    end             1-based end coordinate (inclusive)
    repeat_count    Number of times the motif is consecutively repeated
    strand          '+' (always, since repeat is palindrome-agnostic)
    sequence        The matched repeat sequence"

python probe_repeat_scan.py \
-f "$sc_genome" \
-m "$forward_motif" \
-rc "$reverse_motif" \
-r "$min_rep" \
-o "$output_tsv"

# ------------------------------------------------------------------------------
# MANUAL CORRECTION - on chr9 telomere
# ------------------------------------------------------------------------------
# in terminal, removed 876 bp region preceeding a 22mer telomeric repeat segment in the top arm of chromosome 9 > Scaulophylli_280426.fasta

# ------------------------------------------------------------------------------
# HIC CURATION PART I: YaHS
# ------------------------------------------------------------------------------
# MODULE LOAD
module load bwa/0.7.18
module load bwa-mem2/2.2.1
module load samtools/1.22.1
#yahs manually installed (see above)
yash --version #v1.2.2
module load java/21.0.1

# FILES/PATHS
ASSEMBLY=/home/user/scratch/streptotinia_caulophylli/hic_asm/Scaulophylli_280426.fasta   # your indexed hifiasm fasta
HIC_R1=/home/user/scratch/streptotinia_caulophylli/hic_raw/Scaulophylli_hic/23JWGJLT4_6_1_CAGCGATT-TCAGCCTT_150bp_St.concat.fastq.gz           # Hi-C forward reads
HIC_R2=/home/user/scratch/streptotinia_caulophylli/hic_raw/Scaulophylli_hic/23JWGJLT4_6_2_CAGCGATT-TCAGCCTT_150bp_St.concat.fastq.gz           # Hi-C reverse reads
SAMPLE=Scaulophylli_140526 #output file naming
OUTDIR=yahs_out

#RE="GATC,GANTC" #for Arima v1/2                                 
#RE="AAGCTT" #for HindIII
RE="GATC" #for DpnII

THREADS=${SLURM_CPUS_PER_TASK:-16}
# ------------------------------------------------------------------------------
# STEP 1: BWA-MEM alignment of Hi-C reads

echo "[$(date)] Starting BWA-MEM indexing and alignment..."

mkdir -p ${OUTDIR}

bwa-mem2 index ${ASSEMBLY}
samtools faidx ${ASSEMBLY}

bwa-mem2 mem \
    -t ${THREADS} \
    ${ASSEMBLY} \
    ${HIC_R1} \
    ${HIC_R2} \
    | samtools view \
        -bS \
        -F 0x904 \
        -@ ${THREADS} \
    | samtools sort \
        -n \
        -@ ${THREADS} \
        -o ${OUTDIR}/${SAMPLE}_hic_to_asm.bam

echo "[$(date)] BWA-MEM alignment complete."

# Step 2: Run YAHS scaffolding
	# YAHS outputs:
	#   ${SAMPLE}_scaffolds_final.fa       — final scaffolded assembly
	#   ${SAMPLE}.bin                      — binary contact matrix
	#   ${SAMPLE}_scaffolds_final.agp      — AGP file for coordinate mapping
	
echo "[$(date)] Starting YAHS scaffolding..."

yahs \
    -e ${RE} \
    -o ${OUTDIR}/${SAMPLE} \
    ${ASSEMBLY} \
    ${OUTDIR}/${SAMPLE}_hic_to_asm.bam

echo "[$(date)] YAHS scaffolding complete."

echo -e "\n Counting asm contig lengths after scaffolding..."
touch ${OUTDIR}/${SAMPLE}.seqlen.txt
awk '
/^>/ {
    if (NR > 1) {
        full_len = length(seq)
        gsub(/N/, "", seq)
        noN_len = length(seq)
        print header "\t" full_len "\t" noN_len
        total_full += full_len
        total_noN += noN_len
        count++
    }
    header = substr($0, 2)
    seq = ""
    next
}
{
    seq = seq $0
}
END {
    full_len = length(seq)
    gsub(/N/, "", seq)
    noN_len = length(seq)
    print header "\t" full_len "\t" noN_len
    total_full += full_len
    total_noN += noN_len
    count++
    print "TOTAL\t" count " sequences\t" total_full " bp (no-N: " total_noN " bp)"
}
' '${OUTDIR}/${SAMPLE}_scaffolds_final.fa' > '${OUTDIR}/${SAMPLE}.seqlen.txt'

echo "Counting YAHS-asm contigs length complete."

# Step 3: Generate Hi-C contact map for QC / manual curation (optional)
# This produces a .hic file viewable in Juicebox or PretextView
# Requires juicer_tools or the yahs juicer utilities
# Uncomment if you have juicer_tools available

echo "[$(date)] Generating contact map..."
echo "(Involves attempting juicer)"

samtools faidx ${OUTDIR}/${SAMPLE}_scaffolds_final.fa

yahs juicer pre \
    ${OUTDIR}/${SAMPLE}.bin \
    ${OUTDIR}/${SAMPLE}_scaffolds_final.agp \
    ${ASSEMBLY}.fai \
    > ${OUTDIR}/${SAMPLE}_juicer_pre_stdout.txt \
    2> ${OUTDIR}/${SAMPLE}_juicer_pre_stderr.txt


echo "yahs lines: $(wc -l < ${OUTDIR}/${SAMPLE}_juicer_pre.txt)"
echo "liftover lines: $(wc -l < ${OUTDIR}/${SAMPLE}_juicer_pre.liftover.tsv)"

# Step 2: Sort and filter
sort -k2,2d -k6,6d \
    -T ${OUTDIR} \
    --parallel=${THREADS} \
    -S 32G \
    ${OUTDIR}/${SAMPLE}_juicer_pre.txt \
    | awk 'NF' \
    > ${OUTDIR}/${SAMPLE}_juicer.txt

echo "juicer.txt lines: $(wc -l < ${OUTDIR}/${SAMPLE}_juicer.txt)"

# Step 3: Run juicer_tools with the fully-written liftover file
java -jar /home/user/projects/investigator/user/bin/juicer/CPU/common/juicer_tools.jar pre \
    --threads ${THREADS} \
    ${OUTDIR}/${SAMPLE}_juicer.txt \
    ${OUTDIR}/${SAMPLE}.hic \
    ${OUTDIR}/${SAMPLE}_juicer_pre.liftover.tsv

echo "  Scaffolded assembly : ${OUTDIR}/${SAMPLE}_scaffolds_final.fa"
echo "  AGP file            : ${OUTDIR}/${SAMPLE}_scaffolds_final.agp"
echo "  Contact matrix      : ${OUTDIR}/${SAMPLE}.bin"

# ------------------------------------------------------------------------------
# HIC CURATION PART II: Juicer
	# using juicer v1.2.2
# ------------------------------------------------------------------------------
# Step 1: Rename scaffolds — replace 'scaffold_' with 'chr' to avoid
# juicer_tools misinterpreting underscored names as chromosome combinations

echo "[$(date)] Renaming scaffolds..."
if [[ ! -f ${OUTDIR}/${SAMPLE}_scaffolds_final.fa ]]; then
    echo "ERROR: ${OUTDIR}/${SAMPLE}_scaffolds_final.fa not found" >&2
    exit 1
fi

sed -E 's/>scaffold_/>chr/' \
    ${OUTDIR}/${SAMPLE}_scaffolds_final.fa \
    > ${OUTDIR}/${SAMPLE}_scaffolds_renamed.fa
    
samtools faidx ${OUTDIR}/${SAMPLE}_scaffolds_renamed.fa

echo "Scaffold names after rename:"
cut -f1 ${OUTDIR}/${SAMPLE}_scaffolds_renamed.fa.fai | head -10

# Step 2: Index renamed scaffold assembly

echo "[$(date)] Indexing renamed scaffold assembly..."

bwa index ${OUTDIR}/${SAMPLE}_scaffolds_renamed.fa

# Step 3: Align Hi-C reads to renamed scaffolds
echo "[$(date)] Aligning Hi-C reads to scaffolds..."

bwa mem -5SP -t ${THREADS} \
    ${OUTDIR}/${SAMPLE}_scaffolds_renamed.fa \
    ${HIC_R1} ${HIC_R2} \
    > ${OUTDIR}/${SAMPLE}_hic_to_scaffolds.sam

echo "bwa exit: $?"
echo "sam lines: $(wc -l < ${OUTDIR}/${SAMPLE}_hic_to_scaffolds.sam)"

# Step 4: Convert SAM to BAM and sort by name
echo "[$(date)] Converting and sorting BAM..."

samtools view -@ ${THREADS} -b \
    ${OUTDIR}/${SAMPLE}_hic_to_scaffolds.sam \
    -o ${OUTDIR}/${SAMPLE}_hic_to_scaffolds_unsorted.bam

samtools sort -@ ${THREADS} -n \
    ${OUTDIR}/${SAMPLE}_hic_to_scaffolds_unsorted.bam \
    -o ${OUTDIR}/${SAMPLE}_hic_to_scaffolds.bam

echo "bam read count: $(samtools view -c ${OUTDIR}/${SAMPLE}_hic_to_scaffolds.bam)"

# Step 5: Filter and convert (SORT!) to juicer format

echo "[$(date)] Filtering alignments..."

samtools view -@ ${THREADS} -F 2316 \
    ${OUTDIR}/${SAMPLE}_hic_to_scaffolds.bam \
    > ${OUTDIR}/${SAMPLE}_hic_to_scaffolds_filtered.sam

#ALT
awk 'BEGIN{OFS="\t"}
!/^@/ {
    chr1 = $3
    pos1 = $4
    chr2 = ($7 == "=" ? $3 : $7)
    pos2 = $8
    str1 = and($2, 16) ? 16 : 0
    str2 = and($2, 32) ? 16 : 0
    print str1, chr1, pos1, 0, str2, chr2, pos2, 1
}' ${OUTDIR}/${SAMPLE}_hic_to_scaffolds_filtered.sam \
| awk 'BEGIN{OFS="\t"} {
    if ($2 > $6)
        print $5, $6, $7, $8, $1, $2, $3, $4
    else
        print
}' \
| sort -k2,2d -k6,6d \
> ${OUTDIR}/${SAMPLE}_juicer_almostsorted.txt

#test sorting
awk '{print $2"_"$6}' ${OUTDIR}/${SAMPLE}_juicer_sorted2.txt | uniq | sort | uniq -d
echo "march on comrade. happiness is mandatory"

awk 'NR==FNR{rank[$1]=NR; next}
     {print rank[$2], rank[$6], $0}' chrom.sizes ${OUTDIR}/${SAMPLE}_juicer_almostsorted.txt \
| sort -k1,1n -k2,2n \
| cut -d' ' -f3- > ${OUTDIR}/${SAMPLE}_juicer_sorted.txt

cut -f1,2 $OUTDIR/${SAMPLE}_scaffolds_renamed.fa.fai \
    > $OUTDIR/${SAMPLE}.chrom.sizes

# Step 7: Run juicer_tools pre
echo "[$(date)] Running juicer_tools pre..."

java -jar /home/user/projects/investigator/user/bin/juicer/CPU/common/juicer_tools.jar pre \
    ${OUTDIR}/${SAMPLE}_juicer_sorted.txt \
    ${OUTDIR}/${SAMPLE}.hic \
    ${OUTDIR}/${SAMPLE}_chrom.sizes
echo "hic size: $(ls -lh ${OUTDIR}/${SAMPLE}.hic)"

# Step 9: Remove intermediary files
echo "[$(date)] Cleaning up intermediary files..."

rm ${OUTDIR}/${SAMPLE}_hic_to_scaffolds.sam
rm ${OUTDIR}/${SAMPLE}_hic_to_scaffolds_unsorted.bam
rm ${OUTDIR}/${SAMPLE}_hic_to_scaffolds_filtered.sam
rm ${OUTDIR}/${SAMPLE}_juicer_sorted.txt

echo "Remaining files in ${OUTDIR}:"
ls -lAtrho ${OUTDIR}

# Juicer DONE
bash rename_juicer_contigs.sh

echo "[$(date)] Pipeline complete. Outputs in: ${OUTDIR}/"
echo "  Scaffolded assembly : ${OUTDIR}/${SAMPLE}_scaffolds_final.fa"
echo "  Renamed assembly    : ${OUTDIR}/${SAMPLE}_scaffolds_renamed.fa"
echo "  AGP file            : ${OUTDIR}/${SAMPLE}_scaffolds_final.agp"
echo "  Contact map         : ${OUTDIR}/${SAMPLE}.hic"
