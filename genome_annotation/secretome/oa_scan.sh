#!/bin/bash

module load cd-hit
cd /home/user/scratch/streptotinia_caulophylli/hic_asm/oa_scan

OAPROTEINS="/home/user/scratch/streptotinia_caulophylli/hic_asm/oa_scan/oa_relevant_proteins2.fa"
FASTA1="/home/user/scratch/streptotinia_caulophylli/genespace/run280726/rawGenomes/bcin/bcinerea.protein.fa"
FASTA2="/home/user/scratch/streptotinia_caulophylli/genespace/run280726/rawGenomes/scau/Scaulophylli.proteins.fa"
FASTA3="/home/user/scratch/streptotinia_caulophylli/genespace/run280726/rawGenomes/smin/renamed.Sminor_Sminor.proteins.fa"
FASTA4="/home/user/scratch/streptotinia_caulophylli/genespace/run280726/rawGenomes/sscl/ssclerotiorum.protein.fa"

OUTDIR="oascan_rerun/oascan_bcin"
mkdir -p $OUTDIR
COMBINED="$OUTDIR/combined.fasta"
CLUSTERED="$OUTDIR/clustered"
IDENTITY="0.70"
echo "Combine FASTA files -> $COMBINED ..."
cat $OAPROTEINS $FASTA1 > $COMBINED
echo "Running CD-HIT (identity=$IDENTITY, word_len=$WORD_LEN) against B.cinerea ..."
cd-hit \
  -i "$COMBINED" \
  -o "$CLUSTERED" \
  -c "$IDENTITY" \
  -M 800 \
  -d 0 \
  -g 1
echo "Finished CD-HIT for B.cinerea OA gene scan."

OUTDIR="oascan_rerun/oascan_scau"
mkdir -p $OUTDIR
COMBINED="$OUTDIR/combined.fasta"
CLUSTERED="$OUTDIR/clustered"
IDENTITY="0.70"
echo "Combine FASTA files -> $COMBINED ..."
cat $OAPROTEINS $FASTA2 > $COMBINED
echo "Running CD-HIT (identity=$IDENTITY, word_len=$WORD_LEN) against St. caulophylli ..."
cd-hit \
  -i "$COMBINED" \
  -o "$CLUSTERED" \
  -c "$IDENTITY" \
  -M 800 \
  -d 0 \
  -g 1  
echo "Finished CD-HIT for St. caulophylli OA gene scan."

OUTDIR="oascan_rerun/oascan_smin"
mkdir -p $OUTDIR
COMBINED="$OUTDIR/combined.fasta"
CLUSTERED="$OUTDIR/clustered"
IDENTITY="0.70"
echo "Combine FASTA files -> $COMBINED ..."
cat $OAPROTEINS $FASTA3 > $COMBINED
echo "Running CD-HIT (identity=$IDENTITY, word_len=$WORD_LEN) against S.minor ..."
cd-hit \
  -i "$COMBINED" \
  -o "$CLUSTERED" \
  -c "$IDENTITY" \
  -M 800 \
  -d 0 \
  -g 1    
echo "Finished CD-HIT for S. minor OA gene scan."

OUTDIR="oascan_rerun/oascan_sscl"
mkdir -p $OUTDIR
COMBINED="$OUTDIR/combined.fasta"
CLUSTERED="$OUTDIR/clustered"
IDENTITY="0.70"
echo "Combine FASTA files -> $COMBINED ..."
cat $OAPROTEINS $FASTA4 > $COMBINED
echo "Running CD-HIT (identity=$IDENTITY, word_len=$WORD_LEN) against S.sclerotiorum ..."
cd-hit \
  -i "$COMBINED" \
  -o "$CLUSTERED" \
  -c "$IDENTITY" \
  -M 800 \
  -d 0 \
  -g 1    
echo "Finished CD-HIT for S. sclerotiorum OA gene scan."  

#===============================
echo -e "\nBLASTP FOR A FULL SCAN for SMINOR/STCAU"
#===============================
#blastp
module load blast/2.2.26 blast+/2.17.0

echo "Running BLASTP for Sminor for 2 OA-related GOI"
echo "[1/2] Constructing BLASTP db from SM proteome"
makeblastdb -in $FASTA3 -dbtype prot -out sm_proteome_db
echo "[2/2] Running BLASTP on SM"
blastp -query $OAPROTEINS \
	-db sm_proteome_db \
	-out "oascan_rerun/oascan_smin/sm_blastp_result.tsv" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
  -evalue 1e-5 \
  -max_target_seqs 100
echo "Finished BLASTP on SM"

echo "Running BLASTP for Stcaulophylli for 2 OA-related GOI"
echo "[1/2] Constructing BLASTP db from SC proteome"
makeblastdb -in $FASTA2 -dbtype prot -out sc_proteome_db
echo "[2/2] Running BLASTP on SC"
blastp -query $OAPROTEINS \
	-db sc_proteome_db \
	-out "oascan_rerun/oascan_scau/sc_blastp_result.tsv" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
  -evalue 1e-5 \
  -max_target_seqs 100
echo "Finished BLASTP on SC"

echo "Running BLASTP for Bcinerea for 2 OA-related GOI"
echo "[1/2] Constructing BLASTP db from SC proteome"
makeblastdb -in $FASTA1 -dbtype prot -out bc_proteome_db
echo "[2/2] Running BLASTP on SC"
blastp -query $OAPROTEINS \
	-db bc_proteome_db \
	-out "oascan_rerun/oascan_bcin/bc_blastp_result.tsv" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
  -evalue 1e-5 \
  -max_target_seqs 100
echo "Finished BLASTP on BC"

echo "Running BLASTP for Ssclerotiorum for 3 OA-related GOI"
echo "[1/2] Constructing BLASTP db from SC proteome"
makeblastdb -in $FASTA4 -dbtype prot -out ss_proteome_db
echo "[2/2] Running BLASTP on SC"
blastp -query $OAPROTEINS \
	-db ss_proteome_db \
	-out "oascan_rerun/oascan_sscl/ss_blastp_result.tsv" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
  -evalue 1e-5 \
  -max_target_seqs 100
echo "Finished BLASTP on SS"