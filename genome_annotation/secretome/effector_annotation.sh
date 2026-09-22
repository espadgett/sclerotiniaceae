#=======================================================
# SignalP v6 (I ran locally)
#=======================================================
conda activate signalp6_env
cd /Users/user/bioinformatics_programs/signalp-v.6/signalp6_fast
python -m pip install --upgrade build && python -m build signalp-6-package #optional - for clean wheel build
SIGNALP_DIR=$(python -c "import signalp, os; print(os.path.dirname(signalp.__file__))")  
   cp -r signalp-6-package/models/* "$SIGNALP_DIR/model_weights/"
signalp6 --help

signalp6 --fastafile secretome/Scaulophylli.proteins.fa \
--organism eukarya \
--output_dir /Users/user/postgrad/lab_docx/fungal_asm_data/secretome/signalp_out/scaulophylli_signalp \
--format txt \
--mode fast \
--model_dir /opt/anaconda3/envs/signalp6_env/lib/python3.10/site-packages/signalp/model_weights

#=======================================================
# Phobius (web server) - THIS IS OPTIONAL, SIGNALP THUS FAR WORKED WELL ENOUGH TO REMOVE NEED FOR PHOBIUS
#=======================================================
echo "parse Phobius output for THM detection"
awk 'NR>1 { if ($4 == "i") i++; else more++; } END { print "i_only:", i; print "more_than_i:", more }' phobius_shortout.txt

echo "Rewrite fasta with only verified signal peptides"
awk 'NR>1 && $4=="i"' phobius_shortout.txt > phobius_corrections.txt  
#[nano] column names back in
echo "Cross-reference phobius_corrections.txt with SignalP 3processed_entries.fasta" 
awk 'NR==FNR { keep[$1]=1; next } /^>/ { id=$1; sub(/^>/,"",id) } !drop[id]' phobius_corrections.txt processed_entries.fasta > thmfilt_secretory.fasta

#=======================================================
# NetGPI v1.0 (web server)
#=======================================================
#input <<thmfilt_secretory.fasta>> to web server

#output_mature.fasta identifies GPI-anchor PTModified proteins
#those fasta IDS will be filtered against the previously PTM-filtered secretory fasta output by SignalP/Phobius

awk 'NR==FNR {
        if ($0 ~ /^>/) {
            id=$1
            sub(/^>/,"",id)
            drop[id]=1
        }
        next
     }
     /^>/ {
        id=$1
        sub(/^>/,"",id)
     }
     !drop[id]' secretome/netgpi_out/netgpi/output_mature.fasta \
     secretome/phobius_out/phobius/thmfilt_secretory.fasta \
     > gpifilt_secretory.fasta
     
#=======================================================
# EffectorP v3.0 (web server)
#=======================================================
#input: g<<pifilt_secretory.fasta>> to web server


#=======================================================
# PHI (Pathogen-Host Interaction) Database (web server)
	#cp effectorp output files to HPC
	#PHI-DB is used to create high-confidence effector candidates, but stopping at EffectorP produces quality effector protein estimations
#=======================================================
# MODULE LOAD
module load StdEnv/2023 blast+/2.17.0 
#=======================================================
#make PHI database
mkdir phi-base-raw
cd phi-base-raw
wget https://raw.githubusercontent.com/PHI-base/data/master/releases/phi-base_current.fas
makeblastdb -in phi-base_current.fas -dbtype prot -out phibase_db

psiblast \
  -query effectorpred.fasta \
  -db /home/user/scratch/streptotinia_caulophylli/hic_asm/annotations/phib-blast/phi-base-raw/phibase_db \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" \
  -out phidb_positive.tsv \
  -evalue 1e-5 \
  -num_threads 8

echo "Number of proteins with positive PHI database hits (1e-05): $(cut -f1 phidb_positive.tsv | sort -u | wc -l)"

echo "Writing out positive results to phidb_positive_ids.tsv... $(cut -f1 phidb_positive.tsv | sort -u > phidb_positive_ids.tsv)"

echo "Writing out a multi-fasta file of the PHI-positive results..."
awk 'NR==FNR{ids[$1]=1; next} 
     /^>/{id=substr($1,2); keep=(id in ids)} 
     keep{print}' \
		phidb_positive_ids.tsv effectorpred.fasta > phi_pos.fasta