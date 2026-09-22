#!/bin/bash

#=======================================================
# ORTHOLOG ALIGNMENT
#=======================================================
# MODULE LOAD
module load StdEnv/2023 gcc/14.3 python diamond/2.1.11 r/4.5.0
source ~/envs/genespace_env/bin/activate

export PATH=$PATH:/home/user/scratch/streptotinia_caulophylli/genespace/OrthoFinder_source
echo "Confirm Orthofinder:"
orthofinder -h

# FILES/PATHS

	#input list of genomes to extract orthologs
echo "make file of files of protein fasta for ORTHOFINDER "
input_paths=(
"/home/user/scratch/streptotinia_caulophylli/genespace/run280726/rawGenomes/bcin/bcinerea.protein.fa"
"/home/user/scratch/streptotinia_caulophylli/genespace/run280726/rawGenomes/scau/Scaulophylli.proteins.fa"
"/home/user/scratch/streptotinia_caulophylli/genespace/run280726/rawGenomes/smin/renamed.Sminor_Sminor.proteins.fa"
"/home/user/scratch/streptotinia_caulophylli/genespace/run280726/rawGenomes/sscl/ssclerotiorum.protein.fa"
"/home/user/scratch/streptotinia_caulophylli/ref_genomes/Ncrassa/protein.faa"
)

printf '%s\n' "${input_paths[@]}" > hic_asm/phylogenetic_tree/fof_protein.txt

cd hic_asm/phylogenetic_tree
echo "verify input file of files exists"
file fof_protein.txt

mkdir -p proteomes
while read -r f; do
    name=$(basename "$f" | sed 's/\.[^.]*$//')
    ln -sf "$(realpath "$f")" "proteomes/${name}.fa"
done < fof_protein.txt

# Orthofinder for extraction of orthologs
echo -e "\nInitating OrthoFinder on fungal proteomes"
orthofinder -f proteomes -og -t 8 -a 8
echo "Finished OrthoFinder."

RESULTS_DIR=$(find "OrthoFinder" -maxdepth 1 -type d -name "Results_*" | sort | tail -n 1)
RESULTS_DIR=/scratch/user/streptotinia_caulophylli/hic_asm/phylogenetic_tree/proteomes/OrthoFinder/Results_Jul30_1/
[[ -n "$RESULTS_DIR" ]] || { echo "ERROR: could not locate OrthoFinder Results_* directory." >&2; exit 1; }
echo "    OrthoFinder results: $RESULTS_DIR"

# MAFFT alignment of orthologs
module load mafft/7.526
mkdir -p aligned
for f in "$RESULTS_DIR"/Single_Copy_Orthologue_Sequences/*.fa; do
    og=$(basename "$f" .fa)
    mafft --auto --quiet "$f" > "aligned/${og}_aligned.fasta"
done

# remove poorly aligned regions
module load trimal/1.5.1
mkdir -p trimmed

find aligned -name "*_aligned.fasta" | \
    xargs -P 8 -I{} bash -c \
    'base=$(basename "{}" _aligned.fasta); trimal -in "{}" -out "trimmed/${base}_trimmed.fasta" -automated1'
    
#=======================================================
# TREE CONSTRUCTION
#=======================================================
#MODULE LOAD
module load java/17.0.6

#IQTREE3 and ASTRAL installed
iqtree3 --version
ASTRAL/bin/astral --help


#=======================================================
iqtree3 -S trimmed/ \
 -B 1000 --prefix gene_trees -T AUTO -st AA 
 
ASTRAL/bin/astral \
-i /home/user/scratch/streptotinia_caulophylli/hic_asm/phylogenetic_tree/gene_trees.treefile \
-o species_tree_astral.tre

python /Users/user/trees-to-networks-tutorial/tools/make_nexus_treeset.py \
    --alignments-dir /Users/user/postgrad/lab_docx/fungal_asm_data/phylogenetic/trimmed \
    --treefile /Users/user/postgrad/lab_docx/fungal_asm_data/phylogenetic/gene_trees.treefile \
    --out calculated.nex \
    --outgroup 

#visualization in iTOL web interface
   