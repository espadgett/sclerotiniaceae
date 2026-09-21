#!/bin/bash

#rename_juicer_contigs.sh removes low-confidence contigs/keeps only high-support contigs after Juicer (desktop) visualization

in_fasta="/home/padgett3/scratch/streptotinia_caulophylli/hic_asm/Scaulophylli_280426.fasta"
out_fasta="/home/padgett3/scratch/streptotinia_caulophylli/hic_asm/Scaulophylli_200526.fasta"
keep_contigs=(
    contig1
    contig2
    contig3
    contig4
    contig5
    contig6
    contig7
    contig8
    contig9
    contig10
    contig11
    contig12
    contig13
    contig14
    contig15
    contig16
    contig17
    contig19
    contig30
    contig31
    contig34
    contig66
)
chr_count=16 #number of scaffolds to rename as chr-level sequences
# ------------------------------------------------------------------------------
#safety check:
if [[ ! -f "$in_fasta" ]]; then
    echo "ERROR: Input FASTA not found: $in_fasta. Aborting." >&2
    exit 1
fi

keep=$(printf '%s\n' "${keep_contigs[@]}")
totalseq=${#keep_contigs[@]}
 
awk -v keep_list="$keep" \
    -v chr="$chr_count" \
    -v total="$totalseq" \
'
BEGIN {
    # Load the ordered keep list into an array
    n = split(keep_list, order, "\n")
    for (i = 1; i <= n; i++) {
        wanted[order[i]] = i          # value = position (1-based)
    }
    printing = 0
}

/^>/ {
    printing = 0
    # Extract the contig name (word after ">", strip any description)
    split(substr($0, 2), parts, " ")
    name = parts[1]
 
    if (name in wanted) {
        pos = wanted[name]            # 1-based position in keep list
 
        if (pos <= chr) {
            new_name = "chr" pos
        } else {
            scaffold_num = pos - chr
            new_name = "scaffold" scaffold_num
        }
 
        print ">" new_name
        printing = 1
    }
    next
}

printing { print }
' "$in_fasta" > "$out_fasta"

kept=$(grep -c "^>" "$out_fasta" 2>/dev/null || true )
echo "  Input  : $in_fasta"
echo "  Output : $out_fasta"
echo "  Sequences written: $kept / $totalseq input."
echo -e "\n Sequence filtering done."

echo -e "\n Summarizing new fasta sequence names and lengths..."
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
' "$out_fasta" > "${out_fasta}.seqlen.txt"