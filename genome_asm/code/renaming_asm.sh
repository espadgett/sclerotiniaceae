#!/bin/bash

#renaming_asm.sh reorders the hifiasm-output assembly scaffolds in harmony with the B. cinerea contigs, determined by ortholog synteny

# ---------------------------------------------------------------------
#Renaming Streptotinia caulophylli
# - SORTING ORIGINAL CONTIGS
# - NAMED BY DESC LENGTH
# ---------------------------------------------------------------------

sc_asm="/home/padgett3/scratch/streptotinia_caulophylli/hic_asm/SC_hifihic_020326.asm.hic.p_ctg.fa" #input genome asm (HiC + HiFi, assembled by hifiasm)
sorted_sc_asm="/home/padgett3/scratch/streptotinia_caulophylli/hic_asm/SC_hifihic_200426.sizesorted.fa"
sorted_renamed_sc_asm="/home/padgett3/scratch/streptotinia_caulophylli/hic_asm//SC_hifihic_200426.sizesorted.renamed.fa"
seqlen="/home/padgett3/scratch/streptotinia_caulophylli/hic_asm/sc_contig_len.sorted.renamed.txt"

cd "/home/padgett3/scratch/streptotinia_caulophylli/hic_asm/"
#sort contigs by size
awk '/^>/ {printf("%s%s\t", (N>0?"\n":""), $0); N++; next;} {printf("%s", $0);} END {printf("\n");}' "$sc_asm" \
| awk -F '\t' '{printf("%d\t%s\n", length($2), $0);}' \
| sort -k1,1nr \
| cut -f2- \
| tr "\t" "\n" \
> "$sorted_sc_asm"

#rename contigs
awk '
    /^>/ {
        i++
        print ">contig" i
        next
    }
    { print }
' "$sorted_sc_asm" > "$sorted_renamed_sc_asm"

#if you like to remove intermediary files within scripts, otherwise hash out
rm "$sorted_sc_asm"

#count seq lengths
awk '
    /^>/ {
        # If this is not the first record, process the previous one
        if (NR > 1) {
            full_len = length(seq)
            gsub(/N/, "", seq)
            noN_len = length(seq)
            print header "\t" full_len "\t" noN_len
        }
        header = substr($0, 2)   #get header sequence name
        seq = ""                 
        next
    }
    {
        seq = seq $0             
    }
    END {
        # Process the final record
        full_len = length(seq)
        gsub(/N/, "", seq)
        noN_len = length(seq)
        print header "\t" full_len "\t" noN_len
    }
' "$sorted_renamed_sc_asm" > "$seqlen"

# ---------------------------------------------------------------------
#Performed synteny analysis > using GENESPACE-alignments to inform name changes in harmony with reference B. cinerea
# - NAMED BY REFERENCE SYNTENY
# - REORDER BY DESCENDING CONTIG/SCAFFOLD NUMBER
# ---------------------------------------------------------------------
INPUT="/home/padgett3/scratch/streptotinia_caulophylli/hic_asm/SC_hifihic_200426.sizesorted.renamed.fa"
OUTPUT="Scaulophylli_270426.fa"
SEQLEN="Scaulophylli_270426_contiglen.txt"

echo "Input : $INPUT"
echo "Output: $OUTPUT"
echo ""

# Use a temporary file to avoid overwriting input if output == input
TMP1="TMPFILE1"
TMP2="TMPFILE2"

#rename according to manual GENESPACE curation with hifi-only asm of S. caulophylli
awk '
BEGIN {
    # Mapping: original -> new name
    map["contig3"]  = "contig5"
    map["contig4"]  = "contig3"
    map["contig5"]  = "contig7"
    map["contig7"]  = "contig9"
    map["contig8"]  = "contig10"
    map["contig9"]  = "contig8"
    map["contig10"] = "contig12"
    map["contig12"] = "contig4"
    map["contig13"] = "contig15"
    map["contig14"] = "contig13"
    map["contig15"] = "contig16"
    map["contig16"] = "contig14"

    contig_count = 0
    renamed      = 0
    unchanged    = 0
}

/^>/ {
    contig_count++

    if (contig_count <= 16) {
        # Split header into name + optional description
        header = substr($0, 2)          # strip leading ">"
        split(header, parts, /[ \t]/, seps)
        name = parts[1]

        if (name in map) {
            new_name = map[name]
            # Rebuild the rest of the header (description, if any)
            rest = ""
            for (i = 2; i <= length(parts); i++) {
                rest = rest seps[i-1] parts[i]
            }
            print ">" new_name rest
            printf "  [contig %2d] %s -> %s\n", contig_count, name, new_name > "/dev/stderr"
            renamed++
        } else {
            print $0
            printf "  [contig %2d] %s (no mapping, kept as-is)\n", contig_count, name > "/dev/stderr"
            unchanged++
        }
        next
    }

    # Beyond first 16: write unchanged
    print $0
    next
}

# Sequence lines — always write unchanged
{ print }

END {
    printf "\nDone. %d contig(s) renamed, %d left unchanged (out of first 16). Total records: %d.\n", \
        renamed, unchanged, contig_count > "/dev/stderr"
}
' "$INPUT" > "$TMP1"

#reorder contigs in descending (numeric) name order 
awk '
/^>/ {
    if (seqname != "") {
        # Store completed record
        sequences[seqname] = seq
        headers[seqname] = header
    }
    header = $0
    seqname = substr($0, 2)
    split(seqname, parts, /[ \t]/)
    seqname = parts[1]
    seq = ""
    next
}
{
    seq = seq $0 "\n"
}
END {
    # Store last record
    sequences[seqname] = seq
    headers[seqname] = header
 
    # Collect all names into an array for sorting
    n = 0
    for (name in sequences) {
        names[++n] = name
    }
 
    # Sort numerically by the trailing number in the contig name
    # Use bubble sort (portable, no gawk needed)
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n - i; j++) {
            # Extract numeric suffix for comparison
            a = names[j];   sub(/^[^0-9]*/, "", a); numa = a + 0
            b = names[j+1]; sub(/^[^0-9]*/, "", b); numb = b + 0
            if (numa > numb) {
                tmp = names[j]
                names[j] = names[j+1]
                names[j+1] = tmp
            }
        }
    }
 
    # Print in sorted order
    for (i = 1; i <= n; i++) {
        print headers[names[i]]
        printf "%s", sequences[names[i]]
    }
}
' "$TMP1" > "$TMP2"
 
mv "$TMP2" "$OUTPUT"
rm -f "$TMP1"
 
echo "Reordering complete."
echo ""

#count seq lengths
awk '
    /^>/ {
        # If this is not the first record, process the previous one
        if (NR > 1) {
            full_len = length(seq)
            gsub(/N/, "", seq)
            noN_len = length(seq)
            print header "\t" full_len "\t" noN_len
        }
        header = substr($0, 2)   #get header sequence name
        seq = ""                 
        next
    }
    {
        seq = seq $0             
    }
    END {
        # Process the final record
        full_len = length(seq)
        gsub(/N/, "", seq)
        noN_len = length(seq)
        print header "\t" full_len "\t" noN_len
    }
' "$OUTPUT" > "$SEQLEN"