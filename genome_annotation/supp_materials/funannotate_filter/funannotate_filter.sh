#!/bin/bash
#cd to funannotate output dir
echo "Funannotate final outputs:"
echo "$(ls -F)"

outdir="purgedcontigs"
samplename="Scaulophylli"

echo "Files to directly copy over, without filtering:"
echo -e "\t Gene2Products.need-curating.txt"
echo -e "\t Gene2Products.must-fix.txt"
echo -e "\t Gene2Products.new-names-passed.txt"
echo "These have been manually curated by Esme Padgett."
echo "Note: filters are between the current and older version of fasta - where many contigs were removed during duplication purging.
Don't worry too much about it, I just didn't want to run Funannotate again :)"

echo -e "\n Initiating file filtering to keep only the following (non-duplicate) sequences:"

echo -e "\nSt_cau_sc_1
St_cau_sc_2
St_cau_sc_3
St_cau_sc_4
St_cau_sc_5
St_cau_sc_6
St_cau_sc_7
St_cau_sc_8
St_cau_sc_9
St_cau_sc_10
St_cau_sc_11
St_cau_sc_12
St_cau_sc_13
St_cau_sc_14
St_cau_sc_15
St_cau_sc_16
St_cau_sc_17
St_cau_sc_19
St_cau_sc_29
St_cau_sc_30
St_cau_sc_33
St_cau_sc_51"

KEEP_IDS="St_cau_sc_1 St_cau_sc_2 St_cau_sc_3 St_cau_sc_4 St_cau_sc_5 \
St_cau_sc_6 St_cau_sc_7 St_cau_sc_8 St_cau_sc_9 St_cau_sc_10 St_cau_sc_11 \
St_cau_sc_12 St_cau_sc_13 St_cau_sc_14 St_cau_sc_15 St_cau_sc_16 St_cau_sc_17 \
St_cau_sc_19 St_cau_sc_29 St_cau_sc_30 St_cau_sc_33 St_cau_sc_51"

KEEP_FILE=keep_file.txt
echo "$KEEP_IDS" | tr ' ' '\n' > "$KEEP_FILE"

#gbk
awk 'BEGIN { while ((getline id < ARGV[1]) > 0) keep[id]=1; delete ARGV[1] }
     /^LOCUS/ { record_keep = ($2 in keep) }
     record_keep { print }
     /^\/\/$/ { record_keep = 0 }
    ' "$KEEP_FILE" *.gbk > "${outdir}/${samplename}.gbk"
    
#annotations.txt
awk 'BEGIN { while ((getline id < ARGV[1]) > 0) keep[id]=1; delete ARGV[1] }
     NR==1 { print; next }
     $4 in keep { print }
    ' "$KEEP_FILE" *.annotations.txt > "${outdir}/${samplename}.annotations.txt"

#scaffolds.fa
awk 'BEGIN { while ((getline id < ARGV[1]) > 0) keep[id]=1; delete ARGV[1] }
     /^>/ { seq_id = substr($1, 2); record_keep = (seq_id in keep) }
     record_keep { print }
    ' "$KEEP_FILE" *.scaffolds.fa > "${outdir}/${samplename}.scaffolds.fa"

#tbl
kept=0; skipped=0; current_block=""; current_keep=0
> "${outdir}/${samplename}.tbl"

while IFS= read -r line || [ -n "$line" ]; do
    if [[ "$line" == ">Feature "* ]]; then
        if [ -n "$current_block" ]; then
            if [ "$current_keep" -eq 1 ]; then
                printf '%s\n' "$current_block" >> "${outdir}/${samplename}.tbl"; ((kept++))
            else
                ((skipped++))
            fi
        fi
        seq_id=$(echo "$line" | awk '{print $2}')
        grep -qx "$seq_id" "$KEEP_FILE" && current_keep=1 || current_keep=0
        current_block="$line"
    else
        current_block="${current_block}"$'\n'"${line}"
    fi
done < *.tbl

# Handle last block
if [ -n "$current_block" ]; then
    if [ "$current_keep" -eq 1 ]; then
        printf '%s\n' "$current_block" >> "${outdir}/${samplename}.tbl"; ((kept++))
    else
        ((skipped++))
    fi
fi
echo "TBL: Kept $kept, Skipped $skipped"

#sqn
tbl2asn -i "${outdir}/${samplename}.scaffolds.fa" -f "${outdir}/${samplename}.tbl" -o "${outdir}/${samplename}.sqn" -V vb

#cds-transcripts.fa
awk 'NR>1 { print $2 }' "${outdir}/${samplename}.annotations.txt" > "${outdir}/keep_transcripts.txt"
awk 'BEGIN { while ((getline id < ARGV[1]) > 0) keep[id]=1; delete ARGV[1] }
     /^>/ { seq_id = substr($1, 2); record_keep = (seq_id in keep) }
     record_keep { print }
    ' "${outdir}/keep_transcripts.txt" *cds-transcripts.fa > "${outdir}/${samplename}.cds-transcripts.fa"


#proteins.fa
awk 'BEGIN { while ((getline id < ARGV[1]) > 0) keep[id]=1; delete ARGV[1] }
     /^>/ { seq_id = substr($1, 2); record_keep = (seq_id in keep) }
     record_keep { print }
    ' "${outdir}/keep_transcripts.txt" *proteins.fa > "${outdir}/${samplename}.proteins.fa"
    
#mrna-transcripts.fa
awk 'BEGIN { while ((getline id < ARGV[1]) > 0) keep[id]=1; delete ARGV[1] }
     /^>/ { seq_id = substr($1, 2); record_keep = (seq_id in keep) }
     record_keep { print }
    ' "${outdir}/keep_transcripts.txt" *mrna-transcripts.fa > "${outdir}/${samplename}.mrna-transcripts.fa"

#gff3
awk 'BEGIN { while ((getline id < ARGV[1]) > 0) keep[id]=1; delete ARGV[1] }
     NR==1 { print; next }
     $1 in keep { print }
    ' "$KEEP_FILE" *gff3 > "${outdir}/${samplename}.gff3"

#agp
awk 'BEGIN { while ((getline id < ARGV[1]) > 0) keep[id]=1; delete ARGV[1] }
     NR==1 { print; next }
     $1 in keep { print }
    ' "$KEEP_FILE" *agp > "${outdir}/${samplename}.agp"