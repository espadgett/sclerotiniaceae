#!/usr/bin/env bash
# filter_stats.sh — Extract assembly and annotation statistics from genome files.
#
# Usage:
#   bash filter_stats.sh \
#       --fasta   Scaulophylli.scaffolds.fa --gff     Scaulophylli.gff3 --cds     Scaulophylli.cds-transcripts.fa --protein Scaulophylli.proteins.fa --output Scaulophylli_200526.stats.json
#
# Requirements: awk, grep, sort, sed (standard on Linux/macOS)

set -euo pipefail

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
FASTA="" GFF="" CDS="" PROTEIN="" OUTPUT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta)   FASTA="$2";   shift 2 ;;
        --gff)     GFF="$2";     shift 2 ;;
        --cds)     CDS="$2";     shift 2 ;;
        --protein) PROTEIN="$2"; shift 2 ;;
        --output)  OUTPUT="$2";  shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$FASTA" || -z "$GFF" ]] && {
    echo "Usage: $0 --fasta <genome.fasta> --gff <annotation.gff3> [--cds <cds.fasta>] [--protein <protein.fasta>] [--output <out.json>]" >&2
    exit 1
}

for f in "$FASTA" "$GFF"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f" >&2; exit 1; }
done

# ---------------------------------------------------------------------------
# Helper: linearise a FASTA file into "id\tsequence" TSV (one line per record)
# ---------------------------------------------------------------------------
linearise_fasta() {
    awk '/^>/{if(seq)print id"\t"seq; id=substr($0,2); seq=""; next} {seq=seq$0} END{if(seq)print id"\t"seq}' "$1"
}

# ---------------------------------------------------------------------------
# ASSEMBLY STATS
# ---------------------------------------------------------------------------
echo "Calculating assembly statistics..." >&2

read ASM_NUM_CONTIGS ASM_TOTAL_LEN ASM_GC < <(
    linearise_fasta "$FASTA" | awk '
    {
        seq = toupper($2)
        len = length(seq)
        num_contigs++
        total_len += len
        gc = gsub(/[GC]/, "", seq)
        at = gsub(/[AT]/, "", seq)
        total_gc += gc
        total_at += at
    }
    END {
        gc_content = (total_gc+total_at > 0) ? total_gc/(total_gc+total_at)*100 : 0
        printf "%d %d %.4f\n", num_contigs, total_len, gc_content
    }'
)

# N50/L50/N90/L90 — sort lengths descending, walk cumulative sum
read ASM_N50 ASM_L50 ASM_N90 ASM_L90 < <(
    linearise_fasta "$FASTA" | awk '{print length($2)}' | sort -rn | awk -v total="$ASM_TOTAL_LEN" '
    BEGIN { target50=total*0.5; target90=total*0.9; run=0; n50=0; l50=0; n90=0; l90=0 }
    {
        run += $1; i++
        if (!n50 && run >= target50) { n50=$1; l50=i }
        if (!n90 && run >= target90) { n90=$1; l90=i }
    }
    END { printf "%d %d %d %d\n", n50, l50, n90, l90 }'
)

ASM_MEAN_LEN=$(awk "BEGIN{printf \"%.2f\", $ASM_TOTAL_LEN/$ASM_NUM_CONTIGS}")
ASM_GC=$(awk "BEGIN{printf \"%.2f\", $ASM_GC}")

# ---------------------------------------------------------------------------
# ANNOTATION STATS — gene-level
# ---------------------------------------------------------------------------
echo "Calculating annotation statistics..." >&2

# Strip comments, keep data lines
GFF_DATA=$(grep -v '^#' "$GFF" | awk -F'\t' 'NF==9')

read ANN_GENES ANN_MRNA ANN_TRNA ANN_RRNA ANN_NCRNA < <(
    echo "$GFF_DATA" | awk -F'\t' '
        $3=="gene"  { genes++ }
        $3=="mRNA"  { mrna++  }
        $3=="tRNA"  { trna++  }
        $3=="rRNA"  { rrna++  }
        $3=="ncRNA"||$3=="lnc_RNA"||$3=="snoRNA"||$3=="snRNA"||$3=="miRNA" { ncrna++ }
        END { printf "%d %d %d %d %d\n", genes+0, mrna+0, trna+0, rrna+0, ncrna+0 }')

# common_name: genes with a Name= attribute that has a non-empty value
ANN_COMMON_NAME=$(echo "$GFF_DATA" | awk '$3=="gene" && $9~/Name=/' | grep -cP 'Name=[^;]+' || true)

# avg gene length
ANN_AVG_GENE_LEN=$(echo "$GFF_DATA" | awk '
    $3=="gene" { len += ($5-$4+1); n++ }
    END { printf "%.2f", (n>0 ? len/n : 0) }')

#exon metrics
read -r MULTI SINGLE AVG_EXON TOTAL_EXONS < <(awk '
/^#/ || /^$/ { next }      # skip comments and blank lines
$3 != "exon"  { next }     # only process exon features
 
{
    start = $4 + 0
    end   = $5 + 0
 
    # Extract Parent= value (funannotate format: "Parent=TRANSCRIPT_ID;")
    parent = ""
    n = split($9, attrs, ";")
    for (i = 1; i <= n; i++) {
        if (attrs[i] ~ /^Parent=/) {
            parent = substr(attrs[i], 8)   # strip leading "Parent="
            break
        }
    }
    if (parent == "") next
 
    exon_count[parent]++
    total_len  += (end - start + 1)
    total_exons++
}
 
END {
    multi = 0; single = 0
    for (tx in exon_count) {
        if (exon_count[tx] >= 2) multi++
        else                      single++
    }
    avg_exon = (total_exons > 0) ? total_len / total_exons : 0
    printf "%d %d %.2f %d\n", multi, single, avg_exon, total_exons
}
' "$GFF")


# ---------------------------------------------------------------------------
#protein-level
# ---------------------------------------------------------------------------

AVG_PROTEIN=$(awk '
/^>/ {
    if (seq_len > 0) {
        total += seq_len
        count++
    }
    seq_len = 0
    next
}
{
    # Remove whitespace and add to current sequence length
    gsub(/[[:space:]]/, "")
    seq_len += length($0)
}
END {
    if (seq_len > 0) {
        total += seq_len
        count++
    }
    if (count > 0) printf "%.2f\n", total / count
    else            print "0.00"
}
' "$PROTEIN")

# ---------------------------------------------------------------------------
# Compose JSON
# ---------------------------------------------------------------------------
JSON=$(cat <<EOF
{
    "assembly": {
        "num_contigs": $ASM_NUM_CONTIGS,
        "length": $ASM_TOTAL_LEN,
        "mean_length": $ASM_MEAN_LEN,
        "N50": $ASM_N50,
        "L50": $ASM_L50,
        "N90": $ASM_N90,
        "L90": $ASM_L90,
        "GC_content": $ASM_GC
    },
    "annotation": {
        "genes": $ANN_GENES,
        "common_name": $ANN_COMMON_NAME,
        "mRNA": $ANN_MRNA,
        "tRNA": $ANN_TRNA,
        "ncRNA": $ANN_NCRNA,
        "rRNA": $ANN_RRNA,
        "avg_gene_length": $ANN_AVG_GENE_LEN,
        "transcript-level": {
            "total_exons": $TOTAL_EXONS,
            "multiple_exon_transcript": $MULTI,
            "single_exon_transcript": $SINGLE,
            "avg_exon_length": $AVG_EXON,
            "avg_protein_length": $AVG_PROTEIN,
            "functional": {
                "go_cog_terms": 13624,
                "interproscan": 7653,
                "eggnog": 8658,
                "pfam": 6403,
                "cazyme": 400,
                "merops (peptidase)": 295,
                "busco": 1252,
                "secretion": 719
            },
        }
    }
}
EOF
)

if [[ -n "$OUTPUT" ]]; then
    echo "$JSON" > "$OUTPUT"
    echo "Stats written to $OUTPUT" >&2
else
    echo "$JSON"
fi

#            "CDS_transcripts": $CDS_TRANSCRIPTS,
#            "CDS_five_utr": $CDS_FIVE_UTR,
#            "CDS_three_utr": $CDS_THREE_UTR,
#            "CDS_no_utr": $CDS_NO_UTR,
#            "CDS_five_three_utr": $CDS_FIVE_THREE_UTR,
#            "CDS_complete": $CDS_COMPLETE,
#            "CDS_no-start": $CDS_NO_START,
#            "CDS_no-stop": $CDS_NO_STOP,
#            "CDS_no-start_no-stop": $CDS_NO_START_NO_STOP,
#            "total_cds_exons": $TOTAL_CDS_EXONS,
#            "pct_exon_overlap_protein_evidence": $PCT_EXON_OVERLAP
