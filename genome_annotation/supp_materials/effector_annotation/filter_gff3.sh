#!/bin/bash

set -euo pipefail
# --- Parse arguments ---------------------------------------------------------
#GFF3 from funannotate with GO terms annotated
GFF=[1] #list of PHI-positive (effector) hits
IDS=[2] #extract IDs of effector hits
OUT=[3]

while getopts "g:i:o:" opt; do
    case $opt in
        g) GFF="$OPTARG" ;;
        i) IDS="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "$GFF" || -z "$IDS" || -z "$OUT" ]] && usage
[[ ! -f "$GFF" ]] && { echo "ERROR: GFF3 file not found: $GFF"; exit 1; }
[[ ! -f "$IDS" ]] && { echo "ERROR: IDs file not found: $IDS"; exit 1; }

# --- Build ID sets in a temp directory ---------------------------------------
TMPDIR_LOCAL=$(mktemp -d)
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

MRNA_IDS="$TMPDIR_LOCAL/mrna_ids.txt"
GENE_IDS="$TMPDIR_LOCAL/gene_ids.txt"

# Strip blank lines and comment lines from the IDs file
grep -v '^\s*$\|^#' "$IDS" | sort -u > "$MRNA_IDS"
MRNA_COUNT=$(wc -l < "$MRNA_IDS")
echo "[1/3] Loaded ${MRNA_COUNT} mRNA IDs from '${IDS}'" >&2

# --- Pass 1: resolve parent gene IDs -----------------------------------------
# Strategy A: funannotate convention  GENE_ID-T\d+  ->  GENE_ID
sed 's/-T[0-9]*$//' "$MRNA_IDS" | sort -u > "$GENE_IDS"

# Strategy B: extract Parent= value from mRNA lines whose ID= is in our list
# awk reads the GFF3, skips comments, matches mRNA lines, checks ID= vs list
awk -v mrna_file="$MRNA_IDS" '
BEGIN {
    while ((getline id < mrna_file) > 0) wanted[id] = 1
}
/^#/ { next }
{
    if ($3 != "mRNA") next
    # parse ID= and Parent= from column 9
    n = split($9, attrs, ";")
    id = ""; parent = ""
    for (i = 1; i <= n; i++) {
        if (attrs[i] ~ /^ID=/)     id     = substr(attrs[i], 4)
        if (attrs[i] ~ /^Parent=/) parent = substr(attrs[i], 8)
    }
    if (id in wanted && parent != "") print parent
}
' "$GFF" | sort -u >> "$GENE_IDS"

sort -u "$GENE_IDS" -o "$GENE_IDS"
GENE_COUNT=$(wc -l < "$GENE_IDS")
echo "[2/3] Resolved ${GENE_COUNT} parent gene IDs" >&2

# --- Pass 2: filter GFF3 -----------------------------------------------------
awk -v mrna_file="$MRNA_IDS" -v gene_file="$GENE_IDS" '
BEGIN {
    while ((getline id < mrna_file) > 0) mrna[id] = 1
    while ((getline id < gene_file) > 0) gene[id] = 1
}

# Always keep header / comment lines
/^#/ { print; next }

# Skip blank lines
/^\s*$/ { next }

{
    # parse ID= and Parent= from column 9
    n = split($9, attrs, ";")
    feat_id = ""; parent = ""
    for (i = 1; i <= n; i++) {
        if (attrs[i] ~ /^ID=/)     feat_id = substr(attrs[i], 4)
        if (attrs[i] ~ /^Parent=/) parent  = substr(attrs[i], 8)
    }

    keep = 0

    if ($3 == "gene") {
        if (feat_id in gene) keep = 1

    } else if ($3 == "mRNA") {
        if (feat_id in mrna) keep = 1

    } else {
        # exon / CDS / UTR / etc.
        # Parent= may be comma-separated (multi-parent features)
        np = split(parent, parents, ",")
        for (j = 1; j <= np; j++) {
            if (parents[j] in mrna) { keep = 1; break }
        }
    }

    if (keep) print
}
' "$GFF" > "$OUT"

WRITTEN=$(grep -v '^#' "$OUT" | wc -l)
echo "[3/3] Wrote ${WRITTEN} feature lines to '${OUT}'" >&2      