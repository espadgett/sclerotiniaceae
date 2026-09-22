#!/bin/bash

GFF="${1:-}"
OUT="${2:-/dev/stdout}"

[[ -z "$GFF" ]] && { echo "Usage: $0 <input.gff3> [output.tsv]"; exit 1; }
[[ ! -f "$GFF" ]] && { echo "ERROR: file not found: $GFF"; exit 1; }

awk '
BEGIN { OFS="\t" }

# Skip comment and blank lines
/^#/ { next }
/^\s*$/ { next }

# Only process mRNA features (col 3)
$3 != "mRNA" { next }

{
    col9 = $9
    gene_id = ""
    go_str  = ""

    # --- Extract ID= value ---
    if (match(col9, /ID=([^;]+)/, arr)) {
        gene_id = arr[1]
    }

    # --- Extract Ontology_term= value (everything up to next ;) ---
    if (match(col9, /Ontology_term=([^;]+)/, arr)) {
        go_str = arr[1]
    }

    # Skip lines with no GO terms
    if (gene_id == "" || go_str == "") next

    # --- Split GO terms on comma and print one per line ---
    n = split(go_str, terms, ",")
    for (i = 1; i <= n; i++) {
        # Only print valid GO terms (ignore any non-GO entries)
        if (terms[i] ~ /^GO:[0-9]+$/) {
            print gene_id, terms[i]
        }
    }
}
' "$GFF" > "$OUT"

# Report to stderr so it does not pollute stdout/tsv
if [[ "$OUT" != "/dev/stdout" ]]; then
    LINES=$(wc -l < "$OUT")
    GENES=$(cut -f1 "$OUT" | sort -u | wc -l)
    TERMS=$(cut -f2 "$OUT" | sort -u | wc -l)
    echo "Done: ${LINES} lines | ${GENES} unique genes | ${TERMS} unique GO terms -> ${OUT}" >&2
fi