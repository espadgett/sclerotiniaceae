#!/bin/bash

set -euo pipefail
module load cd-hit # cd-hit/4.8.1

FASTA1="${1:?Usage: $0 <fasta1> <fasta2> [identity] [threads] [outdir]}"
FASTA2="${2:?Usage: $0 <fasta1> <fasta2> [identity] [threads] [outdir]}"
IDENTITY="${3:-0.90}"
THREADS="${4:-4}"
OUTDIR="${5:-cdhit_results}"

LABEL1=$(basename "$FASTA1" | sed 's/\.[^.]*$//')
LABEL2=$(basename "$FASTA2" | sed 's/\.[^.]*$//')

# Word length for cd-hit
ident_int=$(echo "$IDENTITY * 100" | bc | cut -d. -f1)
if   [ "$ident_int" -ge 70 ]; then WORD_LEN=5
elif [ "$ident_int" -ge 60 ]; then WORD_LEN=4
elif [ "$ident_int" -ge 50 ]; then WORD_LEN=3
else                                WORD_LEN=2
fi

# ── Checks ────────────────────────────────────────────────────────────────────
for f in "$FASTA1" "$FASTA2"; do
  [ -f "$f" ] || { echo "ERROR: File not found: $f"; exit 1; }
done
command -v cd-hit >/dev/null || { echo "ERROR: cd-hit not found in PATH"; exit 1; }

mkdir -p "$OUTDIR"

SEP="============================================================"
echo ""
echo "$SEP"
echo "  CD-HIT Proteome Comparison  |  identity >= ${ident_int}%"
echo "$SEP"

# ── 1. Extract & sort IDs ─────────────────────────────────────────────────────
echo "[1/4] Reading FASTA files..."

grep "^>" "$FASTA1" | sed 's/^>//; s/ .*//' | sort > "$OUTDIR/ids1.txt"
grep "^>" "$FASTA2" | sed 's/^>//; s/ .*//' | sort > "$OUTDIR/ids2.txt"

COUNT1=$(wc -l < "$OUTDIR/ids1.txt")
COUNT2=$(wc -l < "$OUTDIR/ids2.txt")
echo "   $LABEL1: $COUNT1 sequences"
echo "   $LABEL2: $COUNT2 sequences"

OVERLAP=$(comm -12 "$OUTDIR/ids1.txt" "$OUTDIR/ids2.txt" | wc -l)
if [ "$OVERLAP" -gt 0 ]; then
  EXAMPLE=$(comm -12 "$OUTDIR/ids1.txt" "$OUTDIR/ids2.txt" | head -1)
  echo "ERROR: $OVERLAP overlapping IDs found (e.g. $EXAMPLE)."
  exit 1
fi

# ── 2. Combine ────────────────────────────────────────────────────────────────
COMBINED="$OUTDIR/combined.fasta"
echo "[2/4] Combining FASTA files -> $COMBINED ..."
cat "$FASTA1" "$FASTA2" > "$COMBINED"

# ── 3. Run CD-HIT ─────────────────────────────────────────────────────────────
CLUSTERED="$OUTDIR/clustered"
echo "[3/4] Running CD-HIT (identity=$IDENTITY, word_len=$WORD_LEN)..."
cd-hit \
  -i "$COMBINED" \
  -o "$CLUSTERED" \
  -c "$IDENTITY" \
  -n "$WORD_LEN" \
  -T "$THREADS" \
  -M 800 \
  -d 0 \
  -g 1 \
  > "$OUTDIR/cdhit.log" 2>&1
echo "   Done. Log: $OUTDIR/cdhit.log"

# ── 4. Parse .clstr and classify ──────────────────────────────────────────────
echo "[4/4] Parsing clusters and writing results..."

CLSTR="${CLUSTERED}.clstr"
TSV="$OUTDIR/shared_clusters.tsv"
UNIQUE1="$OUTDIR/unique_to_${LABEL1}.txt"
UNIQUE2="$OUTDIR/unique_to_${LABEL2}.txt"

printf "cluster_id\trepresentative\t%s_members\t%s_members\n" "$LABEL1" "$LABEL2" > "$TSV"
> "$UNIQUE1"
> "$UNIQUE2"

# mawk-compatible: no gensub, no match with 3-arg array
gawk -v ids1="$OUTDIR/ids1.txt" \
     -v ids2="$OUTDIR/ids2.txt" \
     -v tsv="$TSV" \
     -v u1="$UNIQUE1" \
     -v u2="$UNIQUE2" \
     -v label1="$LABEL1" \
     -v label2="$LABEL2" \
'BEGIN {
    while ((getline line < ids1) > 0) f1[line] = 1
    while ((getline line < ids2) > 0) f2[line] = 1
    sc=0; ss1=0; ss2=0; o1c=0; o1s=0; o2c=0; o2s=0
}
/^>Cluster/ {
    if (clid != "") flush()
    clid = $2; n = 0; rep = ""
    delete members
}
/^[0-9]/ {
    # extract ID between > and ...
    seq = $0; sub(/.*>/, "", seq); sub(/\.\.\..*/, "", seq)
    members[n++] = seq
    if ($0 ~ /\*$/) rep = seq
}
END { if (clid != "") flush() }
function flush(    i, has1, has2, m1, m2, seq) {
    has1=0; has2=0; m1=""; m2=""
    for (i=0; i<n; i++) {
        seq = members[i]
        if (seq in f1) { has1++; m1 = (m1=="" ? seq : m1 ";" seq) }
        if (seq in f2) { has2++; m2 = (m2=="" ? seq : m2 ";" seq) }
    }
    if (has1 && has2) {
        print clid "\t" rep "\t" m1 "\t" m2 >> tsv
        sc++; ss1+=has1; ss2+=has2
    } else if (has1) {
        for (i=0; i<n; i++) if (members[i] in f1) print members[i] >> u1
        o1c++; o1s+=has1
    } else {
        for (i=0; i<n; i++) if (members[i] in f2) print members[i] >> u2
        o2c++; o2s+=has2
    }
}
END {
    print sc "\t" ss1 "\t" ss2 "\t" o1c "\t" o1s "\t" o2c "\t" o2s
}' "$CLSTR" | tail -1 | \
{
  IFS=$'\t' read -r sc ss1 ss2 o1c o1s o2c o2s
  TOTAL=$(grep -c "^>Cluster" "$CLSTR")
  PCT1=$(echo "scale=1; $ss1 * 100 / $COUNT1" | bc)
  PCT2=$(echo "scale=1; $ss2 * 100 / $COUNT2" | bc)

  echo ""
  echo "$SEP"
  echo "  RESULTS SUMMARY"
  echo "$SEP"
  printf "  Identity threshold  : %s%%\n"   "$ident_int"
  printf "  Total clusters      : %s\n"     "$TOTAL"
  echo ""
  echo "  Shared clusters (members from BOTH files):"
  printf "    Clusters          : %s\n"     "$sc"
  printf "    From %-20s: %6s  (%s%%)\n"   "$LABEL1" "$ss1" "$PCT1"
  printf "    From %-20s: %6s  (%s%%)\n"   "$LABEL2" "$ss2" "$PCT2"
  echo ""
  printf "  Unique to %s:\n"               "$LABEL1"
  printf "    Clusters          : %s\n"     "$o1c"
  printf "    Sequences         : %s\n"     "$o1s"
  echo ""
  printf "  Unique to %s:\n"               "$LABEL2"
  printf "    Clusters          : %s\n"     "$o2c"
  printf "    Sequences         : %s\n"     "$o2s"
  echo ""
  echo "  Output files:"
  printf "    %s\n" "$TSV"
  printf "    %s\n" "$UNIQUE1"
  printf "    %s\n" "$UNIQUE2"
  printf "    %s\n" "$CLSTR"
  echo "$SEP"
  echo ""
}