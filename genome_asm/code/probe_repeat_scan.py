#!/usr/bin/env python3

#probe_repeat_scan.py script uses a telomere probe used during wet lab visualization/tracking of chromosomes, applied against genome asm sequences

import re
import argparse
import sys
from pathlib import Path


# ── FASTA parser ──────────────────────────────────────────────────────────────

def parse_fasta(fasta_path: str):
    """Yield (contig_name, sequence) tuples. Handles multi-line sequences."""
    contig_name = None
    seq_parts = []

    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if contig_name is not None:
                    yield contig_name, "".join(seq_parts)
                contig_name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

    if contig_name is not None:
        yield contig_name, "".join(seq_parts)


# ── Repeat finder ─────────────────────────────────────────────────────────────

def find_repeats(sequence: str, motif: str, min_repeats: int,
                 case_sensitive: bool, strand: str):
    """
    Find all non-overlapping runs of `motif` repeated >= `min_repeats` times.

    Returns a list of dicts with keys:
        start, end, repeat_count, matched_seq, strand
    """
    if not case_sensitive:
        search_seq   = sequence.upper()
        search_motif = motif.upper()
    else:
        search_seq   = sequence
        search_motif = motif

    pattern = re.compile(
        r"(?:" + re.escape(search_motif) + r"){" + str(min_repeats) + r",}"
    )

    hits = []
    for m in pattern.finditer(search_seq):
        start = m.start()
        end   = m.end()
        hits.append({
            "start":        start,
            "end":          end,
            "repeat_count": len(m.group()) // len(search_motif),
            "matched_seq":  sequence[start:end],   # original case from FASTA
            "strand":       strand,
        })
    return hits


# ── TSV writer ────────────────────────────────────────────────────────────────

HEADER = "\t".join([
    "contig", "start_1based", "end_1based", "length_bp",
    "repeat_count", "strand", "sequence"
])


def write_hit(fh, contig: str, hit: dict):
    fh.write("\t".join([
        contig,
        str(hit["start"] + 1),          # 0-based → 1-based
        str(hit["end"]),                 # exclusive end = inclusive 1-based end
        str(hit["end"] - hit["start"]),
        str(hit["repeat_count"]),
        hit["strand"],
        hit["matched_seq"],
    ]) + "\n")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Scan a genome FASTA for tandem repeats on both strands.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-f",  "--fasta",    required=True,
                        help="Genome FASTA file")
    parser.add_argument("-m",  "--motif",    required=True,
                        help="Forward (+) strand repeat motif (e.g. CCCTAA)")
    parser.add_argument("-rc", "--rc-motif", default=None,
                        help="Reverse-complement (-) strand motif (e.g. TTAGGG). "
                             "If omitted, only the forward motif is searched.")
    parser.add_argument("-r",  "--repeats",  type=int, required=True,
                        help="Minimum consecutive repeats to report")
    parser.add_argument("-o",  "--output",   default="repeat_hits.tsv",
                        help="Output TSV file (default: repeat_hits.tsv)")
    parser.add_argument("--case-sensitive",  action="store_true", default=False,
                        help="Case-sensitive matching (default: off)")

    args = parser.parse_args()

    # ── Validate ─────────────────────────────────────────────────────────────
    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        print(f"[ERROR] FASTA not found: {fasta_path}", file=sys.stderr)
        sys.exit(1)

    motif = args.motif.strip()
    if not motif:
        print("[ERROR] --motif cannot be empty.", file=sys.stderr)
        sys.exit(1)

    rc_motif = args.rc_motif.strip() if args.rc_motif else None

    if args.repeats < 1:
        print("[ERROR] --repeats must be >= 1.", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ── Info block ───────────────────────────────────────────────────────────
    print(f"[INFO] FASTA          : {fasta_path}")
    print(f"[INFO] Forward motif  : {motif}  (strand +)")
    if rc_motif:
        print(f"[INFO] RC motif       : {rc_motif}  (strand -)")
    else:
        print(f"[INFO] RC motif       : not provided — forward strand only")
    print(f"[INFO] Min repeats    : {args.repeats}")
    print(f"[INFO] Case-sensitive : {args.case_sensitive}")
    print(f"[INFO] Output         : {output_path}")
    print()

    total_hits     = 0
    total_contigs  = 0

    with open(output_path, "w") as out_fh:
        out_fh.write(HEADER + "\n")

        for contig_name, sequence in parse_fasta(str(fasta_path)):
            total_contigs += 1
            seq_len = len(sequence)

            # Forward strand search
            fwd_hits = find_repeats(sequence, motif,    args.repeats,
                                    args.case_sensitive, strand="+")

            # Reverse complement search (searches same forward sequence for RC motif)
            rc_hits = []
            if rc_motif:
                rc_hits = find_repeats(sequence, rc_motif, args.repeats,
                                       args.case_sensitive, strand="-")

            # Merge and sort by start position
            all_hits = sorted(fwd_hits + rc_hits, key=lambda h: h["start"])
            contig_hits = len(all_hits)
            total_hits += contig_hits

            print(
                f"  Contig: {contig_name:<35}  length: {seq_len:>12,} bp  "
                f"hits: {contig_hits:>6,}  "
                f"(+:{len(fwd_hits)}  -:{len(rc_hits)})"
            )

            for hit in all_hits:
                write_hit(out_fh, contig_name, hit)

    print()
    print(f"[DONE] Scanned {total_contigs:,} contig(s), found {total_hits:,} total hit(s).")
    print(f"[DONE] Results written to: {output_path}")


if __name__ == "__main__":
    main()