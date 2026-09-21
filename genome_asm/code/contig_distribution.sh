#!/bin/bash

#contig_distribution.sh to assess how contig length varies across assembly - ideally concentrated primarily within chromosome-level scaffolds

asm="/home/padgett3/scratch/streptotinia_caulophylli/hic_asm/SC_hifihic_020326.asm.hic.p_ctg.fa"

awk '
  /^>/ {
    if (seqlen) {
      names[++n] = name
      lengths[n] = seqlen
      total += seqlen
    }
    name = substr($0,2)
    seqlen = 0
    next
  }
  {
    seqlen += length($0)
  }
  END {
    names[++n] = name
    lengths[n] = seqlen
    total += seqlen

    for (i = 1; i <= n; i++) {
      print names[i], lengths[i], lengths[i] / total
    }
  }
' $asm > sc_hifiasm_contig_len.tsv

#there should be a column for the sequence name, sequence length, and proportional sequence length

#export .tsv to RStudio and continue in cumulative_contig_distribution_plot.R