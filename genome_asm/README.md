
# Genome assembly for S. minor and St. caulophylli

## Genome assembly and manual curation using HiFi and HiC sequences of Streptotinia caulophylli

### Raw sequencing stats used as input
- S. minor
    - 6 Gb (150X) of HiFi sequences
    - 5 Gb (125X) of HiC sequences
- St. caulophylli
    - 6 Gb (150X) of HiFi sequences
    - 2.4 Gb (60X) of HiC sequences

### Code contained in the *genome_asm* folder:

**genome_asm.sh**: primary script for Streptotinia caulophylli genome assembly with hifiasm and some summary QC which helped manually curate assembly
- hifiasm assembly
- optional deduplication (sanity check)
- BUSCO
- QUAST
- repeat annotation -> *renaming_asm.sh*
- telomere annotation -> *probe_repeat_scan.py*
- YaHS + Juicer HiC visualization -> *rename_juicer_contigs.sh*

**supplementary_asm_analyses.sh**
- dot plot -> *dotplot.sh*
- contain distribution plot -> *cumulative_contig_distribution_plot.R*

Embedded codes
**renaming_asm**: reorders the hifiasm-output assembly scaffolds in harmony with the B. cinerea contigs, determined by ortholog synteny
**probe_repeat_scan.py*: script uses a telomere probe used during wet lab visualization/tracking of chromosomes, applied against genome asm sequences
**rename_juicer_contigs.sh*: after Juicer (Desktop) visualization, remove low-confidence contigs
**dotplot.sh*: visualize genome asm (contig-alignment)
**cumulative_contig_distribution_plot.R*: visualizes the distribution of contigs by length, and add the cumulative contribution to genome asm length

END of README

