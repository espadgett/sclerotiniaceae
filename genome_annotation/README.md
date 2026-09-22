
# Genome Annotation

## Analyses of genome annotation (gene and protein features), following complete genome assembly in *genome_asm* folder.  

### Annotations were conducted using Funannotate, and we advise visiting their Github at https://github.com/nextgenusfs/funannotate. Analyses use protein level annotations as input for effector identification, and ortholog genes as input for synteny analyses. 

### Code contained in the *genome_annotation* folder:

####**annotation_data**: supplementary data for the annotation analyses
- effector_overlap_data: details about effector protein similarities between species
- oascan_data: closest protein hits for OAH and ODC2 genes
- funannotate: Funannotate results **and effector protein annotations**, derived from Funannotate protein fasta

####**phylogeny**: scripts to use ortholog gene annotations to create phylogenetic trees and synteny plots across Sclerotinia sclerotiorum, S. minor, Streptotinia caulophylli, and Botrytis cinerea
- synteny analysis -> *genespace_initialize.sh/genespace_run.R/riparian_plot_visualization.R*
- phylogenetic tree -> *riparian_visualization.R*

**figure subfolder**: contains figure file for Figure 1

FIGURE Synteny and phylogenetic tree of Sclerotinia sclerotiorum, S. minor, Streptotinia caulophylli, and Botrytis cinerea

####**secretome**: Detect effector elements from Funannotate proteins annotations and scan for major virulence genes. The effector annotation pipeline largely filters based on structural prediction (SignalP, Phobius (optional), NetGPI) and known effector sequence motifs (EffectorP, PHI-DB). Functional annotation of effector candidates are provided in **annotation_data/funannotate**.
- effector annotation from Funannotate outputs -> *effector_annotation.sh*
- scanning for major virulence genes OAH and OCD2 -> *oa_scan.sh*
- effector similarity across species -> *overlap_parse.sh*

Virulence protein scan covered the following four proteins:
- OAH1 (GenBank: GU932672.1, from https://pmc.ncbi.nlm.nih.gov/articles/PMC2924111/)
- OAHA 1 (AY590264.2, from https://www.jbc.org/article/S0021-9258(19)33612-9/fulltext#references33612-9/fulltext#references)) and 2 (Bcin12p01020.1 = ATZ55511.1)
- OAD2 aka sscle_09g069850 (SsODC2)
```
>SsOAH1
MAPIMDALPIVEDKPTAAVTDFSVSSPQDGAVVAPPTTVYQTGATKLKNMLRDSNELIVCPGVYDGISTRVALQVGFPAL
YMTGAGTTASRLGMADLGIAHLSDMKDHAEMIANLDPFGPPLIADMDTGYGGPLIVDKAVKAYIRAGVAGFHIEDQIQNK
RCGHLAGKKVVPEEEYYMRIRAAKGAKDAMKSDIVLIARTDALQQLGYDECVKRLKVARELGADVGLLEGYTSKEMAAKT
VKEFAPWPILLNMVENGATPIITTKEAQEMGFRIMIFSFAALAPAMLAIQETFVRLKNEGVVGTPKNVTPRALFEVCGLQ
ESIVIDTAAGGGAFADGV
>BcOAHA_1
MASTIAVQEIRLDTTSHSLDRSETPSPIGPLYATGSSSPLKAST
TWLQLNSAAPIYVKVQNGSTAAEDEPFSGAKKLRHLLENTDELIVCPGVYDGLSARTA
MELGFKSLYMTGAGTTASRLGQPDLAIAQLHDMRDNADMIANLDPFGPPLIADMDTGY
GGPIMVARTVEHYIRSGVAGAHLEDQILTKRCGHLSGKKVVSRDEYLVRIRAAVATKR
RLRSDFVLIARTDALQSLGYEECIERLRAARDEGADVGLLEGFRSKEQAAAAVAALAP
WPLLLNSVENGHSPLITVEEAKAMGFRIMIFSFATLAPAYAAIRETLVRLRDHGVVGT
PDGITPVRLFEVCGLQDAMEVDNGAGGKAFSEGV
>BcOAHA_2
MPAYSDKVMLTFNNTVKQTTGFESGATKLKNMLRDSNELIVCPG
VYDGISARVALQVGFPALYMTGAGTTASRLGMADLGIAQLSDMKDHAEMIANLDPYGP
PLIADMDTGYGGPLIIDKAVKSYIRAGVAGFHIEDQIQNKRCGHLQGKKVVPAEEYYM
RIRAAKAAKEAMNSDIVLIARTDALQQLGYDECIKRLKVAREMGADVGLLEGYTSKEM
AAKTVKEFAPWPILLNMVENGSTPIITTKEAQEMGFRIMIFSFAALAPAYLAIQETFL
RLKRDGVVGTPKNLTPKALFDVCGLKDSIVVDTTAGGGAFADGV
>SsODC_2
MHSKTFLLTSLAAGGLLPAVLAAPAPRPQMNYAEQSAAIAGVSAPLPTATAASGSLYGHKSLLGESGKSSSLSKADSAIM
PDTPLVIGQEADSDLGLFLDFNSVDNPQPIRGSTGETDPGPQTYEYQRLNPDLLAPPGTDSGDMPQAQWPMGMSHNRLGV
GKKSGWARQQNTGQMPIATEMAGVDMRLAPGAYRELHWHTANEWSLVLKGSVRISAVNEKGESFIDDISAGDVWFFPAGV
PHSIQALDEGCEFLLIFDDGDFSEDGTFLVSEMFLRNPKSVLSKNFQTPVSSFDKLPSDQLYIFNGTPAPANISTQNQTS
DAGPLPQKDSYTYHFSQQEPYTVPGGSVKIIDTTTFPIAANFAAALVTIEPGAMRELHWHLTSDEWNYFLSGSARITIYS
APSSSRTFDYTAGDVGYIPASNSHYIENTGTEDVVFLEVLQQPVFSDISVAQWLALTPRQVVKDTLGLPDDVIDSLPHEK
TLIKVGNRNMTALAADEKGSGAYE
```

####**supp_materials**: supplementary annotation analyses
- track gene density across genomes -> *gene_density.R*

FIGURE Gene density in Streptotinia caulophylli

**effector_annotation**: additional annotations on the highest-confidence effectors (containing PHI-DB annotations)
- BLASTKoala (web server) KO annotations/KEGG annotations -> *kegg.R*
- assign GO terms to effector candidates -> *filter_gff3.sh*/*filter_gff3to_go.sh*
- test GO enrichment across effector candidates-> *go_enrichment.R*
**funannotate_filter**: To save compute and time resources, existing Funannotate results were filtered as follows to keep only high-confidence contigs. Summary statistics were also produced.
- filter Funannotate results, and summarize some annotation statistics -> *funannotate_filter.sh/filtered_stats.sh*
#### Motive for filtering
*Funannotate was initially run on an older version of the hifiasm genome assembly, still containing duplicates before the **genome_asm/code/rename_juicer_contigs.sh** purging (called in **genome_asm/code/genome_asm.sh**)*

#### Filtering targets:
```(Funannotate_name/sequence_length/new_name)

St_cau_sc_1 3586059 #chr1
St_cau_sc_2 3139222 #chr2
St_cau_sc_3 3014892 #chr5
St_cau_sc_4 2992952 #chr3
St_cau_sc_5 2619870 #chr7
St_cau_sc_6 2551806 #chr6
St_cau_sc_7 2521891 #chr9
St_cau_sc_8 2395583 #chr10
St_cau_sc_9 2392683 #chr8
St_cau_sc_10 2315771 #chr12
St_cau_sc_11 2265534 #chr11
St_cau_sc_12 2213651 #chr4
St_cau_sc_13 1910328 #chr15
St_cau_sc_14  1896953 #chr13
St_cau_sc_15  1837025 #chr16
St_cau_sc_16 1811447 #chr14
St_cau_sc_17 67744 #scaffold1
St_cau_sc_19 57938 #scaffold2
St_cau_sc_29 44251 #scaffold3
St_cau_sc_30 43778 #scaffold4
St_cau_sc_33 42254 #scaffold5
St_cau_sc_51 36738 #scaffold6

```



####Referenced code
** genespace_initialize.sh**: creates environment/input structure for GENESPACE ortholog-based synteny calculation using the outputs of Funannotate (see README.md). For annotation data deposition, see associated publication.
** genespace_run.R**: run GENESPACE, very mem/time intensive. Recommended to run in an interactive R environment on HPC.
** riparian_plot_visualization.R**: still in an interactive R environment, load the GENESPACE output to create a riparian plot. Funannotate reorders contigs, requiring script correction.
** effector_annotation.sh**: annotate effector candidates from structural and sequence motifs, filtering Funannotate protein.fasta
** oa_scan.sh**: scanning protein.fasta of Sclerotinia sclerotiorum, S. minor, Streptotinia caulophylli, and Botrytis cinerea for virulence genes. 
** overlap_parse.sh**: assess the uniqueness of effector sets between species, two at a time, by comparing two protein FASTA files using CD-HIT at 90% identity.
	*Usage*
	```overlap_parse.sh <fasta1> <fasta2> [identity] [threads] [outdir]```
** filter_gff3.sh**: extract effectors from whole genome GFF3. Insert full paths or alternatively:
	*Usage*
	```filter_gff3.sh -g input.gff3 -i ids.txt -o filtered.gff3```
** filter_gff3to_go.sh**: extracts the GO terms from the Funannotate-ouput GFF3 file and output a two-column TSV (#   gene_id <TAB> GO:XXXXXXX).
	*Usage*
	```filter_gff3_to_go.sh phi_pos.gff3 > phi_pos_go.tsv``
** gene_density.R**: taking the funannotate output GFF3, visualize the density of the genes (genes/Mb) across chromosomes. Bin size in image is set as 1 Mb, but is customizable.
** go_enrichment.R**: clusterProfiler GO enrichment plot

FIGURE GO term enrichment for S. minor and S. caulophylli highest-confidence effector candidates

** kegg.R**: KO term enrichment
** funannotate_filter.sh**: see above
** filtered_stats.sh**: summary statistics for filtered annotations

END of README

