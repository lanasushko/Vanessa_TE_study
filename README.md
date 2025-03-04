# Code repository for annotation and analysis of TEs in *Vanessa* butterflies

This repository contains code for TE annotation of 4 buttefly species and the subsequent analyses that have been performed.

### 01_TE_discovery_and_library_construction

De novo TE libraries for the four species were constructed independently using [RepeatModeler2](https://github.com/Dfam-consortium/RepeatModeler), version 2.0.4 for *V. cardui* and *V. atalanta*; and version 2.0.1 for *V. tameamea* and *A. io* with LTR structural module as in `repeatmodeler.sh`. MITE families were detected and added separately using [MITE Tracker](https://github.com/INTABiotechMJ/MITE-Tracker) pipeline [`mitetracker.sh`]. [MCHelper](https://github.com/GonzalezLab/MCHelper) was run in automatic mode on the four TE libraries to extend fragmented consensus sequences and remove false positives based on structural criteria [`mchelper.sh`].

The resulting libraries were combined with the RepBase Arthropoda library, the Danaus plexippus TE library from Baril & Hayward, 2022 and the Leptidea sinapis library generated by Backström lab in 2022. Satellite, tandem repeat non-coding RNA and low complexity repeat consensus sequences were removed from the standard libraries prior to that. A clustering step was produced with cd-hit-est, version 4.8.1, [`cdhit.sh`] to reduce sequence redundancy. The criteria used for collapsing two or more sequences to a single consensus were 95% identity over at least 80% of their total length. 

The collapsed (non-overlapping) versions of the annotation have been generated using `collapse_RM_annotation.py` tool from [collapse_RepeatMasker_annotaion](https://github.com/shohei-kojima/collapse_RepeatMasker_annotation).

### 02_TE_annotation

Repetitive sequences were annotated in the genomes with RepeatMasker, version 4.1.5 [`repeatmasker.sh`]. The annotation was filtered so that only the RepeatMasker hits with more than 300 of score, less than 40% divergency and more than 80 bp in length were retained in the annotation [`filtering_annotation.sh`].

### 03_Analyses
#### 03.1_TE_content_Figure_1
The TE composition of each genome was obtained using modified versions of `buildSummary.pl` that produce a computer readable format output. The Kimura-2-parameter (K2P) distance for each insertion was calculated using the utility calcDivergenceFromAlign.pl from RepeatMasker toolbox. All the code regarding data preprocessing for TE content analysis is located in `preprocessing.sh`. The code to generate Figure 1A is located in `TE_content_plot.R`, Figure 1B in `chromosome_masking.R` and Figure 1C in `divergence_plots.R`.

Correlation between TE content and chromosome size was tested using linear models (`lm()`) function in R. Correlations between TE content - CDS content and TE content - GC content in 10Kb genomic windows was tested in the same way.

#### 03.2_Genomic_synteny_Figure_2
To investigate genomic rearrangements between *V. cardui* (as reference) and *V. atalanta*, *V. tamemea* and *A. io* genomes (as queries) a whole-genome alignment was performed using [minimap2](https://github.com/lh3/minimap2). We used the -x asm20 parameter recommended for alignments with up to 20% sequence divergence [`mapforsynteny.sh`]. The TE annotations for each species were used to generate window-based representation of TE density per 10KB-long genomic window generated with `bedtools makewindows` [`get_TEbp_per_window.sh`]. 

The synteny was visualized using [gggenomes](https://github.com/thackl/gggenomes/) [`gggenomes_synteny_vanessa.R`].

#### 03.3_Mitochondrial_insertion_Figure_3




#### 03.4_GO_enrichment_Figure_4
Enrichment analysis for the presence of MITEs inside genes and in gene flanking regions was performed using `bedtools` package to retrieve genes and flanks with MITE insertions. Flanking regions for protein-coding genes were extracted using `bedtools flank -b 1000`. Gene annotations were downloaded from respective NCBI genome assembly portals. Proteins from were annotated by [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) to retrieve GO terms associated with each gene. The code for protein ortholog mapping and post-processing to obtain comparable protein seeds and GO terms are located in `eggogmapper.sh` , `get_seeds_for_enrichment.sh` , `produce_background_seeds_list.sh`. Afterwards, [clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler) R package was used to perform enrichment analysis using `enricher()` function [`enrichment.R`].

#### 03.5_Gene_expansions_Figure_5
We used the gene annotation produced by [Shipilina et al, 2022](https://doi.org/10.1016/j.ygeno.2022.110481) for this analysis. The genes corresponding to expanded families identified in this study and their 1000bp-long flanking regions were tested for enrichment in transposable element insertions. The number of TE insertions per gene was extracted using `bedtools intersect` both for the whole set of predicted genes in *Vanessa cardui* genome and specifically the set of genes pertaining to expanded families [`get_number_ins_per_gene.sh`]. The gene flanking regions were extracted using bedtools flank followed by the same procedure as with gene bodies. Afterwards, a T-test with 95% significance threshold was conducted to check for enrichment in TE insertion among the expanded gene set [`expanded_genes.R`].
