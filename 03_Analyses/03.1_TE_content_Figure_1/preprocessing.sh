# Generating summaries

## Build summary per TE class
buildSummary_classes_csv.pl -genome_size genome_size rpm.out

## Build summary per genomic sequence
buildSummary_perseq_csv.pl -genome_size genome_size rpm.out

## Generate TE divergence summary
calcDivergenceFromAlign.pl -s output_file align_file


### Exact command list
## A io ##
cd /tmp/global2/ssushko/Vanessa_proj/repeat_annotation/with_NEWtameameagenome/a_io
buildSummary_classes_csv.pl -genome_size 384156914 A_io_GCF_905147045.1_genomic_refseq.fna.filtered.out > a_io_mch_classes.sum
buildSummary_families_csv.pl -genome_size 384156914 A_io_GCF_905147045.1_genomic_refseq.fna.filtered.out > a_io_mch_families.sum
buildSummary_perseq_csv.pl -genome_size 384156914 A_io_GCF_905147045.1_genomic_refseq.fna.filtered.out > a_io_mch_perseq.sum
# source activate /ebio/abt6_projects9/abt6_software/conda/repeat_tools/dfam-tetools_latest
calcDivergenceFromAlign.pl -s a_io_div.sum A_io_GCF_905147045.1_genomic_refseq.fna.align

## V card ##
cd /tmp/global2/ssushko/Vanessa_proj/repeat_annotation/with_NEWtameameagenome/v_card
buildSummary_classes_csv.pl -genome_size 424813635 V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.out > v_card_mch_classes.sum
buildSummary_families_csv.pl -genome_size 424813635 V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.out > v_card_mch_families.sum
buildSummary_perseq_csv.pl -genome_size 424813635 V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.out > v_card_mch_perseq.sum
calcDivergenceFromAlign.pl -s v_card_div.sum V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.align

## V atal ##
cd /tmp/global2/ssushko/Vanessa_proj/repeat_annotation/with_NEWtameameagenome/v_atal
buildSummary_classes_csv.pl -genome_size 370420815 V_atal_GCF_905147765.1_genomic_refseq.fna.filtered.out > v_atal_mch_classes.sum
buildSummary_families_csv.pl -genome_size 370420815 V_atal_GCF_905147765.1_genomic_refseq.fna.filtered.out > v_atal_mch_families.sum
buildSummary_perseq_csv.pl -genome_size 370420815 V_atal_GCF_905147765.1_genomic_refseq.fna.filtered.out > v_atal_mch_perseq.sum
calcDivergenceFromAlign.pl -s v_atal_div.sum V_atal_GCF_905147765.1_genomic_refseq.fna.align

## V tam ##
cd /tmp/global2/ssushko/Vanessa_proj/repeat_annotation/with_NEWtameameagenome/v_tam
# buildSummary_classes_csv.pl -genome_size 357124929 V_tam_GCF_002938995.1_genomic_refseq.fna.out > v_tam_mch_classes.sum

buildSummary_classes_csv.pl -genome_size 363368026 GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna.filtered.out > v_tam_mch_classes.sum
buildSummary_families_csv.pl -genome_size 363368026 GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna.filtered.out > v_tam_mch_families.sum
buildSummary_perseq_csv.pl -genome_size 363368026 GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna.filtered.out > v_tam_mch_perseq.sum
calcDivergenceFromAlign.pl -s v_tam_div.sum GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna.align