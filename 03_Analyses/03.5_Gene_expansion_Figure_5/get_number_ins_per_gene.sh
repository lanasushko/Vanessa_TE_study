# Calculate number of insertions of gene in V. cardui (expanded vs. genomic average)

### GENES ###


# calculate n insertions per expanded gene
cd expanded_genes

cat <(bedtools intersect -a ../gained_genes_filtered_corr.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wa | sort | uniq -c | awk '{print $1,$5,$4-$3}') <(bedtools intersect -a ../gained_genes_filtered_corr.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wao | awk '$13==0{print $0}' | awk '{print "0 "$4,$3-$2}') > expanded_average_per_gene.txt

## Genome-wide Background Rate
# calculate n insertions per gene in the genome
cd expanded_genes

cat <(bedtools intersect -a ../dashas_genes.gff -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wa | sort | uniq -c | awk '{print $1,$10,$6-$5}') <(bedtools intersect -a ../dashas_genes.gff -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wao | awk '$16==0{print $0}' | awk '{print "0 "$9,$5-$4}') > genomic_average_per_gene.txt

## Matched control genes 
# - GC content
# - chromosomal location




## Randomization/Resampling Approach
cd /tmp/global2/ssushko/Vanessa_proj/6_gene_expansion/expanded_genes_cardui
# get random gene sample
#   run get_seed_random script from utils before
shuf -n 77 --random-source=<(get_seeded_random 42) /tmp/global2/ssushko/Vanessa_proj/6_gene_expansion/expanded_genes_cardui/dashas_genes.gff > expanded_genes/random_sample/genes_sample_77_rs42.gff
shuf -n 77 --random-source=<(get_seeded_random 102) /tmp/global2/ssushko/Vanessa_proj/6_gene_expansion/expanded_genes_cardui/dashas_genes.gff > expanded_genes/random_sample/genes_sample_77_rs102.gff
shuf -n 77 --random-source=<(get_seeded_random 12) /tmp/global2/ssushko/Vanessa_proj/6_gene_expansion/expanded_genes_cardui/dashas_genes.gff > expanded_genes/random_sample/genes_sample_77_rs12.gff
shuf -n 77 --random-source=<(get_seeded_random 53) /tmp/global2/ssushko/Vanessa_proj/6_gene_expansion/expanded_genes_cardui/dashas_genes.gff > expanded_genes/random_sample/genes_sample_77_rs53.gff
shuf -n 77 --random-source=<(get_seeded_random 98) /tmp/global2/ssushko/Vanessa_proj/6_gene_expansion/expanded_genes_cardui/dashas_genes.gff > expanded_genes/random_sample/genes_sample_77_rs98.gff

list='expanded_genes/random_sample/genes_sample_77_rs42.gff
expanded_genes/random_sample/genes_sample_77_rs102.gff
expanded_genes/random_sample/genes_sample_77_rs53.gff
expanded_genes/random_sample/genes_sample_77_rs98.gff
expanded_genes/random_sample/genes_sample_77_rs12.gff'

# compute number of insertions
for file in $list;
do
echo 'insertions gene length' > $file.insertions.table
cat <(bedtools intersect -a $file -b V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wa | sort | uniq -c | awk '{print $1,$10,$6-$5}') <(bedtools intersect -a $file -b V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wao | awk '$16==0{print $0}' | awk '{print "0 "$9,$5-$4}') >> $file.insertions.table
done

### GENE FLANKS ###
# expanded
bedtools flank -i gained_genes_filtered_corr.bed -g genomes/ilVanCard2.1/V_card_chrnumbers.genome -b 1000 > flanks/expanded_gene_1000bpflanks.bed
cd flanks
bedtools intersect -a expanded_gene_1000bpflanks.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wa | cut -f4 | sort | uniq -c > expanded_ins_pergeneflanks.txt

cat <(bedtools intersect -a expanded_gene_1000bpflanks.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wa | sort | uniq -c | awk '{print $1,$5}') <(bedtools intersect -a expanded_gene_1000bpflanks.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wao | awk '$13==0{print $0}' | awk '{print "0 "$4}') > expanded_ins_pergeneflanks.txt

## Genome-wide background
# genomic
bedtools flank -i dashas_genes.gff -g genomes/ilVanCard2.1/V_card_chrnumbers.genome -b 1000 > flanks/genomic_gene_1000bpflanks.bed
cd flanks

cat <(bedtools intersect -a genomic_gene_1000bpflanks.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wa | sort | uniq -c | awk '{print $1,$10}') <(bedtools intersect -a genomic_gene_1000bpflanks.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wao | awk '$16==0{print $0}' | awk '{print "0 "$9}') > genomic_ins_pergeneflanks.txt


# random sample
cd /tmp/global2/ssushko/Vanessa_proj/6_gene_expansion/expanded_genes_cardui/flanks/random_sample
awk 'NR%2{printf "%s\t", $0; next}1' ../genome-wide/genomic_ins_pergeneflanks.sorted.txt > temp_flanks.twocols # move each gene to one row

shuf -n 77 --random-source=<(get_seeded_random 42) temp_flanks.twocols > temp_77_random_flanks_rs42.table
shuf -n 77 --random-source=<(get_seeded_random 102) temp_flanks.twocols > temp_77_random_flanks_rs102.table
shuf -n 77 --random-source=<(get_seeded_random 12) temp_flanks.twocols > temp_77_random_flanks_rs12.table
shuf -n 77 --random-source=<(get_seeded_random 53) temp_flanks.twocols > temp_77_random_flanks_rs53.table
shuf -n 77 --random-source=<(get_seeded_random 98) temp_flanks.twocols > temp_77_random_flanks_rs98.table

list='temp_77_random_flanks_rs42.table
temp_77_random_flanks_rs102.table
temp_77_random_flanks_rs12.table
temp_77_random_flanks_rs53.table
temp_77_random_flanks_rs98.table
'

for file in $list;
do
name=$(echo $file | cut -f2- -d '_')
echo 'insertions	gene' > $name
awk -F'\t' '{print $1 "\n" $2}' $file >> $name
done