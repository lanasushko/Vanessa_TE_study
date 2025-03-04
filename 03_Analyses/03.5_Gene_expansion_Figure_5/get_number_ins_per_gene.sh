# Calculate number of insertions of gene in V. cardui (expanded vs. genomic average)

### genes ###
# calculate n insertions per expanded gene
cd expanded_genes

cat <(bedtools intersect -a ../gained_genes_filtered_corr.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wa | sort | uniq -c | awk '{print $1,$5,$4-$3}') <(bedtools intersect -a ../gained_genes_filtered_corr.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wao | awk '$13==0{print $0}' | awk '{print "0 "$4,$3-$2}') > expanded_average_per_gene.txt

# calculate n insertions per gene in the genome
cd expanded_genes

cat <(bedtools intersect -a ../dashas_genes.gff -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wa | sort | uniq -c | awk '{print $1,$10,$6-$5}') <(bedtools intersect -a ../dashas_genes.gff -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wao | awk '$16==0{print $0}' | awk '{print "0 "$9,$5-$4}') > genomic_average_per_gene.txt

### gene flanks ###
# expanded
bedtools flank -i gained_genes_filtered_corr.bed -g genomes/ilVanCard2.1/V_card_chrnumbers.genome -b 1000 > flanks/expanded_gene_1000bpflanks.bed
cd flanks
bedtools intersect -a expanded_gene_1000bpflanks.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wa | cut -f4 | sort | uniq -c > expanded_ins_pergeneflanks.txt

cat <(bedtools intersect -a expanded_gene_1000bpflanks.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wa | sort | uniq -c | awk '{print $1,$5}') <(bedtools intersect -a expanded_gene_1000bpflanks.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wao | awk '$13==0{print $0}' | awk '{print "0 "$4}') > expanded_ins_pergeneflanks.txt


# genomic
bedtools flank -i dashas_genes.gff -g genomes/ilVanCard2.1/V_card_chrnumbers.genome -b 1000 > flanks/genomic_gene_1000bpflanks.bed
cd flanks

cat <(bedtools intersect -a genomic_gene_1000bpflanks.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wa | sort | uniq -c | awk '{print $1,$10}') <(bedtools intersect -a genomic_gene_1000bpflanks.bed -b ../V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed -wao | awk '$16==0{print $0}' | awk '{print "0 "$9}') > genomic_ins_pergeneflanks.txt

