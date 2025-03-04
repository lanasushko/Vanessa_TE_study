# ilVanCard2.2 self-alignment

# aligning with minimap
minimap2 -t 8 -ax asm5 --eqx genomes/ilVanCard2.2/GCA_905220365.2_genomic_genbank.fna genomes/ilVanCard2.2/GCA_905220365.2_genomic_genbank.fna > aln2.2.sam

# obtain a paf version for easier handling
paftools.js sam2paf aln2.2.sam > aln2.2.paf

# sort and index BAM
samtools view -@ 8 -b aln2.2.sam -o aln2.2.bam
samtools sort -@ 8 aln2.2.bam -o aln2.2_sorted.bam
samtools index aln2.2_sorted.bam

rm aln2.2.sam aln2.2.bam