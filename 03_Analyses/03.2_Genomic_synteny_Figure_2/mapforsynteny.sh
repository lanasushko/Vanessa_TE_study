# Mapping genomes for synteny analysis

## cardui, atalanta ##
minimap2 -t 8 -x asm20 genomes/ilVanCard2.2/GCA_905220365.2_genomic_genbank.fna genomes/ilVanAtal1.2/GCF_905147765.1_genomic_refseq.fna > aln_vcard-vatal.paf

## atalanta, tameamea ##
minimap2 -t 8 -x asm20 genomes/ilVanAtal1.2/V_atal_GCF_905147765.1_genomic_refseq.fna genomes/ilVanTame1/GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna > aln_vatal-vtam.paf

## tameamea, aio ##
minimap2 -t 8 -x asm20 genomes/ilVanTame1/GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna genomes/ilAglIoxx1.1/A_io_GCF_905147045.1_genomic_refseq.fna > aln_vtam-aio.paf
