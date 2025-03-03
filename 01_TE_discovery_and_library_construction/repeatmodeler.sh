# RepeatModeler runs to identify transposons in the 4 buttefly genomes

## RepeatModeler v2.0.4
# V. cardui
BuildDatabase -name Vanessa_cardui genomes/ilVanCard2.2/GCA_905220365.2_genomic_genbank.fna
RepeatModeler -database Vanessa_cardui -threads 8 -LTRStruct

# V. atalanta
BuildDatabase -name Vanessa_atalanta genomes/ilVanAtal1.2/V_atal_GCF_905147765.1_genomic_refseq.fna
RepeatModeler -database Vanessa_atalanta -threads 8 -LTRStruct

## RepeatModeler v2.0.1
# V. tameamea
BuildDatabase -name Vanessa_tameamea genomes/ilVanTame1/GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna
RepeatModeler -database Vanessa_tameamea -pa 5 -LTRStruct

# A. io
BuildDatabase -name Aglais_io genomes/ilAglIoxx1.1/GCF_905147045.1_genomic_refseq.fna
RepeatModeler -database Aglais_io -pa 5 -LTRStruct