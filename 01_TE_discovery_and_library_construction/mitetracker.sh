# MITE Tracker runs to detect MITE transposons

# # V. cardui
python3 -m MITETracker -g genomes/ilVanCard2.2/GCA_905220365.2_genomic_genbank.fna -w 3 -j mt-cardui

# V. atalanta
python3 -m MITETracker -g genomes/ilVanAtal1.2/V_atal_GCF_905147765.1_genomic_refseq.fna -w 3 -j mt-atalanta

# V. tameamea
python3 -m MITETracker -g genomes/ilVanTame1/GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna -w 3 -j mt-tameamea

# A. io
python3 -m MITETracker -g genomes/ilAglIoxx1.1/GCF_905147045.1_genomic_refseq.fna -w 3 -j mt-aio
