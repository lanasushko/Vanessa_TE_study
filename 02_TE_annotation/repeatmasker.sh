# Annotate TEs in the 4 genomes using the non-redundant TE library

# V. cardui -- this run was performed without Mitochondrial contig
RepeatMasker -gff -a -no_is -nolow -norna -s -lib RepeatLibraries/4species_finalLibs/withMCHelper_newtame/clustered_full_lib.fa_linearized_rmnonTE_stdnames.fa -pa 5 -dir . genomes/ilVanCard2.2/V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna

# V. atalanta
RepeatMasker -gff -a -no_is -nolow -norna -s -lib RepeatLibraries/4species_finalLibs/withMCHelper_newtame/clustered_full_lib.fa_linearized_rmnonTE_stdnames.fa -pa 5 -dir . genomes/ilVanAtal1.2/V_atal_GCF_905147765.1_genomic_refseq.fna

# V. tameamea
RepeatMasker -gff -a -no_is -nolow -norna -s -lib RepeatLibraries/4species_finalLibs/withMCHelper_newtame/clustered_full_lib.fa_linearized_rmnonTE_stdnames.fa -pa 5 -dir . genomes/ilVanTame1/GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna

# A. io
RepeatMasker -gff -a -no_is -nolow -norna -s -lib RepeatLibraries/4species_finalLibs/withMCHelper_newtame/clustered_full_lib.fa_linearized_rmnonTE_stdnames.fa -pa 5 -dir . genomes/ilAglIoxx1.1/A_io_GCF_905147045.1_genomic_refseq.fna

