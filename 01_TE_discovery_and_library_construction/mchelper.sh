
# V. cardui 
python3 mchelper/MCHelper.py -r A -l RepeatLibraries/vcard_vatal_libraries_raw/v_cardui/vcardui_fulldenovo_forrpmnames.fasta -o RepeatLibraries/vcard_vatal_libraries_raw/MCHelper/card -g _ncbi_downloads/genomes/ilVanCard2.2/V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna -b BUSCOlepidoptera/lepidoptera_odb10.hmm --input_type fasta -a F -t 20


# V. atalanta
python3 mchelper/MCHelper.py -r A -l RepeatLibraries/4species_finalLibs/vatalanta_fulldenovo_forrpmnames.fasta -o RepeatLibraries/vcard_vatal_libraries_raw/MCHelper/atal -g _ncbi_downloads/genomes/ilVanAtal1.2/V_atal_GCF_905147765.1_genomic_refseq.fna -b BUSCOlepidoptera/lepidoptera_odb10.hmm --input_type fasta -a F -t 20


# Aglais io
python3 mchelper/MCHelper.py -r A -l RepeatLibraries/4species_finalLibs/Aglais_io_fulldenovo.fa -o RepeatLibraries/vcard_vatal_libraries_raw/MCHelper/aglaisio -g _ncbi_downloads/genomes/ilAglIoxx1.1/A_io_GCF_905147045.1_genomic_refseq.fna -b BUSCOlepidoptera/lepidoptera_odb10.hmm --input_type fasta -a F -t 20

# V. tameamea
python3 mchelper/MCHelper.py -r A -l RepeatLibraries/Aio_Vtameamae_libraries/V_tameamea_newgenome/Vanessa_tameamea_newg_fulldenovo.fa -o RepeatLibraries/Aio_Vtameamae_libraries/MCHelper/vantam_newg -g _ncbi_downloads/genomes/ilVanTame1/GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna -b BUSCOlepidoptera/lepidoptera_odb10.hmm --input_type fasta -a F -t 20