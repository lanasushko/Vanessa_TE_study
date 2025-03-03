#!/bin/bash

#$ -w e
#$ -l h_vmem=24G
#$ -l h_rt=24:00:00
#$ -o /tmp/global2/ssushko/Vanessa_proj/out
#$ -e /tmp/global2/ssushko/Vanessa_proj/out
#$ -j yes
#$ -m ea
#$ -pe parallel 8 
#$ -M svitlana.sushko@tuebingen.mpg.de

## cardui, aio ##
# cd /tmp/global2/ssushko/Vanessa_proj/align_vcard-aio

# source activate minimap

# minimap2 -t 8 -x asm20 /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanCard2.2/GCA_905220365.2_genomic_genbank.fna /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilAglIoxx1.1/A_io_GCF_905147045.1_genomic_refseq.fna > aln_vcard-aglais.paf

## cardui, tameamea ##
# cd /tmp/global2/ssushko/Vanessa_proj/alignments/align_vcard-vtam

# source activate minimap

# minimap2 -t 8 -x asm20 /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanCard2.2/GCA_905220365.2_genomic_genbank.fna /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanTame1/GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna > aln_vcard-vtam.paf

## atalanta, tameamea ##
# cd /tmp/global2/ssushko/Vanessa_proj/alignments/align_vatal-vtam

# source activate minimap

# minimap2 -t 8 -x asm20 /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanAtal1.2/V_atal_GCF_905147765.1_genomic_refseq.fna /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanTame1/GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna > aln_vatal-vtam.paf

## tameamea, aio ##
# cd /tmp/global2/ssushko/Vanessa_proj/alignments/align_vtam_aio

# source activate minimap

# minimap2 -t 8 -x asm20 /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanTame1/GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilAglIoxx1.1/A_io_GCF_905147045.1_genomic_refseq.fna > aln_vtam-aio.paf

## cardui, atalanta ##

# cd /tmp/global2/ssushko/Vanessa_proj/align_vcard-vatal

# source activate minimap

# minimap2 -t 8 -x asm20 /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanCard2.2/GCA_905220365.2_genomic_genbank.fna /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanAtal1.2/GCF_905147765.1_genomic_refseq.fna > aln3.paf


## Preparing the rest of alignments ##

cd /tmp/global2/ssushko/Vanessa_proj/alignments
genome_vcard=/tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanCard2.2/GCA_905220365.2_genomic_genbank.fna
genome_vatal=/tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanAtal1.2/V_atal_GCF_905147765.1_genomic_refseq.fna
genome_vtam=/tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanTame1/GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna
genome_aio=/tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilAglIoxx1.1/A_io_GCF_905147045.1_genomic_refseq.fna

## atalanta, cardui ##
alname='vatal-vcard'

mkdir align_$alname
cd align_$alname

source activate minimap_gmap

minimap2 -t 8 -x asm20 $genome_vatal $genome_vcard > aln_$alname.paf

cd ../

## atalanta, aio ##
alname='vatal-aio'

mkdir align_$alname
cd align_$alname

minimap2 -t 8 -x asm20 $genome_vatal $genome_aio > aln_$alname.paf

cd ../

## tameamea, cardui ##
alname='vtam-vcard'

mkdir align_$alname
cd align_$alname

minimap2 -t 8 -x asm20 $genome_vtam $genome_vcard > aln_$alname.paf &>minimap.out&

cd ../

## tameamea, atalanta ##
alname='vtam-vatal'

mkdir align_$alname
cd align_$alname

minimap2 -t 8 -x asm20 $genome_vtam $genome_vatal > aln_$alname.paf

cd ../

## aio, cardui ##
alname='aio-vcard'

mkdir align_$alname
cd align_$alname

minimap2 -t 8 -x asm20 $genome_aio $genome_vcard > aln_$alname.paf  &>minimap.out&

cd ../

## aio, atalanta ##
alname='aio-vatal'

mkdir align_$alname
cd align_$alname

minimap2 -t 8 -x asm20 $genome_aio $genome_vatal > aln_$alname.paf  &>minimap.out&

cd ../

## aio, tameamea ##
alname='aio-vtam'

mkdir align_$alname
cd align_$alname

minimap2 -t 8 -x asm20 $genome_aio $genome_vtam > aln_$alname.paf  &>minimap.out&

cd ../