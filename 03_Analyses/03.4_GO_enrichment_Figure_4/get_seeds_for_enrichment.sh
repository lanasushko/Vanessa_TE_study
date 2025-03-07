# Checking for enrichment of certain GO groups among genes where MITEs are inserted inside gene body and flanking regions

# MITE enrichment
## here I took all the genes with MITE insertions, independently of how many insertions there are and how much do they cover the genes or which parts of the gene exon/intron/cds they hit

# Get gene seeds that contain MITE insertions for each species
grep 'CLASSII/MITE' V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.out.bed > V_card_MITEs.bed
bedtools intersect -a cardui_codinggenes.gff -b V_card_MITEs.bed | cut -f9 | cut -f3 -d';' | sed 's/Name=//' | sort | uniq -c | sort -k1.1nr > genes_with_MITEs.uniqc.txt
grep 'ID=cds-' cardui_geneannot.gff | cut -f9 | sed 's/ID=cds-//' | sed -e 's/;.*;gene=/ /' | sed -e 's/;.*//' | sort -u > cardui_proteinid-geneid.txt
grep -f MITE_inside_gene/genes_with_MITEs.idlist cardui_proteinid-geneid.txt | cut -f1 -d' ' > proteins_affected_by_MITEs.list
grep -f proteins_affected_by_MITEs.list geneannot/cardui/cardui.emapper.annotations | cut -f2 | sort -u > seeds_with_MITEs.list

grep 'CLASSII/MITE' V_atal_GCF_905147765.1_genomic_refseq.fna.filtered.nonoverlapping.bed > V_atal_MITE.bed
bedtools intersect -a atalanta_codinggenes.gff -b V_atal_MITE.bed | cut -f9 | cut -f3 -d';' | sed 's/Name=//' | sort | uniq -c | sort -k1.1nr > genes_with_MITEs_inside.txt
bedtools intersect -a atalanta_codinggenes.gff -b V_atal_MITE.bed | cut -f9 | cut -f3 -d';' | sed 's/Name=//' | sort -u > genes_with_MITEs.idlist
grep 'ID=cds-' atalanta_geneannot.gff | cut -f9 | sed 's/ID=cds-//' | sed -e 's/;.*;gene=/ /' | sed -e 's/;.*//' | sort -u > atalanta_proteinid-geneid.txt
grep -f genes_with_MITEs.idlist atalanta_proteinid-geneid.txt | cut -f1 -d' ' > proteins_affected_by_MITEs.list
grep -f proteins_affected_by_MITEs.list geneannot/atalanta/atalanta.emapper.annotations | cut -f2 | sort -u > seeds_with_MITEs.list

grep 'CLASSII/MITE' A_io_GCF_905147045.1_genomic_refseq.fna.filtered.out.nonoverlapping.bed > A_io_MITE.bed
bedtools intersect -a aglaisio_codinggenes.gff -b A_io_MITE.bed | cut -f9 | cut -f3 -d';' | sed 's/Name=//' | sort -u > genes_with_MITEs.idlist
grep 'ID=cds-' aglaisio_geneannot.gff | cut -f9 | sed 's/ID=cds-//' | sed -e 's/;.*;gene=/ /' | sed -e 's/;.*//' | sort -u > aio_proteinid-geneid.txt
grep -f genes_with_MITEs.idlist aio_proteinid-geneid.txt | cut -f1 -d' ' > proteins_affected_by_MITEs.list
grep -f proteins_affected_by_MITEs.list geneannot/aglaisio/aglaisio.emapper.annotations | cut -f2 | sort -u > seeds_with_MITEs.list

grep 'CLASSII/MITE' GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna.filtered.out.nonoverlapping.bed > V_tame_MITE.bed
bash sed_chromosome_names_tameamea.sh V_tame_MITE.bed
bedtools intersect -a tameamea_codinggenes.gff -b V_tame_MITE.bed | cut -f9 | cut -f3 -d';' | sed 's/Name=//' | sort -u > genes_with_MITEs.idlist
grep 'ID=cds-' tameamea_geneannot.gff | cut -f9 | sed 's/ID=cds-//' | sed -e 's/;.*;gene=/ /' | sed -e 's/;.*//' | sort -u > tameamea_proteinid-geneid.txt
grep -f genes_with_MITEs.idlist tameamea_proteinid-geneid.txt -w | cut -f1 -d' ' > proteins_affected_by_MITEs.list
grep -f proteins_affected_by_MITEs.list geneannot/tameamea/tameamea.emapper.annotations | cut -f2 | sort -u > seeds_with_MITEs.list
 