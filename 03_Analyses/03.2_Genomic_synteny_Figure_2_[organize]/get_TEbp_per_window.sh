
w_card=/tmp/global2/ssushko/Vanessa_proj/2_synteny/synteny_analyses/genomic_windows/vcard.windows
w_atal=/tmp/global2/ssushko/Vanessa_proj/2_synteny/synteny_analyses/genomic_windows/vatal.windows
w_tam=/tmp/global2/ssushko/Vanessa_proj/2_synteny/synteny_analyses/genomic_windows/vtam.windows
w_aio=/tmp/global2/ssushko/Vanessa_proj/2_synteny/synteny_analyses/genomic_windows/aio.windows

TE_card=/tmp/global2/ssushko/Vanessa_proj/1_repeat_annotation/with_NEWtameameagenome/v_card/V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.nonoverlapping.chrnamessp.bed
TE_atal=/tmp/global2/ssushko/Vanessa_proj/1_repeat_annotation/with_NEWtameameagenome/v_atal/V_atal_GCF_905147765.1_genomic_refseq.fna.filtered.nonoverlapping.chrnamessp.bed
TE_tam=/tmp/global2/ssushko/Vanessa_proj/1_repeat_annotation/with_NEWtameameagenome/v_tam/GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna.filtered.out.nonoverlapping.chrnamessp.bed
TE_aio=/tmp/global2/ssushko/Vanessa_proj/1_repeat_annotation/with_NEWtameameagenome/a_io/A_io_GCF_905147045.1_genomic_refseq.fna.filtered.out.nonoverlapping.chrnamessp.bed


# get overlaps
bedtools intersect -a $w_card -b $TE_card -wao | awk '{print $1"_"$2"_"$3,$11}' | awk '{arr[$1]+=$2} END {for (i in arr) {print i,arr[i]}}' | sed 's/_/ /g' | sed 's/ /\t/g' | sort -k1,1V -k2,2 | grep -vE 'NW|chrUn|chrWu|JAXCL' > vcard_TE.windows
bedtools intersect -a $w_atal -b $TE_atal -wao | awk '{print $1"_"$2"_"$3,$11}' | awk '{arr[$1]+=$2} END {for (i in arr) {print i,arr[i]}}' | sed 's/_/ /g' | sed 's/ /\t/g' | sort -k1,1V -k2,2 | grep -vE 'NW|chrUn|chrWu|JAXCL' > vatal_TE.windows
bedtools intersect -a $w_tam -b $TE_tam -wao | awk '{print $1"_"$2"_"$3,$11}' | awk '{arr[$1]+=$2} END {for (i in arr) {print i,arr[i]}}' | sed 's/_/ /g' | sed 's/ /\t/g' | sort -k1,1V -k2,2 | grep -vE 'NW|chrUn|chrWu|JAXCL' > vtam_TE.windows
bedtools intersect -a $w_aio -b $TE_aio -wao | awk '{print $1"_"$2"_"$3,$11}' | awk '{arr[$1]+=$2} END {for (i in arr) {print i,arr[i]}}' | sed 's/_/ /g' | sed 's/ /\t/g' | sort -k1,1V -k2,2 | grep -vE 'NW|chrUn|chrWu|JAXCL' > aio_TE.windows
