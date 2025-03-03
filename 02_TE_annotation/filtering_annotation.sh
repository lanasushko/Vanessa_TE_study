# Filtering annotation produced by RepeatMasker

list = 'A_io_GCF_905147045.1_genomic_refseq.fna.out
V_atal_GCF_905147765.1_genomic_refseq.fna.out
V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.out
GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna.out
'


for file in $list;
do
# filter out 'Simple_repeat','Satellite','rRNA','snRNA','Low_complexity'
grep -v -E 'Simple_repeat|Satellite|rRNA|snRNA|Low_complexity' $file > $file.temp
# remove hits with <300 score >40 div <80 length
awk '$1>300 && $2<40 && ($7-$6)>80{print $0}' $file.temp > $file.filtered.out
done
