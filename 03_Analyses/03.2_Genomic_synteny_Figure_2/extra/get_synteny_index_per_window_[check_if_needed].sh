# ### OPTION 1: TAKING INTO ACCOUNT DIVERGENCE (but there are overlaps in alignments)

# # get the regions of mapping in cardui
# cut -f6,8,9 aln_regional_card-atal_chrnames.paf > aln_chrse_cardui.table

# # get divergences
# grep -oP '(?<=dv:f:).{6}' aln_regional_card-atal_chrnames.paf > div_aln_cardui.table 

# # join
# paste aln_chrse_cardui.table div_aln_cardui.table > aln_chrsediv_cardui.table

# # get overlaps
# bedtools intersect -a /tmp/global2/ssushko/Vanessa_proj/coverage_and_enrichment/TEcoverage/cardui.windows.bed -b /tmp/global2/ssushko/Vanessa_proj/alignments/align_vcard-vatal/aln_chrsediv_cardui.table -wao > window_intersect.table

# # get indexes of synteny as %synteny*divergence per each mapping region
# awk '{print $1"_"$2"_"$3,$9/100000*100*(1-$8)}' window_intersect.table | awk '{arr[$1]+=$2} END {for (i in arr) {print i,arr[i]}}' | sed 's/_/ /g' | sed 's/ /\t/g' | sort -k1,1 -k2,2 > window_synteny.table

# # join the synteny and TEbp files
# paste window_synteny.table <(cut -f4 /tmp/global2/ssushko/Vanessa_proj/coverage_and_enrichment/TEcoverage/cardui.TEbp.per.window.table.sorted) > window_synteny_TEbp.table





# ### OPTION 2: REMOVING OVERLAPS (and taking into account divergence)
# cp ../option1/aln_chrsediv_cardui.table .
# sort -k1,1 -k2,2n aln_chrsediv_cardui.table > aln_chrsediv_cardui.table.sorted

# # merge bed intervals and perform mean on divergence
# bedtools merge -i aln_chrsediv_cardui.table.sorted -c 4 -o mean > aln_chrsediv_cardui.table.sorted.merged

# # get overlaps
# bedtools intersect -a /tmp/global2/ssushko/Vanessa_proj/coverage_and_enrichment/TEcoverage/cardui.windows.bed -b aln_chrsediv_cardui.table.sorted.merged -wao > window_intersect.table

# # get indexes of synteny as %synteny*divergence per each mapping region
# awk '{print $1"_"$2"_"$3,$9/100000*100*(1-$8)}' window_intersect.table | awk '{arr[$1]+=$2} END {for (i in arr) {print i,arr[i]}}' | sed 's/_/ /g' | sed 's/ /\t/g' | sort -k1,1 -k2,2 > window_synteny.table


# ### OPTION 3: ONLY "SAME"-CHROMOSOME MATCHES AND DISCARD SMALL ALIGNMENTS (and removing overlaps and taking into account divergence)

# get windows
# bedtools makewindows -g /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilAglIoxx1.1/A_io_GCF_905147045.1.genome -w 100000 -i winnum > aio.windows
# /tmp/global2/ssushko/Vanessa_proj/scripts/change_chr_names/tospecieschr/sed_chromosome_names_aio.sh aio.windows

# bedtools makewindows -g /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanAtal1.2/V_atal.genome -w 100000 -i winnum > vatal.windows
# /tmp/global2/ssushko/Vanessa_proj/scripts/change_chr_names/tospecieschr/sed_chromosome_names_atalanta.sh vatal.windows

# bedtools makewindows -g /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanCard2.2/V_card.genome -w 100000 -i winnum > vcard.windows
# /tmp/global2/ssushko/Vanessa_proj/scripts/change_chr_names/tospecieschr/sed_chromosome_names_2.2cardui.sh vcard.windows

# bedtools makewindows -g /tmp/global2/ssushko/Vanessa_proj/_ncbi_downloads/genomes/ilVanTame1/chr_length_vtam.txt -w 100000 -i winnum > vtam.windows
# /tmp/global2/ssushko/Vanessa_proj/scripts/change_chr_names/tospecieschr/sed_chromosome_names_vtam.shv vtam.windows

paf=$1
windowsfile=$2
name=$(basename $paf)

# remove supplementary alignments and small primary alignments (<3kb)
awk '$13=="tp:A:P" && $11>3000{print $0}' $paf > temp_$name.primary.mt3kb.paf

# get only the synthenic chromosome alignments
python3 /tmp/global2/ssushko/Vanessa_proj/scripts/synteny_analysis/filter_alignment_by_chrsynteny_v2.py temp_$name.primary.mt3kb.paf temp_$name.primary.mt3kb.synchr.paf

# get the regions of mapping in cardui
cut -f6,8,9 temp_$name.primary.mt3kb.synchr.paf > temp_regions_$name.table

# get divergences
grep -oP '(?<=dv:f:).{6}' temp_$name.primary.mt3kb.synchr.paf > temp_div_$name.table 

# join
paste temp_regions_$name.table temp_div_$name.table > regions_div_$name.table

# clean
rm temp_*

# sort 
sort -k1,1 -k2,2n regions_div_$name.table > temp_regions_div_$name.sorted.table

# merge bed intervals and perform mean on divergence
bedtools merge -i temp_regions_div_$name.sorted.table -c 4 -o mean > temp_regions_div_$name.sorted.merged


# get overlaps
bedtools intersect -a $windowsfile -b temp_regions_div_$name.sorted.merged -wao > temp_window_regions_div_$name.sorted.merged_intersect.table

# get indexes of synteny as %synteny*divergence per each mapping region
awk '{print $1"_"$2"_"$3,$9/100000*100*(1-$8)}' temp_window_regions_div_$name.sorted.merged_intersect.table | awk '{arr[$1]+=$2} END {for (i in arr) {print i,arr[i]}}' | sed 's/_/ /g' | sed 's/ /\t/g' | sort -k1,1V -k2,2 | grep -vE 'NW|chrUn|chrWu|JAXCL'> window_$name.synteny.table

# clean 
rm temp_*


