setwd("/tmp/global2/ssushko/Vanessa_proj/2_synteny/synteny_analyses")
library(gggenomes)
library(ggplot2)

#input genomes
vcard_genome=read.table('chr_start_end.cardui.txt')
vcard_genome=vcard_genome[,c(1,2,4)]
colnames(vcard_genome)=c('bin_id','seq_id','length')
vatal_genome=read.table('chr_start_end.atalanta.txt')
vatal_genome=vatal_genome[,c(1,2,4)]
colnames(vatal_genome)=c('bin_id','seq_id','length')
vtam_genome=read.table('chr_length_tameamea.txt')
colnames(vtam_genome)=c('bin_id','seq_id','length')
aio_genome=read.table('chr_length_aglaisio.txt')
colnames(aio_genome)=c('bin_id','seq_id','length')

sequences=rbind(vcard_genome,vatal_genome,vtam_genome,aio_genome)

# synthenic chrs
cardui=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,'Z','W')
atalanta=c(1,2,10,4,8,5,3,9,7,20,6,13,12,11,14,15,16,17,21,18,19,22,23,24,25,26,27,28,29,30,'Z','W')
aio=c(1,2,12,4,9,5,3,10,6,18,8,15,14,13,7,16,11,17,21,19,20,23,25,22,24,26,28,29,27,30,'Z')
tameamea=c(2,5,11,3,9,6,4,10,7,21,8,13,14,12,15,16,17,18,22,19,20,23,24,25,26,27,28,29,30,31,'Z','W')

cardui_chrs=c('vcchr1','vcchr2','vcchr3','vcchr4','vcchr5','vcchr6','vcchr7','vcchr8','vcchr9','vcchr10','vcchr11','vcchr12','vcchr13','vcchr14','vcchr15','vcchr16','vcchr17','vcchr18','vcchr19','vcchr20','vcchr21','vcchr22','vcchr23','vcchr24','vcchr25','vcchr26','vcchr27','vcchr28','vcchr29','vcchr30','vcchrZ','vcchrW')
atalanta_chrs=c('vachr1','vachr2','vachr10','vachr4','vachr8','vachr5','vachr3','vachr9','vachr7','vachr20','vachr6','vachr13','vachr12','vachr11','vachr14','vachr15','vachr16','vachr17','vachr21','vachr18','vachr19','vachr22','vachr23','vachr24','vachr25','vachr26','vachr27','vachr28','vachr29','vachr30','vachrZ','vachrW')
aio_chrs=c('aichr1','aichr2','aichr12','aichr4','aichr9','aichr5','aichr3','aichr10','aichr6','aichr18','aichr8','aichr15','aichr14','aichr13','aichr7','aichr16','aichr11','aichr17','aichr21','aichr19','aichr20','aichr23','aichr25','aichr22','aichr24','aichr26','aichr28','aichr29','aichr27','aichr30','aichrZ')
tameamea_chrs=c('vtchr2','vtchr5','vtchr11','vtchr3','vtchr9','vtchr6','vtchr4','vtchr10','vtchr7','vtchr21','vtchr8','vtchr13','vtchr14','vtchr12','vtchr15','vtchr16','vtchr17','vtchr18','vtchr22','vtchr19','vtchr20','vtchr23','vtchr24','vtchr25','vtchr26','vtchr27','vtchr28','vtchr29','vtchr30','vtchr31','vtchrZ','vtchrW')
order=c(cardui_chrs,atalanta_chrs,tameamea_chrs,aio_chrs)

sequences$seq_id=factor(sequences$seq_id, levels=order)
sequences=sequences[order(sequences$seq_id),]

# input links
library(pafr)
alignment_card_atal=read_paf('aln_regional_card-atal_chrnames.paf')
alignment_atal_tam=read_paf('aln_vatal-vtam_chrnames.paf')
alignment_tam_aio=read_paf('aln_vtam-aio_chrnames.paf')
alignment=rbind(alignment_card_atal,alignment_atal_tam,alignment_tam_aio)
# prim_alignment <- filter_secondary_alignments(alignment_card_atal)
regions=alignment[,c(6,8,9,1,3,4,5)]
colnames(regions)=c('seq_id','start','end','seq_id2','start2','end2','strand')

# filter regions
regions_filtered=regions[(regions$end-regions$start)>10000, ]

# add TE data
te.bpcounts=read.table("cardui.TEbp.per.window.table")
colnames(te.bpcounts)=c('seq_id','start','end','te.bp')
te.bpcounts <- te.bpcounts[order(te.bpcounts$te.bp),]

# add tandem repeat data
repeat.bpcounts=read.table('/tmp/global2/ssushko/Vanessa_proj/satellites/cardui/density/vcard_repeat.windows')
colnames(repeat.bpcounts)=c('seq_id','start','end','repeat.bp')
repeat.bpcounts <- repeat.bpcounts[order(repeat.bpcounts$repeat.bp),]


# plot
for (chromosome in 1:32) {
number=chromosome
chr=subset(sequences, sequences$seq_id==paste0('vcchr',cardui[number]) | sequences$seq_id==paste0('vachr',atalanta[number]) | sequences$seq_id==paste0('vtchr',tameamea[number]) | sequences$seq_id==paste0('aichr',aio[number]))
plot=gggenomes(seqs=chr,feats=te.bpcounts) +
  geom_seq() + geom_bin_label() + 
  # geom_seq_label(vjust=-0.5) + 
  xlim(-4*10^6,18*10^6) + theme_void() + ggtitle(chromosome)
# plot=plot + geom_coverage(aes(z=te.bp, color=te.bp), height=0.5, geom = "linerange",linewidth=0.3)  +
  # scale_color_gradient(limits=c(700,100000),guide="none")
plot=plot %>% add_links(regions_filtered) %>%
  sync() + geom_link()
name=paste0("c",chromosome)
assign(name, plot)
}

# chr28=subset(sequences, sequences$seq_id==paste0('vcchr',cardui[28]) | sequences$seq_id==paste0('vachr',atalanta[28]) | sequences$seq_id==paste0('vtchr',tameamea[28]) | sequences$seq_id==paste0('aichr',aio[28]))
# c28=gggenomes(seqs=chr28,feats=te.bpcounts) +
#   geom_seq() + geom_bin_label() + 
#   # geom_seq_label(vjust=-0.5) + 
#   xlim(-4*10^6,18*10^6) + theme_void() + ggtitle(chromosome)
# # plot=plot + geom_coverage(aes(z=te.bp, color=te.bp), height=0.5, geom = "linerange",linewidth=0.3)  +
# # scale_color_gradient(limits=c(700,100000),guide="none")
# c28=c28 %>% add_links(regions_filtered) %>%
#   flip(3:4) + geom_link()
# c28

library(grid)
library(gridExtra)
grid.arrange(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,ncol=4)

# Assume you have the following 32 plots
plots <- list(
  c1, c2, c3, c4, c5, c6, c7, c8, c9, c10,
  c11, c12, c13, c14, c15, c16, c17, c18, c19, c20,
  c21, c22, c23, c24, c25, c26, c27, c28, c29, c30,
  c31, c32
)

# Reorder plots in a column-major order
#2 columns reorder
reordered_plots <- plots[c(
  1, 17, 2, 18, 3, 19, 4, 20, 5, 21,
  6, 22, 7, 23, 8, 24, 9, 25, 10, 26,
  11, 27, 12, 28, 13, 29, 14, 30, 15, 31,
  16, 32
)]

#4 column reorder
reordered_plots <- plots[c(
  1, 9, 17, 25, 2, 10, 18, 26, 3, 11, 19, 27,
  4, 12, 20, 28, 5, 13, 21, 29, 6, 14, 22, 30,
  7, 15, 23, 31, 8, 16, 24, 32
)]

grid.arrange(grobs=reordered_plots, ncol=4)

#### individual density plots
te.bpcounts$seq_id=factor(te.bpcounts$seq_id, levels=cardui_chrs)
repeat.bpcounts$seq_id=factor(repeat.bpcounts$seq_id, levels=cardui_chrs)

# TEs
ggplot(data=te.bpcounts, aes(x=start,y=te.bp, group=seq_id, fill=te.bp)) +
  geom_col() + facet_wrap(~seq_id,ncol=4) + theme_classic() + scale_fill_gradient(low = "#b50707", high = "yellow", na.value = NA)

# Repeats
# stable y
ggplot(data=repeat.bpcounts, aes(x=start,y=repeat.bp, group=seq_id, fill=repeat.bp)) +
  geom_col() + facet_wrap(~seq_id,ncol=4) + theme_classic() + scale_fill_gradient(low = "#b50707", high = "yellow", na.value = NA)

# free y
ggplot(data=repeat.bpcounts, aes(x=start,y=repeat.bp, group=seq_id, fill=repeat.bp)) +
  geom_col() + facet_wrap(~seq_id,ncol=4, scales='free_y') + theme_classic() + scale_fill_gradient(low = "#b50707", high = "yellow", na.value = NA)



