setwd("~/Documents/PROJECTS/3_IBB_vanessa/New_repeat_annotations/commonlib4sp_withnewtameameagenome")

# data input
aio=read.csv(file = 'summaries/a_io_mch_classes.sum', header=T, sep=',',na.strings=c("","NA"))
vcard=read.csv(file = 'summaries/v_card_mch_classes.sum', header=T, sep=',',na.strings=c("","NA"))
vatal=read.csv(file = 'summaries/v_atal_mch_classes.sum', header=T, sep=',',na.strings=c("","NA"))
vtam=read.csv(file = 'summaries/v_tam_mch_classes.sum', header=T, sep=',',na.strings=c("","NA"))

# make genome size vectors for a more precise % calculation
aiogs=384156914
vcardgs=424813635
vatalgs=370420815
vtamgs=363368026

#join data
aio$X.masked=aio$bpMasked/aiogs*100
vcard$X.masked=vcard$bpMasked/vcardgs*100
vatal$X.masked=vatal$bpMasked/vatalgs*100
vtam$X.masked=vtam$bpMasked/vtamgs*100

datafull=rbind(aio, vcard, vatal, vtam)

# remove unused rows
datafull=datafull[!(datafull$Class %in% "total_interspersed"),]
datafull=datafull[!(datafull$Class %in% "Total"),]

# add TEtype column with both class and superfamily information
datafull$TEtype <- paste(datafull$Class,datafull$Superfamily)

# change sample names to plottable names
datafull$Sample=gsub("A_io_GCF_905147045.1_genomic_refseq.fna.filtered.out", "A_io", datafull$Sample)
datafull$Sample=gsub("V_card_GCA_905220365.2_genomic_genbank_noMTchr.fna.filtered.out", "V_cardui", datafull$Sample)
datafull$Sample=gsub("V_atal_GCF_905147765.1_genomic_refseq.fna.filtered.out", "V_atalanta", datafull$Sample)
datafull$Sample=gsub("GCA_037043105.1_ilVanTame1_primary_haplotype_genomic.fna.filtered.out", "V_tameamea", datafull$Sample)

# species as factors for plotting in this order
datafull$Sample <- factor(datafull$Sample, levels=c("A_io", "V_cardui", "V_atalanta","V_tameamea"))

# add additional columns to show order and superfamily information separately
datafull=within(datafull, TE<-data.frame(do.call('rbind', strsplit(as.character(datafull$Superfamily), '/', fixed=TRUE))))

# to aggregate the % values of each TE type
df1=aggregate(X.masked~Sample+TE$X1, datafull, 'sum')

#only if making pie charts
# datafull$X.masked=as.numeric(datafull$X.masked)
# datafull[nrow(datafull) + 1,] = c("A_io","Non-TE","","","",100-sum(datafull[datafull$Sample=='A_io',]$X.masked),"Non-TE")
# datafull$X.masked=as.numeric(datafull$X.masked)
# datafull[nrow(datafull) + 1,] = c("V_card","Non-TE","","","",100-sum(datafull[datafull$Sample=='V_card',]$X.masked),"Non-TE")
# datafull$X.masked=as.numeric(datafull$X.masked)
# datafull[nrow(datafull) + 1,] = c("V_atal","Non-TE","","","",100-sum(datafull[datafull$Sample=='V_atal',]$X.masked),"Non-TE")
# datafull$X.masked=as.numeric(datafull$X.masked)
# datafull[nrow(datafull) + 1,] = c("V_tam","Non-TE","","","",100-sum(datafull[datafull$Sample=='V_tam',]$X.masked),"Non-TE")
# datafull$X.masked=as.numeric(datafull$X.masked)

# open libraries and define colors
library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
col_vector=col_vector[7:60]



###### PLOTS #######

# general overview all TE classes in all species
ggplot(data=df1,mapping = aes(y=Sample, x=X.masked, fill=`TE$X1`)) +
  geom_col(position='stack') +
  theme_bw() +
  scale_fill_manual(values=col_vector) +
  geom_text(label=df1$X.masked,check_overlap = F, position = position_stack(vjust = 0.5)) +
  ylab('% masked') +
  theme_classic()

# LINEs

# the next lines should only be ran when making the LINE plots
# datafull$TE$X2[datafull$TE$X2=="LINE"]<-"Unknown"
# datafull$TE$X2[datafull$TE$X2=="CRE"]<-"Other"
# datafull$TE$X2[datafull$TE$X2=="Dong-R4"]<-"Other"
# datafull$TE$X2[datafull$TE$X2=="Jockey-I"]<-"Other"
# datafull$TE$X2[datafull$TE$X2=="L1"]<-"Other"
# datafull$TE$X2[datafull$TE$X2=="R2"]<-"Other"
# datafull$TE$X2[datafull$TE$X2=="RTE-X"]<-"Other"
# datafull$TE$X2[datafull$TE$X2=="Proto2"]<-"Other"
# datafull$TE$X2[datafull$TE$X2=="Penelope"]<-"Other"
# datafull$TE$X2[datafull$TE$X2=="R2-NeSL"]<-"Other"
# datafull$TE$X2 <- factor(datafull$TE$X2, levels=c("CR1", "I", "L2","R1","R2","RTE","Unknown","Other"))

library(ggbreak)
vcardline=ggplot(data=subset(datafull, grepl("LINE", Superfamily) & grepl("V_cardui", Sample)) ,mapping = aes(y=TE$X2, x=X.masked, label=X.masked)) +
  geom_col(position='stack') +
  theme_classic() +
  scale_fill_manual(values=col_vector) +  xlab('% masked') +
  scale_y_discrete(limits=rev) + ylab('Superfamily') +
  scale_x_break(c(4.2, 6)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())
  #geom_text(check_overlap = F, position = position_stack(vjust = 0.5))

aioline=ggplot(data=subset(datafull, grepl("LINE", Superfamily) & grepl("A_io", Sample)) ,mapping = aes(y=TE$X2, x=X.masked, label=X.masked)) +
  geom_col(position='stack') +
  theme_classic() +
  scale_fill_manual(values=col_vector) +  xlab('% masked') +
  scale_y_discrete(limits=rev) + ylab('Superfamily') +
  scale_x_break(c(3, 8.5)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())

vatalline=ggplot(data=subset(datafull, grepl("LINE", Superfamily) & grepl("V_atalanta", Sample)) ,mapping = aes(y=TE$X2, x=X.masked, label=X.masked)) +
  geom_col(position='stack') +
  theme_classic() +
  scale_fill_manual(values=col_vector) +  xlab('% masked') +
  scale_y_discrete(limits=rev) + ylab('Superfamily') +
  scale_x_break(c(1.3, 3.5)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())

vtamline=ggplot(data=subset(datafull, grepl("LINE", Superfamily) & grepl("V_tameamea", Sample)) ,mapping = aes(y=TE$X2, x=X.masked, label=X.masked)) +
  geom_col(position='stack') +
  theme_classic() +
  scale_fill_manual(values=col_vector) +  xlab('% masked') +
  scale_y_discrete(limits=rev) + ylab('Superfamily') +
  scale_x_break(c(1.3, 3.5)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())

# LTRs

# the next lines should only be ran when making the LINE plots
datafull$TE$X2[datafull$TE$X2=="LTR"]<-"Unknown"

ggplot(data=subset(datafull, grepl("LTR", Superfamily) & grepl("V_cardui", Sample)) ,mapping = aes(y=TE$X2, x=X.masked, label=X.masked)) +
  geom_col(position='stack') +
  theme_classic() +
  scale_fill_manual(values=col_vector) +  xlab('% masked') +
  scale_y_discrete(limits=rev) + ylab('Superfamily') +
  scale_x_break(c(0.25, 3)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())
#geom_text(check_overlap = F, position = position_stack(vjust = 0.5))

ggplot(data=subset(datafull, grepl("LTR", Superfamily) & grepl("A_io", Sample)) ,mapping = aes(y=TE$X2, x=X.masked, label=X.masked)) +
  geom_col(position='stack') +
  theme_classic() +
  scale_fill_manual(values=col_vector) +  xlab('% masked') +
  scale_y_discrete(limits=rev) + ylab('Superfamily') +
  scale_x_break(c(0.2, 1.7)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())

ggplot(data=subset(datafull, grepl("LTR", Superfamily) & grepl("V_atalanta", Sample)) ,mapping = aes(y=TE$X2, x=X.masked, label=X.masked)) +
  geom_col(position='stack') +
  theme_classic() +
  scale_fill_manual(values=col_vector) +  xlab('% masked') +
  scale_y_discrete(limits=rev) + ylab('Superfamily') +
  scale_x_break(c(0.2, 1.7)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())

ggplot(data=subset(datafull, grepl("LTR", Superfamily) & grepl("V_tameamea", Sample)) ,mapping = aes(y=TE$X2, x=X.masked, label=X.masked)) +
  geom_col(position='stack') +
  theme_classic() +
  scale_fill_manual(values=col_vector) +  xlab('% masked') +
  scale_y_discrete(limits=rev) + ylab('Superfamily') +
  scale_x_break(c(0.2, 1.5)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())

# col_vector=col_vector[8:30]
# ggplot(data=subset(datafull, grepl("LTR", Superfamily)),mapping = aes(x=Sample, y=X.masked, fill=Superfamily, label=X.masked)) +
#   geom_col(position='stack') +
#   theme_bw() +
#   scale_fill_manual(values=col_vector) +
#   geom_text(check_overlap = F, position = position_stack(vjust = 0.5)) +
#   ylab('% masked') 

# TIRs

# the next lines should only be ran when making the LINE plots
datafull$TE$X2[datafull$TE$X2=="TIR"]<-"Unknown"
datafull$TE$X2[datafull$TE$X2=="CMC"]<-"Other"
datafull$TE$X2[datafull$TE$X2=="CMC-Chapaev"]<-"Other"
datafull$TE$X2[datafull$TE$X2=="HAT-Pegasus"]<-"Other"
datafull$TE$X2[datafull$TE$X2=="Kolobok-Hydra"]<-"Other"
datafull$TE$X2[datafull$TE$X2=="Kolobok-T2"]<-"Other"
datafull$TE$X2[datafull$TE$X2=="MULE"]<-"Other"
datafull$TE$X2[datafull$TE$X2=="MULE-NOF"]<-"Other"
datafull$TE$X2[datafull$TE$X2=="TcMar-Tc4"]<-"Other"
datafull$TE$X2[datafull$TE$X2=="Zator"]<-"Other"
datafull$TE$X2 <- factor(datafull$TE$X2, levels=c("Harbinger", "HAT", "MERLIN","P","PIFHARBINGER","PIGGYBAC","TcMar","Unknown","Other"))


ggplot(data=subset(datafull, grepl("TIR", Superfamily) & grepl("V_cardui", Sample)) ,mapping = aes(y=TE$X2, x=X.masked, label=X.masked)) +
  geom_col(position='stack') +
  theme_classic() +
  scale_fill_manual(values=col_vector) +  xlab('% masked') +
  scale_y_discrete(limits=rev) + ylab('Superfamily') +
  scale_x_break(c(0.2, 1.2)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())
#geom_text(check_overlap = F, position = position_stack(vjust = 0.5))

ggplot(data=subset(datafull, grepl("TIR", Superfamily) & grepl("A_io", Sample)) ,mapping = aes(y=TE$X2, x=X.masked, label=X.masked)) +
  geom_col(position='stack') +
  theme_classic() +
  scale_fill_manual(values=col_vector) +  xlab('% masked') +
  scale_y_discrete(limits=rev) + ylab('Superfamily') +
  scale_x_break(c(0.2, 1)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())

ggplot(data=subset(datafull, grepl("TIR", Superfamily) & grepl("V_atalanta", Sample)) ,mapping = aes(y=TE$X2, x=X.masked, label=X.masked)) +
  geom_col(position='stack') +
  theme_classic() +
  scale_fill_manual(values=col_vector) +  xlab('% masked') +
  scale_y_discrete(limits=rev) + ylab('Superfamily') +
  scale_x_break(c(0.1, 0.8)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())

ggplot(data=subset(datafull, grepl("TIR", Superfamily) & grepl("V_tameamea", Sample)) ,mapping = aes(y=TE$X2, x=X.masked, label=X.masked)) +
  geom_col(position='stack') +
  theme_classic() +
  scale_fill_manual(values=col_vector) +  xlab('% masked') +
  scale_y_discrete(limits=rev) + ylab('Superfamily') +
  scale_x_break(c(0.1, 1)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())



# col_vector=col_vector[10:30]
# ggplot(data=subset(datafull, grepl("TIR", Superfamily)),mapping = aes(x=Sample, y=X.masked, fill=Superfamily, label=X.masked)) +
#   geom_col(position='stack') +
#   theme_bw() +
#   scale_fill_manual(values=col_vector) +
#   geom_text(check_overlap = F, position = position_stack(vjust = 0.5)) +
#   ylab('% masked') 

ggplot(data=datafull,mapping = aes(y=X.masked, x=TE$X1, fill=TE$X1)) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(values=col_vector) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(rows = vars(Sample))  +
  ylab('% masked') 

#theme(legend.position = "none")


