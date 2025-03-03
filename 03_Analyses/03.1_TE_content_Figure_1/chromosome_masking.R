setwd("/Users/ssushko/Documents/PROJECTS/3_IBB_vanessa/New_repeat_annotations/commonlib4sp_withnewtameameagenome/")

library(ggplot2)
library(tidyr)
library(dplyr)
library(forcats)
library(ggpubr)


vcardchrs=read.csv('./summaries/v_card_mch_perseq.sum', sep=',')
vcardchrs_sorted=vcardchrs[order(vcardchrs$Seq),]
vcardseqreport=read.table('../mchelper_common4species/genome_seq_reports_ncbi/sequence_report_vcardui.tsv', sep='\t', header=T)
vcardseqreport=vcardseqreport[1:36,]
vcardseqreport_sorted=vcardseqreport[order(vcardseqreport$GenBank.seq.accession),]
vcarddf=cbind(vcardchrs_sorted, vcardseqreport_sorted)
vcarddf_=vcarddf[!grepl("^CAJM.+", vcarddf$Seq),]
vcarddf_$Chromosome.name=factor(vcarddf_$Chromosome.name, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,'Z','W'))

vatalchrs=read.csv('./summaries/v_atal_mch_perseq.sum', sep=',')
vatalchrs_sorted=vatalchrs[order(vatalchrs$Seq),]
vatalseqreport=read.table('../mchelper_common4species/genome_seq_reports_ncbi/sequence_report_vatalanta.tsv', sep='\t', header=T)
vatalseqreport_sorted=vatalseqreport[order(vatalseqreport$RefSeq.seq.accession),]
vataldf=cbind(vatalchrs_sorted, vatalseqreport_sorted)
vataldf_=vataldf[!grepl("^NW_+", vataldf$Seq),]
#vataldf_$Chromosome.name=factor(vataldf_$Chromosome.name, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,'Z','W'))
vataldf_$Chromosome.name=factor(vataldf_$Chromosome.name, levels=c(1,2,10,4,8,5,3,9,7,20,6,13,12,11,14,15,16,17,21,18,19,22,23,24,25,26,27,28,29,30,'Z','W'))


aiochrs=read.csv('./summaries/a_io_mch_perseq.sum', sep=',')
aiochrs_sorted=aiochrs[order(aiochrs$Seq),]
aioseqreport=read.table('../mchelper_common4species/genome_seq_reports_ncbi/sequence_report_aglais_io.tsv', sep='\t', header=T)
aiodf=cbind(aiochrs_sorted, aioseqreport)
aiodf_=aiodf[!grepl("^NW_+", aiodf$Seq),]
# aiodf_$Chromosome.name=factor(aiodf_$Chromosome.name, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,'Z'))
aiodf_$Chromosome.name=factor(aiodf_$Chromosome.name, levels=c(1,2,12,4,9,5,3,10,6,18,8,15,14,13,7,16,11,17,21,19,20,23,25,22,24,26,28,29,27,30,'Z'))

vtamechrs=read.csv('./summaries/v_tam_mch_perseq.sum', sep=',')
vtamechrs_sorted=vtamechrs[order(vtamechrs$Seq),]
vtameseqreport=read.table('../mchelper_common4species/genome_seq_reports_ncbi/sequence_report_vtameamea.tsv', sep='\t', header=T)
vtamechrs_sorted=vtamechrs_sorted[!grepl("^JAXCL+", vtamechrs_sorted$Seq),]
vtameseqreport=vtameseqreport[!grepl("^JAXCL+", vtameseqreport$GenBank.seq.accession),]
vtamdf=cbind(vtamechrs_sorted, vtameseqreport)
# vtamdf$Chromosome.name=factor(vtamdf$Chromosome.name, levels=c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,'Z','W'))
vtamdf$Chromosome.name=factor(vtamdf$Chromosome.name, levels=c(2,5,11,3,9,6,4,10,7,21,8,13,14,12,15,16,17,18,22,19,20,23,24,25,26,27,28,29,30,31,'Z','W'))

# bpmasked order
ggplot(vcarddf_) + 
  geom_bar(aes(x=fct_reorder(Chromosome.name, bpMasked), y = Seq.length),stat='identity', fill='white', color='black')+
  geom_col(aes(x=fct_reorder(Chromosome.name, bpMasked), y = bpMasked,fill= Count),position="dodge") +
  scale_fill_gradient(limits=c(7000,27000)) +
  theme_classic() + ylab('Length (bp)') + xlab('Chromosome') + labs(fill='Number of copies')

# chromosome number order
vcard=ggplot(vcarddf_) + 
  geom_bar(aes(x=Chromosome.name, y = Seq.length),stat='identity', fill='white', color='black')+
  geom_col(aes(x=Chromosome.name, y = bpMasked,fill= Count),position="dodge") +
  scale_fill_gradient(limits=c(7000,27000)) +
  theme_classic() + ylab('Length (bp)') + xlab('Chromosome') + labs(fill='Number of copies')

# bpmasked order
ggplot(vataldf_) + 
  geom_bar(aes(x=fct_reorder(Chromosome.name, bpMasked), y = Seq.length),stat='identity', fill='white', color='black')+
  geom_col(aes(x=fct_reorder(Chromosome.name, bpMasked), y = bpMasked,fill= Count),position="dodge") +
  scale_fill_gradient(limits=c(3000,27000)) +
  theme_classic() + ylab('Length (bp)') + xlab('Chromosome') + labs(fill='Number of copies')

# chromosome number order
vatal=ggplot(vataldf_) + 
  geom_bar(aes(x=Chromosome.name, y = Seq.length),stat='identity', fill='white', color='black')+
  geom_col(aes(x=Chromosome.name, y = bpMasked,fill= Count),position="dodge") +
  scale_fill_gradient(limits=c(3000,27000)) +
  theme_classic() + ylab('Length (bp)') + xlab('Chromosome') + labs(fill='Number of copies')

# bpmasked order
ggplot(aiodf_) + 
  geom_bar(aes(x=fct_reorder(Chromosome.name, bpMasked), y = Seq.length),stat='identity', fill='white', color='black')+
  geom_col(aes(x=fct_reorder(Chromosome.name, bpMasked), y = bpMasked,fill= Count),position="dodge") +
  scale_fill_gradient(limits=c(3000,27000)) +
  theme_classic() + ylab('Length (bp)') + xlab('Chromosome') + labs(fill='Number of copies')

# chromosome number order
vtam=ggplot(vtamdf) + 
  geom_bar(aes(x=Chromosome.name, y = Seq.length),stat='identity', fill='white', color='black')+
  geom_col(aes(x=Chromosome.name, y = bpMasked,fill= Count),position="dodge") +
  scale_fill_gradient(limits=c(3000,27000)) +
  theme_classic() + ylab('Length (bp)') + xlab('Chromosome') + labs(fill='Number of copies')

# bpmasked order
ggplot(vtamdf) + 
  geom_bar(aes(x=fct_reorder(Chromosome.name, bpMasked), y = Seq.length),stat='identity', fill='white', color='black')+
  geom_col(aes(x=fct_reorder(Chromosome.name, bpMasked), y = bpMasked,fill= Count),position="dodge") +
  scale_fill_gradient(limits=c(3000,27000)) +
  theme_classic() + ylab('Length (bp)') + xlab('Chromosome') + labs(fill='Number of copies')

# chromosome number order
aio=ggplot(aiodf_) + 
  geom_bar(aes(x=Chromosome.name, y = Seq.length),stat='identity', fill='white', color='black')+
  geom_col(aes(x=Chromosome.name, y = bpMasked,fill= Count),position="dodge") +
  scale_fill_gradient(limits=c(3000,27000)) +
  theme_classic() + ylab('Length (bp)') + xlab('Chromosome') + labs(fill='Number of copies')


ggarrange(vtam,vatal,vcard,aio, ncol=1)

## Linear models to test correlation with chromosome size ##
model=lm(bpMasked ~ Seq.length, data=vcarddf_)
summary(model)
predict <- cbind(vcarddf_, predict(model, interval = 'confidence'))

ggplot(data=predict, aes(x=Seq.length, y=bpMasked)) +
  geom_point() +
  geom_line(aes(Seq.length,fit)) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  geom_text(label=predict$Chromosome.name,
            nudge_x = 0.25, nudge_y = 500000, 
            check_overlap = T)


model=lm(bpMasked ~ Seq.length, data=vataldf_)
summary(model)
predict <- cbind(vataldf_, predict(model, interval = 'confidence'))

ggplot(data=predict, aes(x=Seq.length, y=bpMasked)) +
  geom_point() +
  geom_line(aes(Seq.length,fit)) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  geom_text(label=predict$Chromosome.name,
            nudge_x = 0.25, nudge_y = 500000, 
            check_overlap = T)       

model=lm(bpMasked ~ Seq.length, data=aiodf_)
summary(model)
predict <- cbind(aiodf_, predict(model, interval = 'confidence'))

ggplot(data=predict, aes(x=Seq.length, y=bpMasked)) +
  geom_point() +
  geom_line(aes(Seq.length,fit)) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  geom_text(label=predict$Chromosome.name,
            nudge_x = 0.25, nudge_y = 500000, 
            check_overlap = T) 

model=lm(bpMasked ~ Seq.length, data=vtamdf)
summary(model)
predict <- cbind(vtamdf, predict(model, interval = 'confidence'))

ggplot(data=predict, aes(x=Seq.length, y=bpMasked)) +
  geom_point() +
  geom_line(aes(Seq.length,fit)) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  geom_text(label=predict$Chromosome.name,
            nudge_x = 0.25, nudge_y = 500000, 
            check_overlap = T) 







# Linear regression TE cov vs TE copy number
vcarddf_noW=vcarddf_[!grepl("W", vcarddf_$Chromosome.name),]
model=lm(bpMasked ~ Count, data=vcarddf_noW)
summary(model)
predict <- cbind(vcarddf_noW, predict(model, interval = 'confidence'))
ggplot(data=predict, aes(x=Count, y=bpMasked)) +
  geom_point() +
  geom_line(aes(Count,fit)) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  geom_text(label=predict$Chromosome.name,
            nudge_x = 0.25, nudge_y = 500000, 
            check_overlap = T)

vataldf_noW=vataldf_[!grepl("W", vataldf_$Chromosome.name),]
model=lm(bpMasked ~ Count, data=vataldf_noW)
summary(model)
predict <- cbind(vataldf_noW, predict(model, interval = 'confidence'))
ggplot(data=predict, aes(x=Count, y=bpMasked)) +
  geom_point() +
  geom_line(aes(Count,fit)) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  geom_text(label=predict$Chromosome.name,
            nudge_x = 0.25, nudge_y = 500000, 
            check_overlap = T)

model=lm(bpMasked ~ Count, data=aiodf_)
summary(model)
predict <- cbind(aiodf_, predict(model, interval = 'confidence'))
ggplot(data=predict, aes(x=Count, y=bpMasked)) +
  geom_point() +
  geom_line(aes(Count,fit)) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  geom_text(label=predict$Chromosome.name,
            nudge_x = 0.25, nudge_y = 500000, 
            check_overlap = T)

vtamdf_noW=vtamdf[!grepl("W", vtamdf$Chromosome.name),]
model=lm(bpMasked ~ Count, data=vtamdf_noW)
summary(model)
predict <- cbind(vtamdf_noW, predict(model, interval = 'confidence'))
ggplot(data=predict, aes(x=Count, y=bpMasked)) +
  geom_point() +
  geom_line(aes(Count,fit)) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  geom_text(label=predict$Chromosome.name,
            nudge_x = 0.25, nudge_y = 500000, 
            check_overlap = T)
