setwd("/Users/ssushko/Documents/PROJECTS/3_IBB_vanessa/New_repeat_annotations/commonlib4sp_withnewtameameagenome/")

#V_cardui
landvcard=read.table('summaries/v_card_div_rtype.txt', sep='\t', header=T)

#V_atalanta
landvatal=read.csv('summaries/v_atal_div_rtype.txt', sep='\t', header=T)

#A_io
landaio=read.table('summaries/a_io_div_rtype.txt', sep='\t', header=T)

#A_io
landvtam=read.table('summaries/v_tam_div_rtype.txt', sep='\t', header=T)

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

#colors2=c('#4472c4','#ed7d31','#a5a5a5','#ffc000','#5b9bd5','#70ad47','#997300','#264478')
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
col_vector=col_vector[7:60]

# Separating the columns by delimiter
landvcard= landvcard %>% separate(Type, c('Class', 'Superfam','Order'), sep='/')
landvatal= landvatal %>% separate(Type, c('Class', 'Superfam','Order'), sep='/')
landaio= landaio %>% separate(Type, c('Class', 'Superfam','Order'), sep='/')
landvtam= landvtam %>% separate(Type, c('Class', 'Superfam','Order'), sep='/')

# landvcard$Type = factor(landvcard$Type, levels=c('CLASSI','CRYPTON','DIRS','DNA','HELITRON','LINE','LTR','MAVERICK','MITE','SINE'))

## Cardui ##

ggplot(data=subset(landvcard, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU'), aes(x=Div, y=bp, fill=Superfam)) +
  geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) + ylim(0,16000000) #+ guides(fill="none")

#ggplot(data=landvcard, aes(x=Div, y=bp, fill=Type)) +
#  geom_col(width=1) + theme_classic() + labs(fill='Superfam')

# ggplot(data=subset(landvcard, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU' & Superfam =='LINE'), aes(x=Div, y=bp, fill=Order)) +
#   geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) #+ ylim(0,16000000)
# 
# ggplot(data=subset(landvcard, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU' & Superfam =='LTR'), aes(x=Div, y=bp, fill=Order)) +
#   geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) #+ ylim(0,16000000)
# 
# ggplot(data=subset(landvcard, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU' & Superfam =='TIR'), aes(x=Div, y=bp, fill=Order)) +
#   geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) #+ ylim(0,16000000)


## Atalanta ###

ggplot(data=subset(landvatal, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU'), aes(x=Div, y=bp, fill=Superfam)) +
  geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) + ylim(0,16000000) + guides(fill="none")

# ggplot(data=subset(landvatal, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU' & Superfam =='LINE'), aes(x=Div, y=bp, fill=Order)) +
#   geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) #+ ylim(0,16000000)
# 
# ggplot(data=subset(landvatal, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU' & Superfam =='LTR'), aes(x=Div, y=bp, fill=Order)) +
#   geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) #+ ylim(0,16000000)
# 
# ggplot(data=subset(landvatal, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU' & Superfam =='TIR'), aes(x=Div, y=bp, fill=Order)) +
#   geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) #+ ylim(0,16000000)


### Aio ###

ggplot(data=subset(landaio, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU'), aes(x=Div, y=bp, fill=Superfam)) +
  geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) + ylim(0,16000000) + guides(fill="none")

# ggplot(data=subset(landaio, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU' & Superfam =='LINE'), aes(x=Div, y=bp, fill=Order)) +
#   geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) #+ ylim(0,16000000)
# 
# ggplot(data=subset(landaio, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU' & Superfam =='LTR'), aes(x=Div, y=bp, fill=Order)) +
#   geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) #+ ylim(0,16000000)
# 
# ggplot(data=subset(landaio, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU' & Superfam =='TIR'), aes(x=Div, y=bp, fill=Order)) +
#   geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) #+ ylim(0,16000000)

### Vtame ###

ggplot(data=subset(landvtam, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU'), aes(x=Div, y=bp, fill=Superfam)) +
  geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) + ylim(0,16000000) + guides(fill="none")

# ggplot(data=subset(landvtam, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU' & Superfam =='LINE'), aes(x=Div, y=bp, fill=Order)) +
#   geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) #+ ylim(0,16000000)
# 
# ggplot(data=subset(landvtam, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU' & Superfam =='LTR'), aes(x=Div, y=bp, fill=Order)) +
#   geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) #+ ylim(0,16000000)
# 
# ggplot(data=subset(landvtam, !is.na(Superfam) & Superfam !='Dada' & Superfam !='IS3EU' & Superfam =='TIR'), aes(x=Div, y=bp, fill=Order)) +
  # geom_col(width=1) + theme_classic() + labs(fill='Superfam') + scale_fill_manual(values=col_vector) #+ ylim(0,16000000)

