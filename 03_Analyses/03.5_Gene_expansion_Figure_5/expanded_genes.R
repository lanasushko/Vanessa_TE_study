setwd("~/Documents/PROJECTS/3_IBB_vanessa/enrichment/expanded_genes")

#### INSIDE GENES ####

expanded=read.table('expanded_average_per_gene.txt', header=T)
genomic.av=read.table('genomic_average_per_gene.txt', header=T)

expanded$ins.norm=expanded$insertions/expanded$length
genomic.av$ins.norm=genomic.av$insertions/genomic.av$length

mean.exp=mean(expanded$ins.norm)
mean.gen.av=mean(genomic.av$ins.norm)

expanded$type='expanded'
genomic.av$type='genomic-av'

df=rbind(expanded,genomic.av)

library(ggplot2)

ggplot(df, mapping=aes(y=ins.norm, x=type)) +
  geom_boxplot()


t.test(expanded$ins.norm, genomic.av$ins.norm,
       alternative = c("greater"),
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95)


##### GENE FLANKS #####

expanded.FL=read.table('expanded_ins_pergeneflanks.txt', header=T)
genomic.av.FL=read.table('genomic_ins_pergeneflanks.txt', header=T)

mean.exp.fl=mean(expanded.FL$insertions)
mean.gen.av.fl=mean(genomic.av.FL$insertions)

expanded.FL$type='expanded'
genomic.av.FL$type='genomic-av'

df=rbind(expanded.FL,genomic.av.FL)

library(ggplot2)

ggplot(df, mapping=aes(y=insertions, x=type)) +
  geom_boxplot()


t.test(expanded.FL$insertions, genomic.av.FL$insertions,
       alternative = c("greater"),
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95)
