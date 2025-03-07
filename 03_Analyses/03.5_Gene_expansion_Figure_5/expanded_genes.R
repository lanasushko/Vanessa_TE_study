setwd("/tmp/global2/ssushko/Vanessa_proj/6_gene_expansion/expanded_genes_cardui")

#### INSIDE GENES ####

## Genome-wide
expanded=read.table('expanded_genes/expanded_average_per_gene.txt', header=T)
genomic.av=read.table('expanded_genes/genome-wide/genomic_average_per_gene.txt', header=T)

# normalize by length -- insertions per gene bp
expanded$ins.norm=expanded$insertions/expanded$length
genomic.av$ins.norm=genomic.av$insertions/genomic.av$length

mean.exp=mean(expanded$ins.norm)
mean.gen.av=mean(genomic.av$ins.norm)

expanded$type='expanded'
genomic.av$type='genomic-av'

df=rbind(expanded,genomic.av)

library(ggplot2)
library(ggstatsplot)

ggplot(df, mapping=aes(y=ins.norm, x=type)) +
  geom_boxplot()

ggbetweenstats(
  data  = df,
  x     = type,
  y     = ins.norm
)

t.test(expanded$ins.norm, genomic.av$ins.norm,
       alternative = c("greater"),
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95)

## Random sample

sample1=read.table('expanded_genes/random_sample/genes_sample_77_rs102.gff.insertions.table', header=T)
sample2=read.table('expanded_genes/random_sample/genes_sample_77_rs12.gff.insertions.table', header=T)
sample3=read.table('expanded_genes/random_sample/genes_sample_77_rs42.gff.insertions.table', header=T)
sample4=read.table('expanded_genes/random_sample/genes_sample_77_rs53.gff.insertions.table', header=T)
sample5=read.table('expanded_genes/random_sample/genes_sample_77_rs98.gff.insertions.table', header=T)

sample1$ins.norm=sample1$insertions/sample1$length
sample2$ins.norm=sample2$insertions/sample2$length
sample3$ins.norm=sample3$insertions/sample3$length
sample4$ins.norm=sample4$insertions/sample4$length
sample5$ins.norm=sample5$insertions/sample5$length


sample1$type='random_sample'
sample2$type='random_sample'
sample3$type='random_sample'
sample4$type='random_sample'
sample5$type='random_sample'

df1=rbind(expanded,sample1)
df2=rbind(expanded,sample2)
df3=rbind(expanded,sample3)
df4=rbind(expanded,sample4)
df5=rbind(expanded,sample5)

ggbetweenstats(
  data  = df1,
  x     = type,
  y     = ins.norm
)

ggbetweenstats(
  data  = df2,
  x     = type,
  y     = ins.norm
)

ggbetweenstats(
  data  = df3,
  x     = type,
  y     = ins.norm
)

ggbetweenstats(
  data  = df4,
  x     = type,
  y     = ins.norm
)

ggbetweenstats(
  data  = df5,
  x     = type,
  y     = ins.norm
)


## Matched control genes




##### GENE FLANKS #####

## Genome-wide
expanded.FL=read.table('flanks/expanded_ins_pergeneflanks.txt')
genomic.av.FL=read.table('flanks/genome-wide/genomic_ins_pergeneflanks.txt')
colnames(expanded.FL) = c('insertions', 'gene')
colnames(genomic.av.FL) = c('insertions', 'gene')

mean.exp.fl=mean(expanded.FL$insertions)
mean.gen.av.fl=mean(genomic.av.FL$insertions)

expanded.FL$type='expanded'
genomic.av.FL$type='genomic-av'

df=rbind(expanded.FL,genomic.av.FL)

library(ggplot2)

ggplot(df, mapping=aes(y=insertions, x=type)) +
  geom_boxplot()

ggbetweenstats(
  data  = df,
  x     = type,
  y     = insertions
)


t.test(expanded.FL$insertions, genomic.av.FL$insertions,
       alternative = c("greater"),
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95)


## Random sample

fl_sample1=read.table('flanks/random_sample/77_random_flanks_rs42.table', header=T)
fl_sample2=read.table('flanks/random_sample/77_random_flanks_rs102.table', header=T)
fl_sample3=read.table('flanks/random_sample/77_random_flanks_rs12.table', header=T)
fl_sample4=read.table('flanks/random_sample/77_random_flanks_rs53.table', header=T)
fl_sample5=read.table('flanks/random_sample/77_random_flanks_rs98.table', header=T)

fl_sample1$type='random_sample'
fl_sample2$type='random_sample'
fl_sample3$type='random_sample'
fl_sample4$type='random_sample'
fl_sample5$type='random_sample'

df_fl1=rbind(expanded.FL,fl_sample1)
df_fl2=rbind(expanded.FL,fl_sample2)
df_fl3=rbind(expanded.FL,fl_sample3)
df_fl4=rbind(expanded.FL,fl_sample4)
df_fl5=rbind(expanded.FL,fl_sample5)

ggbetweenstats(
  data  = df_fl1,
  x     = type,
  y     = insertions
)

ggbetweenstats(
  data  = df_fl2,
  x     = type,
  y     = insertions
)

ggbetweenstats(
  data  = df_fl3,
  x     = type,
  y     = insertions
)

ggbetweenstats(
  data  = df_fl4,
  x     = type,
  y     = insertions
)

ggbetweenstats(
  data  = df_fl5,
  x     = type,
  y     = insertions
)


## Matched control gene flanks

