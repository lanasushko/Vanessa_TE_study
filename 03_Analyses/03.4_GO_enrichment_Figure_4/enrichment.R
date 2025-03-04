setwd("~/Documents/PROJECTS/3_IBB_vanessa/manuscript/FIG4_enrichment/")

library(clusterProfiler)

### CARDUI ###

seeds=read.table('cardui.seeds_with_MITEs.list')
# seeds=read.table('card_seeds_with_flankMITEs.list')
seeds=seeds[['V1']]
card.seeds=seeds
seeduniverse=read.table('cardui.seeds.list')
seeduniverse=seeduniverse[['V1']]
card.seeduniverse=seeduniverse
TERM2seed=read.table('cardui.GO-seed.2cols.list')
card.TERM2seed=TERM2seed

### ATALANTA ###
seeds=read.table('atalanta.seeds_with_MITEs.list')
# seeds=read.table('atal_seeds_with_flankMITEs.list')
seeds=seeds[['V1']]
atal.seeds=seeds
seeduniverse=read.table('atalanta.seeds.list')
seeduniverse=seeduniverse[['V1']]
atal.seeduniverse=seeduniverse
TERM2seed=read.table('atalanta.GO-seed.2cols.list')
atal.TERM2seed=TERM2seed

### AGLAIS IO ###
seeds=read.table('aglaisio.seeds_with_MITEs.list')
# seeds=read.table('aio_seeds_with_flankMITEs.list')
seeds=seeds[['V1']]
aio.seeds=seeds
seeduniverse=read.table('aglaisio.seeds.list')
seeduniverse=seeduniverse[['V1']]
aio.seeduniverse=seeduniverse
TERM2seed=read.table('aglaisio.GO-seed.2cols.list')
aio.TERM2seed=TERM2seed

### TAMEAMEA ###
seeds=read.table('tameamea.seeds_with_MITEs.list')
# seeds=read.table('tam_seeds_with_flankMITEs.list')
seeds=seeds[['V1']]
tam.seeds=seeds
seeduniverse=read.table('tameamea.seeds.list')
seeduniverse=seeduniverse[['V1']]
tam.seeduniverse=seeduniverse
TERM2seed=read.table('tameamea.GO-seed.2cols.list')
tam.TERM2seed=TERM2seed

library(GO.db)
# extract a named vector of all terms
goterms <- Term(GOTERM)
#convert into a data frame
TERM2NAME <- data.frame("GOID"=names(goterms),"term"=goterms )

card.enrichseed=enricher(gene=card.seeds, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                universe=card.seeduniverse, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, 
                TERM2GENE=card.TERM2seed, TERM2NAME=TERM2NAME)
atal.enrichseed=enricher(gene=atal.seeds, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                         universe=atal.seeduniverse, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, 
                         TERM2GENE=atal.TERM2seed, TERM2NAME=TERM2NAME)
tam.enrichseed=enricher(gene=tam.seeds, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                         universe=tam.seeduniverse, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, 
                         TERM2GENE=tam.TERM2seed, TERM2NAME=TERM2NAME)
aio.enrichseed=enricher(gene=aio.seeds, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                         universe=aio.seeduniverse, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, 
                         TERM2GENE=aio.TERM2seed, TERM2NAME=TERM2NAME)

library(DOSE)
library(pheatmap)
library(enrichplot)
library(ggupset)
library(ggplot2)

# Barplot
barplot(card.enrichseed, showCategory = 20) 
mutate(enrichseed, qscore = -log(p.adjust, base = 10)) %>% 
  barplot(x = "qscore")

# Dotplot
dotplot(enrichseed, showCategory = 15)

# Cnetplot
cnetplot(enrichseed)

# Heatplot
heatplot(enrichseed, showCategory = 5)

# Treeplot
card.enrichres <- pairwise_termsim(card.enrichseed) 
atal.enrichres <- pairwise_termsim(atal.enrichseed) 
tam.enrichres <- pairwise_termsim(tam.enrichseed) 
aio.enrichres <- pairwise_termsim(aio.enrichseed) 


treeplot(card.enrichres) #+ guides(size = guide_legend(override.aes = list(size = c(100,200,300))))
treeplot(atal.enrichres)
treeplot(tam.enrichres)
treeplot(aio.enrichres)


























# random sample test

foe=sample(universe, 3897)
foenrich=enricher(gene=foe, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                universe=universe, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, 
                TERM2GENE=TERM2GENE, TERM2NAME=TERM2NAME)

barplot(foenrich, showCategory = 20) 
# ok this shows 


##################################
# for 4 species comparative plot #
##################################

########### combine enricher results ############
df_cardui=card.enrichseed@result
df_atalanta=atal.enrichseed@result
df_tameamea=tam.enrichseed@result
df_aio=aio.enrichseed@result
# Combine the data frames
combined_df <- rbind(df_cardui, df_atalanta,df_tameamea,df_aio)

#
compare_result <- new("compareClusterResult",
                      compareClusterResult = combined_df,
                      geneClusters = list(
                        cardui = card.seeds,
                        atalanta = atal.seeds,
                        tameamea = tam.seeds,
                        aio = aio.seeds
                      ),
                      fun = "enricher",
                      readable=FALSE ,
                      keytype='UNKNOWN'
)

compare <- pairwise_termsim(compare_result)

treeplot(compare_result)

df_cardui=na.omit(card.enrichres@result)
df_cardui$species=sample('cardui',nrow(df_cardui),replace=T)
df_atalanta=na.omit(atal.enrichres@result)
df_atalanta$species=sample('atal',nrow(df_atalanta),replace=T)
df_tameamea=na.omit(tam.enrichres@result)
df_tameamea$species=sample('tameamea',nrow(df_tameamea),replace=T)
df_aio=na.omit(aio.enrichres@result)
df_aio$species=sample('aio',nrow(df_aio),replace=T)

combined_df <- rbind(df_cardui,df_atalanta,df_tameamea,df_aio)


meanGOpvalues=data.frame(matrix(ncol = 2))
colnames(meanGOpvalues)=c('GO','mean.p.adj')

GOsin4=subset(as.data.frame(table(combined_df$ID)),as.data.frame(table(combined_df$ID))[,2]==4)
allGOin4=as.vector(unique(GOsin4$Var1))


for (go in 1:length(allGOin4)) {
  mean=mean(c(df_cardui[allGOin4[go],]$p.adjust,df_atalanta[allGOin4[go],]$p.adjust,df_tameamea[allGOin4[go],]$p.adjust,df_aio[allGOin4[go],]$p.adjust),na.rm = TRUE)
  output=c(allGOin4[go], mean)
  meanGOpvalues = rbind(meanGOpvalues, output)
}

meanGOpvalues$mean.p.adj=as.numeric(as.character(meanGOpvalues$mean.p.adj))
top20=meanGOpvalues[order(as.numeric(as.character(meanGOpvalues$mean.p.adj)),decreasing = F),][c(1:20),]
top30=meanGOpvalues[order(as.numeric(as.character(meanGOpvalues$mean.p.adj)),decreasing = F),][c(1:30),]

combined_df_top20=combined_df[combined_df$ID %in% top20$GO, ]
combined_df_top30=combined_df[combined_df$ID %in% top30$GO, ]

ggplot(combined_df_top30, aes(y=Description,x=species, fill=p.adjust)) +
  geom_tile() + scale_fill_continuous(trans = 'reverse')

ggplot(combined_df_top30, aes(y=Description,x=species, color=p.adjust, size=Count)) +
  geom_point() + scale_color_continuous(trans = 'reverse') + theme_classic()



########### COMPARECLUSTER ############
term2seed=unique(rbind(card.TERM2seed,atal.TERM2seed,tam.TERM2seed,aio.TERM2seed))

gene_clusters <- list(
  cardui = card.seeds,
  atalanta = atal.seeds,
  tameamea = tam.seeds,
  aio = aio.seeds
)

compare_result_compcl <- compareCluster(geneClusters = gene_clusters, fun = "enricher", TERM2GENE = term2seed, TERM2NAME=TERM2NAME)

compare <- pairwise_termsim(compare_result)
treeplot(compare, hclust_method = "average", showCategory=20)

