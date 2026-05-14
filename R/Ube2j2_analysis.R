# Analysis of ER stress-induced changes in protein abundance associated with Ube2j2 deficiency in cell culture

#### WGCNA ####
require(WGCNA); require(flashClust)

options(stringsAsFactors = FALSE);

# Allow up to 4 processes
#enableWGCNAThreads(2);

data<-read.table("Data/processed_minprob.tsv", sep="\t", header=T, row.names=1)

datExp0<-data[,-c(1:2)]
row.names(datExp0)<-data[,2]

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=40, by=4))

# Call the network topology analysis function
sft = pickSoftThreshold(datExp0, powerVector = powers, verbose = 5)

# Plot the results
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red");

# This line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.91,col="red")

softPower = 10;

adjacency = adjacency(datExp0, power = softPower, type = "unsigned");

# Turn adjacency into topological overlap
dissTOM = 1-TOMsimilarity(adjacency, TOMType = "unsigned")

# Remove large objects and clear environment memory
rm(adjacency)
collectGarbage()

# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# set the minimum module size relatively high to make bigger modules:
minModuleSize = 100;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)

print(table(dynamicColors))

# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExp0, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Choose a height cut of 0.20, corresponding to correlation of 0.80
MEDissThres = 0.20

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(datExp0, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors;

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

#pdf(file = "geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Save module colors and labels for use in subsequent parts
# You can skip this. It takes up space and isn't necessary.
#save(MEs, moduleLabels, moduleColors, geneTree, file = "networkConstruction.RData")

# Save gene - module assignments
gene_modules<-data.frame(gene=colnames(datExp0),module=moduleColors, uniprot=sub("(.+)\\.{1}\\d+$", "\\1", colnames(datExp0), perl=T))
gene_modules$uniprot<-sub(".+\\.(\\w+)$", "\\1", gene_modules$uniprot, perl=T)

# Relating modules to experimental factors
sample_table<-read.csv("sample_info.csv")

row.names(sample_table)<-sample_table$sample

# Get the condition info and sort based on library names
conditions<-sample_table[row.names(datExp0),]

conditions_ME<-cbind(conditions,MEs)

head(conditions_ME)

require(plyr); require(ggplot2); require(gridExtra); require(ggdendro)

# Calculate distance and cluster
dissimME=(1-t(cor(MEs, method="pearson")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )

# Plot the eigengene dendrogram
plot(hclustdatME, main="Clustering tree based of the module eigengenes")

# Plot a simple dendrogram with ggdendro package
plot1 <- ggdendrogram(hclustdatME, rotate = T, labels=F)+ theme(axis.text.y=element_blank())
plot1

#### ME trait analysis ####
# Create an empty variable that we can use to store the summary statistics
require(emmeans)
summarize_ME<-NULL
pairwise_ME=NULL
ME.test<-NULL

for(j in 9:(ncol(conditions_ME)-1))
{
  # Summarize each eigengene (EG) based on condition  
  EG<-conditions_ME[,c(7,j)]
  names(EG)[2]<-"EG"
  tmp.avg<-ddply(EG, .(condition), summarise, tmp.avg=mean(EG), tmp.se=sd(EG)/sqrt(length(EG)))
  
  ME.t<-anova(lm(conditions_ME[,j]~conditions_ME$genotype + conditions_ME$drug + conditions_ME$time))
  
  model<-lm(conditions_ME[,j]~conditions_ME$condition)
  
  em2<-as.data.frame(pairs(emmeans(model,  ~ condition), adjust="none"))
  
  pairwise_ME<-rbind(pairwise_ME,c(em2[17, 6], # "C6 Tu 18hrs - WT Tu 18hrs"
                                   em2[11, 6], #"C6 no time6 - WT no time6"
                                   em2[1, 6], #"C6 no time6 - C6 no time18"
                                   em2[23, 6])) #"WT no time6 - WT no time18"
  
  
  # Append the results to the bottom of the table 
  summarize_ME<-rbind(summarize_ME,cbind(names(conditions_ME)[j],tmp.avg))
  ME.test<-rbind(ME.test, c(names(conditions_ME)[j],ME.t$`Pr(>F)`[1:3]))
}

#output main effects p-value
names(summarize_ME)[1]<-"module"
colnames(ME.test)<-c("module", "genotype_pvalue", "drug_pvalue", "time_pvalue")
write.table(ME.test, file="Data/Module_main_effects.txt", sep="\t", quote=F, row.names = F)

#output ME pairwise p-values and adjusted p-values
colnames(pairwise_ME)<-c("C6_Tu_time18-WT_Tu_time18",
                         "C6_no_time6-WT_no_time6",
                         "C6_no_time6-C6_no_time18",
                         "WT_no_time6-WT_no_time18")
pairwise_ME_out<-as.data.frame(cbind(ME.test$module,pairwise_ME))
names(pairwise_ME_out)[1]<-"module"
write.table(pairwise_ME_out, file="Data/Module_pairwise_pvalue.txt",sep="\t", quote=F,row.names=F)


pairwise_ME_padj<-pairwise_ME_out
pairwise_ME_padj[,2]<-p.adjust(pairwise_ME_padj[,2])
pairwise_ME_padj[,3]<-p.adjust(pairwise_ME_padj[,3])
pairwise_ME_padj[,4]<-p.adjust(pairwise_ME_padj[,4])
pairwise_ME_padj[,5]<-p.adjust(pairwise_ME_padj[,5])
write.table(pairwise_ME_padj, file="Data/Module_pairwise_adjusted_pvalue.txt",sep="\t", quote=F,row.names=F)

my_palette <- colorRampPalette(c("white", "red"))(n = 20)

pdf(file="Data/WGCNA_trait_heatmap.pdf",width=4,height=8)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = matrix(-log10(as.numeric(ME.test[,2:4])), ncol=3, nrow=12),
               xLabels = colnames(ME.test)[2:4],
               yLabels = ME.test[,1],
               yColorLabels = TRUE,
               yColorWidth = 0.05,
               ySymbols = ME.test[,1],
               colors = my_palette,
               #textMatrix = txt,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.y = 0.8,
               zlim = c(0,32),
               #main = paste("Module-trait relationships")
)
dev.off()

pdf(file="Data/WGCNA_pairwise_heatmap.pdf",width=4,height=8)
my_palette <- colorRampPalette(c("white", "red"))(n = 20)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = -log10(pairwise_ME_padj[,-c(1)]),
               xLabels = colnames(pairwise_ME_padj)[-1],
               yLabels = ME.test[,1],
               yColorLabels = TRUE,
               yColorWidth = 0.05,
               ySymbols = ME.test[,1],
               colors = my_palette,
               #textMatrix = txt,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.y = 0.8,
               zlim = c(0,32),
               #main = paste("Module-trait relationships")
)
dev.off()


#### GO analysis ####
require(simplifyEnrichment); require(clusterProfiler); require(org.Hs.eg.db)

gene_modules$genes<-sub("(\\w+)\\..+","\\1",gene_modules$gene, perl=T)
gene_modules$uniprot<-sub("\\w+\\.(.+)","\\1",gene_modules$gene, perl=T)
gene_modules$uniprot<-sub("(\\w+)\\..+","\\1",gene_modules$uniprot, perl=T)
write.table(gene_modules, file="Data/Protein_modules.txt", sep="\t")

GOresults_new<-list()

# Run GO enrichment on each set of module genes
for(g in unique(gene_modules[,2]))
{
  genes<-gene_modules[which(gene_modules[,2]==g),3]
  
  GOresults_new[[g]]<-       enrichGO(gene          = genes,
                                      OrgDb         = org.Hs.eg.db,
                                      keyType  = 'UNIPROT',
                                      ont           = "BP",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff  = 0.01,
                                      qvalueCutoff  = 0.05,
                                      readable      = TRUE)
  print(g)
}

#save(GOresults_new,file="012725_GO_results_P10.rdata")

GOresults_CC<-list()
GOresults_MF<-list()

# Run GO enrichment on each set of module genes
for(g in unique(gene_modules[,2]))
{
  genes<-gene_modules[which(gene_modules[,2]==g),3]
  
  GOresults_MF[[g]]<-       enrichGO(gene          = genes,
                                     OrgDb         = org.Hs.eg.db,
                                     keyType  = 'UNIPROT',
                                     ont           = "MF",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.01,
                                     qvalueCutoff  = 0.05,
                                     readable      = TRUE)
  
  GOresults_CC[[g]]<-       enrichGO(gene          = genes,
                                     OrgDb         = org.Hs.eg.db,
                                     keyType  = 'UNIPROT',
                                     ont           = "CC",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.01,
                                     qvalueCutoff  = 0.05,
                                     readable      = TRUE)
  print(g)
}

#save(GOresults_CC, GOresults_MF,file="013025_GO_CC_MF_P10.rdata")


require(simplifyEnrichment); library(tm); library(ComplexHeatmap)
GOresults_05 = lapply(GOresults_new, function(x) x@result[which(x@result$p.adjust < 0.05),])

M_order<-sub("ME","",ME_order)

pdf(file="Data/GOsim.pdf",w=10,h=12)
simplifyGOFromMultipleLists(GOresults_05[M_order],padj_cutoff = 0.05,ont = "BP",db = org.Hs.eg.db, measure = "Wang", method = "mclust", heatmap_param = list(breaks=c(5e-02, 5e-04), col = c("transparent", "blue")))
dev.off()


####save GO results to excel
require(openxlsx)
wb <- createWorkbook()

lapply(seq_along(GOresults_05), function(i){
  addWorksheet(wb=wb, sheetName = names(GOresults_05[i]))
  writeData(wb, sheet = i, GOresults_05[[i]][-length(GOresults_05[[i]])])
})

#Save Workbook
#saveWorkbook(wb, "012825_GO_Results_P10.xlsx", overwrite = TRUE)

#MF
wb_MF <- createWorkbook()

lapply(seq_along(GOresults_MF), function(i){
  addWorksheet(wb=wb_MF, sheetName = names(GOresults_MF[i]))
  writeData(wb_MF, sheet = i, GOresults_MF[[i]][-length(GOresults_MF[[i]])])
})

#Save Workbook
#saveWorkbook(wb_MF, "012825_GO_Results_P10_MF.xlsx", overwrite = TRUE)

#CC
wb_CC <- createWorkbook()

lapply(seq_along(GOresults_CC), function(i){
  addWorksheet(wb=wb_CC, sheetName = names(GOresults_CC[i]))
  writeData(wb_CC, sheet = i, GOresults_CC[[i]][-length(GOresults_CC[[i]])])
})

#Save Workbook
#saveWorkbook(wb_CC, "012825_GO_Results_P10_CC.xlsx", overwrite = TRUE)

#### DEG analysis ####

sample_table<-read.csv("Data/sample_info.csv")
row.names(sample_table)<-sample_table$sample
lib_names<-row.names(sample_table)

#build the statistical model
genotype<-as.factor(sample_table[lib_names,"genotype"])
drug<-as.factor(sample_table[lib_names,"drug"])
time<-sample_table[lib_names,"time"]
group<-as.factor(sample_table[lib_names,"group"])

DE_results<-data.frame(gene=NA, c1=NA, c2=NA, c3=NA, Estimate3=NA, c4=NA, c5=NA)

require(emmeans)

for(g in 1:ncol(datExp0))
{
  tmp<-data.frame(gene=datExp0[lib_names,g], genotype=genotype,drug=drug,time=time, group=group)
  
  m1<-glm(gene ~ genotype * drug * time , data=tmp)
  em1<-as.data.frame(pairs(emmeans(m1,  ~ genotype * drug * time), adjust="none"))
  #emmip(m1, drug~time | genotype)
  
  DE_results[g,1]<-names(datExp0)[g]
  DE_results[g,2]<-em1[which(em1$contrast == "C6 Tu time18 - WT Tu time18"), 6]
  DE_results[g,3]<-em1[which(em1$contrast == "C6 no time6 - WT no time6"), 6]
  DE_results[g,4]<-em1[which(em1$contrast == "C6 no time6 - WT Tu time18"), 6]
  DE_results[g,5]<-em1[which(em1$contrast == "C6 no time6 - WT Tu time18"), 2]
  DE_results[g,6]<-em1[which(em1$contrast == "C6 no time6 - C6 no time18"), 6]
  DE_results[g,7]<-em1[which(em1$contrast == "WT no time6 - WT no time18"), 6]
  
}


DE_results[,8]<-p.adjust(DE_results[,2])
DE_results[,9]<-p.adjust(DE_results[,3])
DE_results[,10]<-p.adjust(DE_results[,6])
DE_results[,11]<-p.adjust(DE_results[,7])


names(DE_results)[2:11]<-c("C6_Tu_time18-WT_Tu_time18", "C6_no_time6-WT_no_time6", "C6_no_time6-WT_Tu_time18", "C6_no_time6-WT_Tu_time18_ESTIMATE", "C6_no_time6-C6_no_time18", "WT_no_time6-WT_no_time18", "C6_Tu_time18-WT_Tu_time18_padj", "C6_no_time6-WT_no_time6_padj", "C6_no_time6-C6_no_time18_padj", "WT_no_time6-WT_no_time18_padj")
row.names(DE_results)<-DE_results[,1]
#write.table(DE_results,file="012125_DEP_results.txt", sep="\t")


#### HUB gene analysis for modules####

ModuleMembership<-as.data.frame(cor(datExp0, MEs, use = "p"));

HubGene=list()
HubData=list()

mods=c("blue", "black", "purple", "brown", "midnightblue", "magenta")

for(i in mods)
{
  gm<-gene_modules[which(gene_modules$module==i), "gene"]
  GeneSig1<-as.numeric(DE_results[gm,"C6_no_time6-WT_no_time6_padj"])
  GeneSig2<-as.numeric(DE_results[gm,"C6_Tu_time18-WT_Tu_time18_padj"])
  GeneSig3<-as.numeric(DE_results[gm,"C6_no_time6-WT_Tu_time18"])
  GeneEst3<-as.numeric(DE_results[gm,"C6_no_time6-WT_Tu_time18_ESTIMATE"])
  GeneSig4<-as.numeric(DE_results[gm,"C6_no_time6-C6_no_time18_padj"])
  GeneSig5<-as.numeric(DE_results[gm,"WT_no_time6-WT_no_time18_padj"])
  mg<-ModuleMembership[gm,paste0("ME",i)]
  HubData[[i]]<-data.frame(gene=gm,"C6_no_time6-WT_no_time6_padj"=GeneSig1,"C6_Tu_time18-WT_Tu_time18_padj"=GeneSig2, "C6_no_time6-WT_Tu_time18_pvalue"=GeneSig3, "C6_no_time6-WT_Tu_time18_Est"=GeneEst3, "C6_no_time6-C6_no_time18_padj"=GeneSig4, "WT_no_time6-WT_no_time18_padj"=GeneSig5, MM=mg)
  
  if(i=="black")
  { 
    HubGene[[i]]<-HubData[[i]][which(HubData[[i]]$C6_no_time6.WT_Tu_time18_pvalue>0.7  & abs(HubData[[i]]$MM)>=0.85),]
    HubGene[[i]]<-merge(HubGene[[i]], gene_modules, by.x="gene", by.y="gene", all.x=T)
  } else if ( i=="purple"| i=="brown" | i=="midnightblue")
  {
    HubGene[[i]]<-HubData[[i]][which((HubData[[i]]$C6_no_time6.C6_no_time18_padj<0.01 | HubData[[i]]$WT_no_time6.WT_no_time18_padj<0.01) & abs(HubData[[i]]$MM)>=0.85),]
    HubGene[[i]]<-merge(HubGene[[i]], gene_modules, by.x="gene", by.y="gene", all.x=T)
  } else if ( i=="blue" | i=="magenta")
  {
    HubGene[[i]]<-HubData[[i]][which((HubData[[i]]$C6_no_time6.WT_no_time6<0.01 | HubData[[i]]$C6_Tu_time18.WT_Tu_time18<0.01) & abs(HubData[[i]]$MM)>=0.85),]
    HubGene[[i]]<-merge(HubGene[[i]], gene_modules, by.x="gene", by.y="gene", all.x=T)
  }
  
}

for(u in mods)
{
  
  if(u=="blue"| u=="magenta")
  {
    tp<-reshape2::melt(HubData[[u]],id.vars=c("gene","MM"), meassure.vars=c("C6_no_time6.WT_no_time6_padj", "C6_Tu_time18.WT_Tu_time18_padj"), variable.name="contrast", value.name = "p.value")
    
    ggplot(tp[which(tp$contrast %in% c("C6_no_time6.WT_no_time6_padj", "C6_Tu_time18.WT_Tu_time18_padj")),],aes(x=abs(MM),y=-log10(p.value), group=contrast, colour=contrast)) + geom_point() +
      #facet_grid(.~contrast) + 
      xlab(paste0("Module Membership in ",u)) + geom_vline(aes(xintercept = 0.85), colour="red") + geom_hline(aes(yintercept = 2),colour="red") + theme_bw()
    ggsave(paste0("HubGenes/",u,"_HubGenes.pdf"),h=4,w=8)
    
  } else if (u=="black")
  {
    
    tp<-reshape2::melt(HubData[[u]],id.vars=c("gene","MM"), meassure.vars=c("C6_no_time6-WT_Tu_time18_pvalue"), variable.name="contrast", value.name = "p.value")
    
    ggplot(tp[which(tp$contrast == "C6_no_time6.WT_Tu_time18_pvalue"),],aes(x=abs(MM),y=p.value)) + geom_point() + facet_grid(.~contrast) + xlab(paste0("Module Membership in ",u)) + geom_vline(aes(xintercept = 0.85), colour="red") + geom_hline(aes(yintercept = 0.70),colour="red")
    ggsave(paste0("HubGenes/",u,"_HubGenes.pdf"),h=4,w=8)
    
  } else if (u=="purple"|u=="brown"|u=="midnightblue")
  {
    
    tp<-reshape2::melt(HubData[[u]],id.vars=c("gene","MM"), meassure.vars=c("C6_no_time6-C6_no_time18_padj", "WT_no_time6-WT_no_time18_padj"), variable.name="contrast", value.name = "p.value")
    
    ggplot(tp[which(tp$contrast %in% c("C6_no_time6.C6_no_time18_padj", "WT_no_time6.WT_no_time18_padj")),],aes(x=abs(MM),y=-log10(p.value))) + geom_point() + facet_grid(.~contrast) + xlab(paste0("Module Membership in ",u)) + geom_vline(aes(xintercept = 0.85), colour="red") + geom_hline(aes(yintercept = 2),colour="red")
    ggsave(paste0("HubGenes/",u,"_HubGenes.pdf"),h=4,w=8) 
  }
  
}


require(openxlsx)
wb <- createWorkbook()

for(i in 1:length(HubGene))
{
  addWorksheet(wb, sheetName = names(HubGene)[i], gridLines = FALSE)
  writeDataTable(wb, sheet = i, x = HubGene[[i]])
  
}
#saveWorkbook(wb, "012425_HubGenes.xlsx", overwrite = TRUE) 

gene_modules_DEG<-merge(gene_modules, DE_results, by.x="gene", by.y="gene", all.x=T)
#write.table(gene_modules_DEG, file="012425_GeneModules_FullResults.txt", sep="\t")



#### magenta hubgenes GO ####
require(openxlsx)
magenta_hub<-read.xlsx(xlsxFile = "Data/012425_HubGenes.xlsx",sheet = "magenta")
magenta_hub_genes<-sub("(\\w+)\\..+","\\1",magenta_hub$gene, perl=T)
magenta_hub_uniprot<-sub("\\w+\\.(.+)","\\1",magenta_hub$gene, perl=T)
magenta_hub_uniprot<-sub("(\\w+)\\..+","\\1",magenta_hub_uniprot, perl=T)

magenta_hub_GO<-list()
magenta_hub_GO$BP<-       enrichGO(gene          = magenta_hub_uniprot,
                                OrgDb         = org.Hs.eg.db,
                                keyType  = 'UNIPROT',
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)
magenta_hub_GO$CC<-       enrichGO(gene          = magenta_hub_uniprot,
                                   OrgDb         = org.Hs.eg.db,
                                   keyType  = 'UNIPROT',
                                   ont           = "CC",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
magenta_hub_GO$MF<-       enrichGO(gene          = magenta_hub_uniprot,
                                   OrgDb         = org.Hs.eg.db,
                                   keyType  = 'UNIPROT',
                                   ont           = "MF",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
magenta_hub_GO$ALL<-       enrichGO(gene          = magenta_hub_uniprot,
                                   OrgDb         = org.Hs.eg.db,
                                   keyType  = 'UNIPROT',
                                   ont           = "ALL",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)


View(magenta_hub_GO$ALL@result)
#save(magenta_hub_GO ,file="Magenta_hub_GO_results.rdata")
#write.table(magenta_hub_GO$ALL@result, file="Magenta_HubGenes_GO.txt", sep="\t")

#### Add gene plots ####
require(ggplot2); require(stringr)

# load metadata
sample_table<-read.csv("Data/sample_info.csv")
row.names(sample_table)<-sample_table$sample
lib_names<-row.names(sample_table)

  sample_table$genotype[which(sample_table$genotype=="C6")]="Ube2j2 KO"
  sample_table$drug[which(sample_table$drug=="no")]="Ctrl"
  
  sample_table$genotype<-factor(sample_table$genotype, levels=c("WT", "Ube2j2 KO"))
  sample_table$drug<-factor(sample_table[lib_names,"drug"], levels=c("Ctrl","Tu"))
  sample_table$time<-factor(sample_table[lib_names,"time"], levels=c("6", "18"))
  sample_table$group<-factor(paste0(sample_table$time, " ", sample_table$drug), levels=c("6 Ctrl", "18 Ctrl", "6 Tu", "18 Tu"))

#scale data
data<-read.table("Data/processed_minprob.tsv", sep="\t", header=T, row.names=1)
row.names(data)<-data$sample
data.scale<-merge(sample_table[,c(3,1)], scale(data[,-c(1:2)]), by.x="row.names", by.y="row.names")
row.names(data.scale)<-data.scale$Row.names
data.scale<-data.scale[,-1]

#load module data
gene_modules<-read.table("Data/Protein_modules.txt", sep="\t", header=T)

PlotGeneExp<-function(gene="COPA.P53621", data=data.scale, modules=gene_modules)
{
  ind<-grep(gene, names(data))
  mod<-modules[which(modules$gene==gene),2]
  
  tmp<-data[,c(1:2,ind)]
  #tmp[,1]<-factor(tmp[,1], levels=c("WT 6hrs","WT 18hrs","WT Tu 6hrs", "WT Tu 18hrs","C6 6hrs","C6 18hrs","C6 Tu 6hrs","C6 Tu 18hrs"))
  #tmp[,2]<-as.numeric(tmp[,2])

  return(ggplot(tmp, aes(x=tmp[,2], y=tmp[,3], color=tmp[,1])) + geom_boxplot()+ ylab(paste0(gene," Expression")) + scale_x_discrete(labels = function(x) str_wrap(x, width = 2)) + xlab("Treatments") + ggtitle(paste0(gene," occurs in the ", mod, " module")) + facet_grid(.~tmp[,1]) + scale_color_manual(values = c("#7ac4e3", "#9f3b75")) + theme_bw() + guides(color = "none")) 
}

#ube2j2 color #9f3b75
#WT color #7ac4e3

# plot individual gene
PlotGeneExp(gene="COPA.P53621", data=data.scale, modules = gene_modules)

# plot all genes and save to folder
for(i in names(data.scale)[-c(1:2)])
{
  G<-PlotGeneExp(i,data.scale)
  ggsave(file=paste0("Data/plots/",i,".pdf"))
}

#### Add module plots ####
require(ggplot2); require(reshape2); require(plyr); require(stringr)

#scale data
data<-read.table("Data/processed_minprob.tsv", sep="\t", header=T, row.names=1)
data.scale<-cbind(data[,1:2],scale(data[,-c(1:2)]))

# load module data
load("Data/networkConstruction.RData")
gene_modules<-read.table("Data/Protein_modules.txt", sep="\t", header=T)

PlotModule<-function(mod="blue", data=data.scale, gene_m=gene_modules, merMEs=MEs) 
{
  indx<-which(gene_m[,2]==mod)
  
  module_exp<-data[,c(1,2,which(names(data) %in% gene_m$gene[indx]))]
  
  mod_ME<-data.frame(merMEs[,which(names(merMEs)==paste0("ME",mod))])
  names(mod_ME)<-paste0("ME",mod)
  row.names(mod_ME)<-row.names(merMEs)
  
  mod.cor<-cor(module_exp[,-c(1:2)], mod_ME)
  neg_cor<-row.names(mod.cor)[which(mod.cor<0)]
  
  mod_exp_long<-melt(module_exp, id.vars=c("group", "sample"), variable.name = "gene", value.name = "exp")
  
  mod_exp_long$MEcorrelation<-"positive"
  
  #invert the negative correlations genes
  mod_exp_long[which(mod_exp_long$gene %in% neg_cor),"exp"]<--mod_exp_long[which(mod_exp_long$gene %in% neg_cor),"exp"]
  mod_exp_long[which(mod_exp_long$gene %in% neg_cor),"MEcorrelation"]<-"negative"
  
  mod_exp_long[,1]<-factor(mod_exp_long[,1], levels=c( "WT 6hrs","WT 18hrs","WT Tu 6hrs", "WT Tu 18hrs", "C6 6hrs","C6 18hrs","C6 Tu 6hrs","C6 Tu 18hrs"))
  
  mod_exp_long_mean<-ddply(mod_exp_long, .(group, gene, MEcorrelation), summarize, mean=mean(exp))
  mod_exp_long_mean$genotype<-sub("(\\w{2}) .+", "\\1",mod_exp_long_mean$group,  perl=T)
  mod_exp_long_mean$treatment<-sub("\\w{2} (.+)", "\\1",mod_exp_long_mean$group,  perl=T)
  
  
  modAV<-ddply(mod_exp_long, .(group), summarize, mean=mean(exp))
  modAV$gene<-paste0("ME",mod)
  modAV$genotype<-sub("(\\w{2}) .+", "\\1",modAV$group,  perl=T)
  modAV$treatment<-sub("\\w{2} (.+)", "\\1",modAV$group,  perl=T)
  
  p<-ggplot(mod_exp_long_mean, aes(x=treatment, y=mean, group=gene, color=MEcorrelation)) + geom_line(alpha=.2)+ ylab(paste0(mod," Expression")) + theme_bw() +  geom_line(data=modAV, aes(x=treatment, group=gene, y=mean), color="black") + facet_grid(genotype~.)
  return(p)
}

PlotModule()
