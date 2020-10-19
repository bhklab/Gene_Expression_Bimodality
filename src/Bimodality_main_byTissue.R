library(Biobase)
library(PharmacoGx)
library(wCI)
#library(ellipse)
library(ggplot2)

# library(mixtools)
# library(mclust)
# library(scales)
# require(caret)
# library(Hmisc)
# library(psych)

setwd("~/Projects/github/bimodality_expression/")

source("normalizeCellLine.R")
source("getBiModalScore.R")
source("AzureRun/getPredictions_LOBICO.R")
source("AzureRun/crossValidate_Lobico_updated_Rversion.R")
source("AzureRun/resample_Cindex_PValue.R")


############
# Data needed
############

CCLE <- readRDS("~/Projects/PSets/final_psets/CCLE.rds")
CTRPv2 <- readRDS("~/Projects/PSets/final_psets/CTRPv2.rds")

genes_mappings <- featureInfo(CCLE,"rnaseq")


############
# Processsing Data
############
ccleRnaSeq_ALL <- t(exprs(summarizeMolecularProfiles(CCLE,mDataType = "rnaseq",fill.missing = F)))
ccleTissues <- CCLE@cell[rownames(ccleRnaSeq_ALL),"tissueid",drop=F]
#ibx <- which(ccleTissues$tissueid=="haematopoietic_and_lymphoid_tissue")
ibx <- which(ccleTissues$tissueid=="lung")
ccleTissues_WithoutHAEMATO <- ccleTissues[ibx,,drop=F]
ccleLung <- rownames(ccleTissues_WithoutHAEMATO)
ccleRnaSeq <- ccleRnaSeq_ALL[rownames(ccleTissues_WithoutHAEMATO),]
CCLE_tissues <- data.frame(table(ccleTissues[rownames(ccleRnaSeq),]),stringsAsFactors = F)

colors <- c("#3d4e2e", "#be44d5", "#6bd651", "#6142c8", "#c3e142", "#583685", "#d3ba3d", "#657cd0", "#db492a", "#65d89c", "#d44398", "#568b3a", "#bf72cd", "#bed781", "#802c57", "#9bdbca", "#d53f59", "#4b897e", "#d88334", "#73aecf", "#8d3a24", "#cba5d6", "#8a712e", "#434d6b", "#ddae75", "#966b88", "#969372", "#d37c7e", "#673e37", "#dfc0bc")
colors <- colors[1:length(table(ccleTissues_WithoutHAEMATO))]
names(colors) <- names(table(ccleTissues_WithoutHAEMATO))

set.seed(12354)
CCLE_tissues <- CCLE_tissues[sample(CCLE_tissues$Var1,size = dim(CCLE_tissues)[1]),]

pdf("Plots/pieChart_CCLE_tissues.pdf",height = 9,width = 9)
pie(CCLE_tissues$Freq,labels = NA,main = "CCLE",col = colors[CCLE_tissues$Var1])
legend("bottom",legend = names(colors[CCLE_tissues$Var1]),fill = colors[CCLE_tissues$Var1],bty="n",ncol = 4)
dev.off()

pdf("Plots/pieChart_CCLE_tissues_legend.pdf",height = 3.5,width = 13)
#pie(CCLE_tissues$Freq,labels = NA,main = "CCLE",col = colors[CCLE_tissues$Var1])
plot.new()
legend("center",legend = names(colors[CCLE_tissues$Var1]),fill = colors[CCLE_tissues$Var1],bty="n",ncol = 4)
dev.off()


library(parallel)
cl <- makeCluster(3)

BiModalScores_ccleRnaSeq_ALL <- parApply(cl = cl,FUN = function(x){source("getBiModalScore.R");getBiModalScore_Updated(x)},MARGIN = 2,X = ccleRnaSeq)
BiModalScores_ccleRnaSeq <- unlist(lapply(BiModalScores_ccleRnaSeq_ALL, "[[",1))
BiModalScores_ccleRnaSeq <- BiModalScores_ccleRnaSeq[order(BiModalScores_ccleRnaSeq,decreasing = T)]

stopCluster(cl)



protCodingGenes <- genes_mappings[names(BiModalScores_ccleRnaSeq),"gene_type"]
names(protCodingGenes) <- names(BiModalScores_ccleRnaSeq)

protCodingGenes <- protCodingGenes[protCodingGenes=="protein_coding"]
protCodingGenes <- names(protCodingGenes)


percentile <- 0.80
ccle_RS_th <- quantile(BiModalScores_ccleRnaSeq[protCodingGenes],probs = percentile)

pdf("Plots/CCLE_bimodality_dist_Lung.pdf")
hist(BiModalScores_ccleRnaSeq[protCodingGenes],main = "Dist. of bimodality score in CCLE [Lung]",xlab = "BI score",breaks = 100,col="skyblue")
#abline(v=ccle_RS_th,col="red",lty=2)
#legend("topright",legend = paste(percentile*100,"th percentile",sep = ""),lty = 2,col = "red",bty = "n")
dev.off()





ccle2 <- ccleRnaSeq#[,commonBimodalGenes]

cutoffs_ccle2 <- unlist(lapply(colnames(ccle2), function(x){
  val <- BiModalScores_ccleRnaSeq_ALL[[x]][["mix"]][["m.step"]][["mu"]]
  return(mean(val))
}))
names(cutoffs_ccle2) <- colnames(ccle2)

ccle2_binary <- getBinaryValues(ccle2,cutoffs_ccle2)


topBimodal <- names(sort(BiModalScores_ccleRnaSeq[protCodingGenes],decreasing = T)[1:100])
head(featureInfo(CCLE,"rnaseq")[topBimodal,])

pdf("Plots/CCLE_Top_bimodal_genes_Lung.pdf",height = 4,width = 6)
lapply(topBimodal[c(61,62,87)], function(x){
  
  hist(ccle2[,x],main = featureInfo(CCLE,"rnaseq")[x,"Symbol"],xlab = "Expression",breaks = 20,col="skyblue",probability = T)
  lines(density(ccle2[,x]))
  
  
})
dev.off()



#ccle_RS_bimodalGenes <- names(BiModalScores_ccleRnaSeq[protCodingGenes])[BiModalScores_ccleRnaSeq[protCodingGenes]>=ccle_RS_th]
#names(ccle_RS_bimodalGenes) <- ccle_RS_bimodalGenes
#commonBimodalGenes <- ccle_RS_bimodalGenes

#finalSetOfGenes <- commonBimodalGenes
#names(finalSetOfGenes) <-  genes_mappings[finalSetOfGenes,"Symbol"]
#finalSetOfGenes <- finalSetOfGenes[!duplicated(names(finalSetOfGenes))]
#write.table(names(finalSetOfGenes),file = paste("finalBimodalGenesSymbols_",percentile,"_CCLE.csv",sep = ""),quote = F,row.names = F,col.names = F)



#save(list = c("ccle2","ccle2_binary","finalSetOfGenes"),file = paste("AzureRun/dataNeededForLOBICO_CTRPv2_",percentile,".RData",sep = ""))

drugsCTRPv2 <- drugNames(CTRPv2)

names(drugsCTRPv2) <- 1:length(drugsCTRPv2)
names(drugsCTRPv2) <- paste0("drug",names(drugsCTRPv2))

saveRDS(drugsCTRPv2,"AzureRun/CTRPv2Drugs_dictionary.rds")

#################
# TCGA constraints
TCGA_mat <- t(readRDS("TCGA/lungMat_TCGA.rds"))


library(parallel)
cl <- makeCluster(3)

BiModalScores_TCGARnaSeq_ALL <- parApply(cl = cl,FUN = function(x){source("getBiModalScore.R");getBiModalScore_Updated(x)},MARGIN = 2,X = TCGA_mat)
BiModalScores_TCGARnaSeq <- unlist(lapply(BiModalScores_TCGARnaSeq_ALL, "[[",1))
BiModalScores_TCGARnaSeq <- BiModalScores_TCGARnaSeq[order(BiModalScores_TCGARnaSeq,decreasing = T)]

stopCluster(cl)


hist(BiModalScores_TCGARnaSeq)
TCGA_threshold <- quantile(BiModalScores_TCGARnaSeq[protCodingGenes],probs = 0.8)
ccle_Threshold <- quantile(BiModalScores_ccleRnaSeq[protCodingGenes],probs = 0.8)


ibx <- which(BiModalScores_ccleRnaSeq[names(BiModalScores_TCGARnaSeq[protCodingGenes])] >= ccle_Threshold & BiModalScores_TCGARnaSeq[protCodingGenes] >= TCGA_threshold  )

df1 <- data.frame(data.frame("TCGA"=BiModalScores_TCGARnaSeq[protCodingGenes],"CCLE"=BiModalScores_ccleRnaSeq[protCodingGenes],"group"=(BiModalScores_TCGARnaSeq[protCodingGenes]>TCGA_threshold & BiModalScores_ccleRnaSeq[protCodingGenes]>ccle_Threshold)))



#A <- 0.4
#df_ell <- NULL
#df_ell <- rbind(df_ell, cbind(as.data.frame(with(df1[df1$group==1,], ellipse(A, 
#                                                                             scale=c(sd(TCGA)*1.5,sd(CCLE)*1.5), 
#                                                                             centre=c(quantile(TCGA,probs = 0.8),quantile(CCLE,probs = 0.8)),t=2.3
#)
#)
#),
#group=1))

myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
myColor_scale_fill <- scale_fill_gradientn(colours = myColor)

pdf("Plots/Bimodality_CCLE_TCGA_Lung.pdf",width = 11,height = 10)
ggplot(df1, aes(TCGA, CCLE)) + stat_binhex(aes(fill=log(..count..)))  + labs(y= "CCLE bimodality score", x = "TCGA bimodality score") +
  xlim(0, 5) + ylim(0, 5) +
  myColor_scale_fill +
  geom_hline(yintercept = ccle_Threshold, colour = "red",lty=2,show.legend = F) +
  geom_vline(xintercept = TCGA_threshold, colour = "red",lty=2, show.legend = F) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank()) 
#geom_path(data=df_ell, aes(x=x, y=y,colour=group), size=1, linetype=1,show.legend = F) 
dev.off()

pdf("Plots/bimodality_legend_Lung.pdf")
plot.new()
legend("center",legend = "80th percentile",lty = 2,col = "red",bty = "n")
finalSetOfGenes <- names(ibx)
names(finalSetOfGenes) <- genes_mappings[finalSetOfGenes,"Symbol"]

save(list = c("ccle2","ccle2_binary","finalSetOfGenes"),file = "AzureRun/dataNeededForLOBICO_CTRPv2_TCGA_0.8_LUNG.RData")





#load("dataGeneratedToThisPoint_22Jul2019_Blood.RData")


################################
#

finalSetOfGenes_Lung <- finalSetOfGenes

common_features <- intersect(finalSetOfGenes,finalSetOfGenes_Lung)
names(common_features) <- names(finalSetOfGenes)[match(common_features,finalSetOfGenes)]
LungFeatures <- setdiff(finalSetOfGenes_Lung,finalSetOfGenes)
names(LungFeatures) <- names(finalSetOfGenes_Lung)[match(LungFeatures,finalSetOfGenes_Lung)]
write.table(names(LungFeatures),file = "Lung_Specific_features.csv",sep = ",",quote = F,col.names = F,row.names = F)

save.image("dataGeneratedToThisPoint_31Oct2019_Lung.RData")


load("dataGeneratedToThisPoint_31Oct2019_Lung.RData")



ccle2_final <- ccle2[,finalSetOfGenes]
colnames(ccle2_final) <- finalSetOfGenes





A <- cor(ccle2_final)

library(ComplexHeatmap)

pdf("Plots/correlation of best bimodal genes_final.pdf")
Heatmap(A,show_row_names = F,show_column_names = F,cluster_columns = T,cluster_rows = T,name = "correlation",column_title = "Correlations between final list of bimodal genes")
dev.off()

png("Plots/correlation of best bimodal genes_final.png")
Heatmap(A,show_row_names = F,show_column_names = F,cluster_columns = T,cluster_rows = T,name = "correlation",column_title = "Correlations between final list of bimodal genes")
dev.off()

ccle2_bin_final <- ccle2_binary[,finalSetOfGenes]
colnames(ccle2_bin_final) <- finalSetOfGenes

mccCCLE <- matrix(nrow = dim(ccle2_bin_final)[2],ncol = dim(ccle2_bin_final)[2])
for (i in 1:(dim(ccle2_bin_final)[2]-1)) {
  for (j in (i+1):dim(ccle2_bin_final)[2]) {
    tmpEst <- mccr(as.factor(ccle2_bin_final[,i]),as.factor(ccle2_bin_final[,j]))
    mccCCLE[i,j] <- tmpEst
    mccCCLE[j,i] <- tmpEst
  }
}


for (i in 1:(dim(ccle2_bin_final)[2])) {
  for (j in 1:dim(ccle2_bin_final)[2]) {
    if(i==j){
      mccCCLE[j,i] <- 1
    }
  }
}

row.clusters = hclust(dist(mccCCLE))
clusters <- cutree(row.clusters,h = 7)
names(clusters) <- colnames(ccle2_bin_final)

bimodality_groups_of_genes <- list()

l <- length(table(clusters))
for (i in 1:13) {
  ibx <- clusters[clusters==i]
  bimodality_groups_of_genes[[i]] <- names(ibx)
}

mccCCLE_tmp <- mccCCLE[row.clusters$order,row.clusters$order]
for (i in 1:(dim(ccle2_bin_final)[2]-1)) {
  for (j in (i+1):dim(ccle2_bin_final)[2]) {
    mccCCLE_tmp[j,i] <- NA
  }
}


library(ComplexHeatmap)
library(circlize)



pdf("Plots/MCC correlation of best bimodal genes_final2.pdf")
Heatmap(mccCCLE_tmp,show_row_names = F,show_column_names = F,cluster_columns = F,cluster_rows = F,show_row_dend = F,name = "correlation",column_title = "Correlations between final list of bimodal genes"
        ,col = colorRamp2(c(-1, 0, 1),colors = c("red","white","blue")),na_col = "white")
dev.off()

png("MCC correlation of best bimodal genes_final.png")
Heatmap(mccCCLE_tmp,show_row_names = F,show_column_names = F,cluster_columns = F,cluster_rows = F,name = "correlation",column_title = "Correlations between final list of bimodal genes"
        ,col = colorRamp2(c(-1, 0, 1),colors = c("red","white","blue")),na_col = "white")
dev.off()




A2 <- apply(ccle2[,finalSetOfGenes], 2, scale)

samples <- data.frame("samples"=CCLE@cell[rownames(ccleTissues_WithoutHAEMATO),"cellid"],stringsAsFactors = F)
samples$tissue = CCLE@cell[rownames(ccleTissues_WithoutHAEMATO),"tissueid"]
rownames(samples) <- samples$samples

#clrs <- c("#ff6fdc","#49d457","#c675ff","#bdca12","#b20085","#96d94f","#ff3386","#02d695","#f93650","#008c5d","#bc0065","#3c8100","#a797ff","#cb9a00","#0278d3","#ff9240","#1ca8ff","#bf0028","#77a3ff","#948600","#6c427b","#c6cd61","#017eab","#ff6358","#a2d38d","#a3162e","#606f00","#ff8b9c","#7b6837","#fab0a7","#8b3a06","#8b3640","#7b6837","#fab0a7","#8b3a06","#8b3640")

colrs <- data.frame("tissue"=names(colors), "Tcolr"=colors,stringsAsFactors = F)


samples <- cbind(samples,"colrs"=colrs[match(samples[,"tissue"],colrs$tissue),"Tcolr"])

colrs1 <- colrs[,"Tcolr"]
names(colrs1) <- colrs$tissue
ha <- HeatmapAnnotation(samples[rownames(ccle2_binary),"tissue",drop=F],which = "column",col = list(tissue=colrs1))

pdf("Plots/clustering of samples based on best bimodal genes_final.pdf",width = 25,height = 15)
Heatmap(t(A2[,finalSetOfGenes]),show_row_names = F,show_column_names = F,cluster_columns = T,cluster_rows = T,top_annotation = ha
        ,name = "expression",column_title = "Clustering of samples based on final list of bimodal genes") 
dev.off()

png("clustering of samples based on best bimodal genes_final.png",width = 25,height = 15)
Heatmap(t(A2[,finalSetOfGenes]),show_row_names = F,show_column_names = F,cluster_columns = T,cluster_rows = T,top_annotation = ha
        ,name = "expression",column_title = "Clustering of samples based on final list of bimodal genes") 
dev.off()



pdf("Plots/clustering of samples based on best bimodal genes_final_binary.pdf",width = 25,height = 15)
Heatmap(t(ccle2_binary[,finalSetOfGenes]),show_row_names = F,show_column_names = F,cluster_columns = T,cluster_rows = T,top_annotation = ha
        ,name = "expression",column_title = "Clustering of samples based on final list of bimodal genes") 
dev.off()



library(GSA) 
library(piano)

gsc1 <- loadGSC("~/Projects/general_functions/GeneSets_MsigDB_20180823/h.all.v6.2.symbols.gmt")
gsc1 <- loadGSC("~/Projects/general_functions/GeneSets_MsigDB_20180823/c5.bp.v6.2.symbols.gmt")
gsc1 <- loadGSC("~/Projects/general_functions/GeneSets_MsigDB_20180823/c2.cp.reactome.v6.2.symbols.gmt")
gsc1 <- loadGSC("~/Projects/general_functions/GeneSets_MsigDB_20180823/c2.cp.kegg.v6.2.symbols.gmt")





finalSetOfGenes_tmp <- finalSetOfGenes[!duplicated(names(finalSetOfGenes))]

gsea_out <- piano::runGSAhyper(genes = names(finalSetOfGenes_tmp), gsc=gsc1,universe = unique(genes_mappings[protCodingGenes,"Symbol"]))

gsea_out_final <- as.data.frame(gsea_out$resTab)

gsea_out_final$fdr <- p.adjust(gsea_out_final$`p-value`,method = "fdr")

gsea_out_final <- gsea_out_final[order(gsea_out_final$`p-value`),]

(gsea_out_final[1:50,])


gsea_out_final1 <- gsea_out_final[gsea_out_final$fdr<0.05,]

gsea_out_final1 <- gsea_out_final1[order(gsea_out_final1$fdr,decreasing = T),]

pdf("Bimodal_genes_GSEA.pdf",width = 11,height = 8)
par(mai=c(1,7,1,1))
#barplot(-log10(gsea_out_final1$fdr),horiz = T,names.arg = gsub("REACTOME_","",rownames(gsea_out_final1)),las=2,xlab = "-log10(FDR)")
barplot(-log10(gsea_out_final1$fdr),horiz = T,names.arg =rownames(gsea_out_final1),las=2,xlab = "-log10(FDR)")
dev.off()


write.table(rownames(gsea_out_final1),file = "Bimdal_genes_GSEA_REACTOM.txt",sep = ",",col.names = F,row.names = F,quote = F)


write.table(names(finalSetOfGenes_tmp),"finalGenes_bimodal_forPathways.csv",col.names = F,row.names = F)


###
# per cluster of genes
library(GSA) 
library(piano)

gsc1 <- loadGSC("~/Projects/general_functions/GeneSets_MsigDB_20180823/h.all.v6.2.symbols.gmt")
gsc1 <- loadGSC("~/Projects/general_functions/GeneSets_MsigDB_20180823/c5.bp.v6.2.symbols.gmt")
gsc1 <- loadGSC("~/Projects/general_functions/GeneSets_MsigDB_20180823/c2.cp.reactome.v6.2.symbols.gmt")
gsc1 <- loadGSC("~/Projects/general_functions/GeneSets_MsigDB_20180823/c2.cp.kegg.v6.2.symbols.gmt")


gsea_out_final1_clusters <- lapply(bimodality_groups_of_genes, function(genesSet){
  
  ibx <- match(genesSet,finalSetOfGenes)
  names(genesSet) <- names(finalSetOfGenes)[ibx]
  
  finalSetOfGenes_tmp <- genesSet[!duplicated(names(genesSet))]
  
  gsea_out <- piano::runGSAhyper(genes = names(finalSetOfGenes_tmp), gsc=gsc1,universe = unique(genes_mappings[protCodingGenes,"Symbol"]))
  
  gsea_out_final <- as.data.frame(gsea_out$resTab)
  
  gsea_out_final$fdr <- p.adjust(gsea_out_final$`p-value`,method = "fdr")
  
  gsea_out_final <- gsea_out_final[order(gsea_out_final$`p-value`),]
  
  gsea_out_final1 <- gsea_out_final[gsea_out_final$fdr<0.05,]
  
  return(gsea_out_final1)
  
})


saveRDS(gsea_out_final1_clusters,"gsea_out_final1_clusters_GOBP.rda")

gsea_out_final1_clusters <- readRDS("gsea_out_final1_clusters_GOBP.rda")
gsea_out_final1_clusters <- readRDS("gsea_out_final1_clusters_REACTOM.rda")

#############################################
#############################################
#############################################
# Models over CTRPv2 

load("AzureRun/dataNeededForLOBICO_CTRPv2_TCGA_0.8_LUNG.RData",verbose = T)

AAC_CTRPv2 <- summarizeSensitivityProfiles(CTRPv2,sensitivity.measure = "aac_recomputed",fill.missing = F)

keptDrugs_Lung <- apply(AAC_CTRPv2, 1, function(x){
  x <- x[!is.na(x)]
  x <- x[intersect(names(x),ccleLung)]
  ((sum(x>0.2,na.rm = T)/sum(x>=0,na.rm = T)) >0.1)
})


#load("AzureRun/dataNeededForLOBICO_CTRPv2_TCGA.RData")
#load("AzureRun/dataNeededForLOBICO_CTRPv2_TCGA_0.75.RData")
#load("AzureRun/dataNeededForLOBICO_CTRPv2_0.95.RData")



#featuresType <- "output_mCI_CCLE_TCGA_90_F"
#featuresType <- "output_mCI_CCLE_95_F"
#featuresType <- "output_mCI_CCLE_TCGA_75_F"

featuresType <- "output_mCI_CCLE_TCGA_75_LUNG"



listOfdrugsFiels <- dir(paste("AzureRun/",featuresType,"/",sep = ""))

output_ensemble_CTRPv2_5CV3NS_LUNG <- lapply(listOfdrugsFiels, function(x){
  print(x)
  output <- readRDS(paste("AzureRun/",featuresType,"/",x,sep = ""))
  return(output)
})


drugsProcessed <- unlist(lapply(strsplit(listOfdrugsFiels,"_CTRPv2"),"[[",1))
drugsProcessed <- gsub("%2C",",",drugsProcessed)
drugsProcessed <- gsub("%5B","[",drugsProcessed)
drugsProcessed <- gsub("%5D","]",drugsProcessed)
drugsProcessed <- gsub("%20"," ",drugsProcessed)
drugsProcessed <- gsub("%7B","{",drugsProcessed)
drugsProcessed <- gsub("%7D","}",drugsProcessed)

names(output_ensemble_CTRPv2_5CV3NS_LUNG) <- drugsProcessed



output_ensemble_CTRPv2_5CV3NS_LUNG_final <- output_ensemble_CTRPv2_5CV3NS_LUNG[intersect(names(keptDrugs_Lung)[keptDrugs_Lung],names(output_ensemble_CTRPv2_5CV3NS_LUNG))]
#output_ensemble_CTRPv2_5CV3NS_LUNG_final <- output_ensemble_CTRPv2_5CV3NS_LUNG_ML_single[intersect(names(keptDrugs)[keptDrugs],names(output_ensemble_CTRPv2_5CV3NS_LUNG_ML_single))]
mciTh = 0.6
#output_ensemble_CTRPv2_5CV3NS_LUNG <- readRDS("output_ensemble_CTRPv2_5CV3NS_LUNG.rda")

output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation <- lapply(output_ensemble_CTRPv2_5CV3NS_LUNG_final,function(x){
  return(x[[1]])
})


output_ensemble_CTRPv2_5CV3NS_LUNG_CIs <- unlist(lapply(1:length(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation), function(x){output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[[x]][2]}))

#output_ensemble_CTRPv2_5CV3NS_LUNG_CIs <- unlist(lapply(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation, "[[",2))

output_ensemble_CTRPv2_5CV3NS_LUNG_CIs <- as.numeric(output_ensemble_CTRPv2_5CV3NS_LUNG_CIs)
names(output_ensemble_CTRPv2_5CV3NS_LUNG_CIs) <- names(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation)
output_ensemble_CTRPv2_5CV3NS_LUNG_CIs <- output_ensemble_CTRPv2_5CV3NS_LUNG_CIs[!is.na(output_ensemble_CTRPv2_5CV3NS_LUNG_CIs)]

output_ensemble_CTRPv2_5CV3NS_LUNG_final <- output_ensemble_CTRPv2_5CV3NS_LUNG_final[names(output_ensemble_CTRPv2_5CV3NS_LUNG_CIs)]

output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation <- output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[names(output_ensemble_CTRPv2_5CV3NS_LUNG_CIs)]
output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation <- do.call(rbind,output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation)
rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation) <- names(output_ensemble_CTRPv2_5CV3NS_LUNG_CIs)
class(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation) <- "numeric"


output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation <- output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[order(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[,"mCI"],decreasing = T),,drop=F]
output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation <- cbind(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation,"fdr"=p.adjust(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[,"Pvalue"]))
#output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[,"fdr"] <- p.adjust(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[,"Pvalue.mCI"],method = "bonferroni")

sum(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[,"mCI"] > 0.6,na.rm = T)/dim(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation)[1]


library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

drugClasses <- read.csv("CTRPv2Drugs_BasedOnGDSC1000_Nehme.csv",sep = "\t", header = T,row.names = 1,stringsAsFactors = F)


drugClasses$PossibleGDSC1000Class_final <- (apply(drugClasses[,c("PossibleGDSC1000Class","X","X.1")],1,function(x){ return(trimws(paste(x[c(1,2,3)],collapse  = " "), "r"))}))

drugClasses[intersect(grep("romat",drugClasses$PossibleGDSC1000Class_final),grep(":",drugClasses$PossibleGDSC1000Class_final,invert = T)),"PossibleGDSC1000Class_final" ] <- "chromatin"
drugClasses[grep("chromain",drugClasses$PossibleGDSC1000Class_final),"PossibleGDSC1000Class_final" ] <- "chromatin"


druginfo_CTRPv2 <- read.csv("~/Projects/general_functions/drugs_with_ids.csv",header = T,row.names = 1,sep = ",",stringsAsFactors = F)

ibx_final <- lapply(1:dim(drugClasses)[1], function(x){
  drug <- rownames(drugClasses)[x]
  ibx <- apply(druginfo_CTRPv2, 1, function(y){
    if(toupper(drug) %in% toupper(y)){
      return(T)
    }else{
      return(F)
    }
  })
  
  return(which(ibx==T))
})

sum(unlist(ibx_final)==0)
unlist(ibx_final)
which(duplicated(unlist(ibx_final)))
ibx_final_uniq <- unlist(lapply(ibx_final, function(x){ if(length(x)!=0){return(x[1])}else{return(NA)}}))
drugClasses$new_drug_id <- druginfo_CTRPv2[ibx_final_uniq,"unique.drugid"]
drugClasses[c(293),"new_drug_id"] <- "Pazopanib"
drugClasses[c(472),"new_drug_id"] <- rownames(drugClasses)[472] 

output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses <- as.data.frame(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation,stringsAsFactors=F)
rownames(drugClasses) <- drugClasses$new_drug_id
write.table(drugClasses,file = "CTRPv2Drugs_BasedOnGDSC1000_Nehme_UPDATED.csv",quote = F,col.names = NA,row.names = T,sep = ",")
output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses <- cbind(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses
                                                              ,"Class"=drugClasses[rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses),"PossibleGDSC1000Class_final"],stringsAsFactors=F)





#set.seed(1252)
cols_DrugClasses <- c(` other` = "#004700", `ABL signaling` = "#004743", `apoptosis regulation` = "#b1ff64", 
                      `cell cycle` = "#f20040", chromatin = "#0185ff", cytoskeleton = "#707d00", 
                      `DNA replication` = "#b900a3", `EGFR signaling` = "#002562", 
                      `ERK MAPK signaling` = "#f8c500", `Genome integrity` = "#01e3a8", 
                      `IGFR signaling` = "#8c26d9", metabolism = "#01b631", mitosis = "#36000e", 
                      other = "#d7d3ff", `p53 pathway` = "#a1004d", `PI3K signaling` = "#483400", 
                      `RTK signaling` = "#00d9e5", `TOR signaling` = "#ffeea1")

output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses <- cbind(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses,"ClassCol"=rep(NA,dim(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses)[1]),stringsAsFactors=F)

for (i in 1:dim(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses)[1]) {
  output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses[i,"ClassCol"] <- cols_DrugClasses[output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses[i,"Class"]]
}

ibx <- output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses$mCI>0.6
table(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses[ibx,"Class"])

pdf("Plots/CTRPv2 - cross validation results - mCI LUNG.pdf",width = 12)
data1 <- sort(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[,"mCI"])
plotdata_tmp <- data1 - 0.5
p <- barplot(plotdata_tmp,las=2,col = ifelse(data1>mciTh,"#5ab4ac",ifelse(data1>0.6,"#5ab4ac","#d8b365"))
             ,main = "CTRPv2 - cross validation results [mCI]",border = NA,names.arg = NA,xlab = "Drugs",ylab="mCI",axes=FALSE
             ,ylim = c(min(plotdata_tmp,na.rm = T)-0.05,max(plotdata_tmp,na.rm = T)+0.05))

axis(side=2,at=seq(0,1,0.1)-0.5,labels=(seq(0,1,0.1)))

legend("topleft",legend = c("mCI > 0.6","mCI <= 0.6","FDR < 0.05")
       ,col = c("#5ab4ac","#d8b365","black"),bty = "n",pch=c(15,15,16))

points(p,plotdata_tmp+0.01,pch=ifelse(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[names(data1),"fdr"]<0.05,16,NA))
points(p,rep(-0.005,dim(p)[1]),pch=ifelse(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses[names(data1),"mCI"]>0.6,15,NA)
       ,col=output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses[names(data1),"ClassCol"])

dev.off()

pdf("Plots/drugClasses_CTRPv2_LUNG.pdf",width = 13,height = 4)
plot.new()
legend("center",legend = names(table(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses[ibx,"Class"]))
       ,fill = cols_DrugClasses[names(table(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation_drugClasses[ibx,"Class"]))],ncol=5,bty="n")
dev.off()

sum(data1>=mciTh,na.rm = T)/  sum(data1>-1,na.rm = T)

sum(data1>=0.65,na.rm = T)/  sum(data1>-1,na.rm = T)

topHits <- rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation)[output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[,"mCI"] >= 0.7]

# Top hits and biology:
# pdf("topHits_CTRPv2.pdf")
# topHits <- rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation)[output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[,"mCI"] >= 0.6]
# 
# lapply(topHits, function(x){
#   data <- output_ensemble_CTRPv2_5CV3NS_LUNG[[x]][[4]]
#   highS <- rownames(data)[data[,"Vote"]==1]
#   lowS <- rownames(data)[data[,"Vote"]==0]
#   
#   boxplot(list(AUC[x,highS],AUC[x,lowS]),main=x)
# })
# dev.off()





results_final_common <- lapply(topHits, function(x,output,experssion,AUCs_Drugs){
  
  AUCs_Drug <- AUCs_Drugs[x,]
  AUCs_Drug <- AUCs_Drug[!is.na(AUCs_Drug)]
  
  common <- intersect(names(AUCs_Drug),rownames(experssion))
  
  
  final_results <- lapply(1:length(output[[x]][[3]]), function(y){
    print(y)
    SolMat <- output[[x]][[3]][[y]][[2]]
    K <- as.numeric(output[[x]][[3]][[y]][[3]])
    M <- as.numeric(output[[x]][[3]][[y]][[4]])
    print((match(colnames(SolMat),colnames(experssion))))
    results <- getPredictions_LOBICO_outer(SolMat,cbind("drug_AUC"=AUCs_Drug[common], experssion[common,colnames(SolMat)]),K,M)
  })
  
  final_results <- do.call(cbind,final_results)
  
  final_results <- cbind(final_results[,seq(1,length(output[[x]][[3]])*2,2)],final_results[,2])
  
  finalPredictions <- apply(final_results[,1:length(output[[x]][[3]]),drop=F], 1, function(x){
    round(mean(x))
  })
  
  return(cbind("Vote"=finalPredictions,"Obs"=final_results[,dim(final_results)[2]]))
},output=output_ensemble_CTRPv2_5CV3NS_LUNG_final,experssion=ccle2_binary[,finalSetOfGenes],AUCs_Drugs=AAC_CTRPv2)


names(results_final_common) <- topHits


pdf("Plots/topHits_CTRPv2_barPlots_Lung.pdf",width = 12,height = 6)
lapply(topHits, function(x){
  results <- results_final_common[[x]]
  results <- results[order(results[,2],decreasing = T),]
  results_tmp <- results
  results_tmp[,2] <- results_tmp[,2]-0.2
  p <- barplot(results_tmp[,2],col = ifelse(results_tmp[,1]==1,"#75d2ff","#f24720"),border = NA,axes=FALSE,names.arg = NA,ylab = "AAC",ylim = c(-0.3,max(results_tmp[,2])))
  axis(side=2,at=seq(0,1,0.1)-0.2,labels=(seq(0,1,0.1)))
  legend("topright",legend = c("Sensitive","Resistant"),fill = c("#75d2ff","#f24720"),bty="n")
  title(main = paste(x,"[CTRPv2]\n","mCI:",sprintf("%.2f",as.numeric(output_ensemble_CTRPv2_5CV3NS_LUNG_final[[x]][[1]][,"mCI"]))))
  
  #  tmpTissues <- ccleTissues[rownames(results_tmp),1]
  #  ibx <- match(tmpTissues,names(colors))
  #  points(p,rep(-.25,dim(results_tmp)[1]),pch=16,col=colors[ibx],cex=0.5)
  #  legend("right",legend = names(colors),fill = colors,bty="n")
})

dev.off()




rules <- lapply(topHits,function(drug){
  
  
  rule <- unique(unlist(output_ensemble_CTRPv2_5CV3NS_LUNG[[drug]][[2]]))
  #    rule <- unique(unlist(lapply(rule, function(y){
  #      y <- gsub("[&~|]","",y)
  #      r <- unlist(strsplit(y,split = " "))
  #      r <- r[r!=""]
  #      r_symbols <- genes_mappings[r,"Symbol"]
  #    })))
  
  return(rule)
})
names(rules) <- topHits

rules_symbols <- lapply(rules, function(r){
  rules_symbols_tmp <- c()
  for (rule in r) {
    genes <- strapplyc(rule, "ENSG[0-9]+\\.[0-9]+")[[1]]
    
    rules_symbols_tmp <- c(rules_symbols_tmp,lapply(genes, function(s){
      rule <<- gsub(s,names(finalSetOfGenes)[match(s,finalSetOfGenes)],rule)
    })[[length(genes)]]
    )
    
  }
  rules_symbols_tmp <- gsub("~","\U00AC",rules_symbols_tmp)
  return(rules_symbols_tmp)
})



rules_mapping <- unique(unlist(rules))
rules_mapping <- gsub("[&~|]","",rules_mapping)
rules_mapping <- unlist(strsplit(rules_mapping,split = " "))
rules_mapping <- rules_mapping[rules_mapping!=""]
rules_mapping <- unique(rules_mapping)
rules_mapping_symbols <- genes_mappings[rules_mapping,"Symbol"]
cbind(rules_mapping,rules_mapping_symbols)



pdf("Plots/topHits_CTRPv2_barPlots_withRules_LUNG.pdf",width = 12,height = 6)
lapply(topHits, function(x){
  print(x)
  results <- results_final_common[[x]]
  results <- results[order(results[,2],decreasing = T),]
  results_tmp <- results
  results_tmp[,2] <- results_tmp[,2]-0.2
  barplot(results_tmp[,2],col = ifelse(results_tmp[,1]==1,"#75d2ff","#f24720"),border = NA,axes=FALSE,names.arg = NA,ylab = "AAC")
  axis(side=2,at=seq(0,1,0.1)-0.2,labels=(seq(0,1,0.1)))
  legend("topright",legend = c("Sensitive","Resistant"),fill = c("#75d2ff","#f24720"),bty="n")
  title(main = paste(x,"[CTRPv2]\n","mCI:",sprintf("%.2f",as.numeric(output_ensemble_CTRPv2_5CV3NS_LUNG_final[[x]][[1]][,"mCI"])),
                     "\n",paste("[",rules_symbols[[x]],"]",collapse = " | ")))
  #  tmpTissues <- ccleTissues[rownames(results_tmp),1]
  #  ibx <- match(tmpTissues,names(colors))
  #  points(p,rep(-.25,dim(results_tmp)[1]),pch=16,col=colors[ibx],cex=0.5)
  #  legend("right",legend = names(colors),fill = colors,bty="n")
})

dev.off()



pdf("Plots/topHits_CTRPv2_boxPlots_withRules.pdf",width = 12,height = 6)
lapply(topHits[1:10], function(x){
  print(x)
  results <- results_final_common[[x]]
  results <- results[order(results[,2],decreasing = T),]
  
  tmpTissues <- ccleTissues[rownames(results),1]
  ibx <- match(tmpTissues,names(colors))
  
  results <- cbind(as.data.frame(results),as.data.frame(tmpTissues),colors[ibx])
  results$Vote <- as.factor(results$Vote)
  ggplot(results,mapping = aes(tmpTissues,Obs)) + geom_violin(trim = T) + 
    geom_boxplot(width = 0.1,position = position_dodge(0.8)) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
    geom_jitter(
      aes(fill = Vote,col=Vote),
      position =position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
    ) +
    ggtitle(paste(x,"[CTRPv2]\n","mCI:",sprintf("%.2f",as.numeric(output_ensemble_CTRPv2_5CV3NS_LUNG_final[[x]][[1]][,"mCI"])),
                  "\n",rules_manulally_curated[[x]]))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5))
  
  
})

dev.off()


library(circlize)
library(ComplexHeatmap)

pdf("topHits_CTRPv2_Fancy.pdf")

lapply(topHits[1:10], function(x){
  print(x)
  results <- as.data.frame(results_final_common[[x]])
  results <- results[order(results[,2],decreasing = T),]
  
  col_fun = colorRamp2(c(0, max(results$Obs)), c("blue", "red"))
  
  ha1 <- HeatmapAnnotation(Obs = anno_barplot(results[,"Obs"], axis = TRUE,border = F,axis_side = "left",gp = gpar(at = 1:length(results$Obs),col = NA, fill = col_fun(results$Obs))),annotation_height = unit(c(3), "cm")
  )
  
  
  tissuesSub <- names(table(ccleTissues_WithoutHAEMATO[rownames(results),,drop=F]))
  
  colrs <- colors[tissuesSub]
  
  ha2 <- HeatmapAnnotation(cbind(results[,"Vote",drop=F],ccleTissues_WithoutHAEMATO[rownames(results),,drop=F])
                           , col = list(Vote = c("1" = "red", "0" = "blue"),
                                        tissueid = colrs)
                           ,annotation_legend_param = list(
                             Vote = list(labels=c("Sensitive", "Resistant"),ncol=1,title="Prediction"),
                             tissueid = list(ncol=4)
                           )
  )
  
  #  ha3 <- HeatmapAnnotation(ccleTissues_WithoutHAEMATO[rownames(results),,drop=F])
  
  
  zero_row_mat = matrix(c(1:dim(results)[1]),nrow = 1, ncol = dim(results)[1])
  zero_row_mat = matrix(nrow = 0, ncol = dim(results)[1])
  
  colnames(zero_row_mat) <- rownames(results)
  
  lgd = Legend(at = seq(0,max(results$Obs),0.1), col_fun = col_fun, title = "AAC",direction = "horizontal")
  
  ht <- Heatmap(zero_row_mat,top_annotation = ha1,bottom_annotation = ha2,na_col = "gray",show_row_names = F,show_column_names = F,cluster_columns = F,cluster_rows = F,column_title = x
  )
  
  
  draw(ht, padding = unit(c(2, 20, 20, 10), "mm"),annotation_legend_side = "bottom",heatmap_legend_side = "bottom",heatmap_legend_list = list(lgd))
  
  decorate_annotation("Obs", {
    grid.text("AAC", unit(-10, "mm"), just = "bottom", rot = 90)
    grid.lines(c(0, 1), unit(c(0.2, 0.2), "native"), gp = gpar(lty = 2, col = "pink"))
  })
  
  decorate_annotation("Vote", {
    grid.text("Prediction", unit(-10, "mm"), just = "bottom", rot = 0)
    #   grid.lines(c(0, 1), unit(c(0.2, 0.2), "native"), gp = gpar(lty = 2, col = "pink"))
  })
  decorate_annotation("tissueid", {
    grid.text("Tissue", unit(-10, "mm"), just = "bottom", rot = 0)
    #   grid.lines(c(0, 1), unit(c(0.2, 0.2), "native"), gp = gpar(lty = 2, col = "pink"))
  })
  
})
dev.off()


#################
# Drugs classes and pathways

drugClasses <- read.csv("CTRPv2Drugs_BasedOnGDSC1000_Nehme.csv",sep = "\t", header = T,row.names = 1,stringsAsFactors = F)

drugClasses$PossibleGDSC1000Class[intersect(grep("romat",drugClasses$PossibleGDSC1000Class),grep(":",drugClasses$PossibleGDSC1000Class,invert = T)) ] <- "chromatin"
drugClasses$PossibleGDSC1000Class[grep("chromain",drugClasses$PossibleGDSC1000Class) ] <- "chromatin"


ibx <- rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation)[output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[,"mCI"] >= 0.6]

drugClassesInTopHits <- names(table(drugClasses[ibx,"PossibleGDSC1000Class"]))

drugClassesInTopHits_associated_rules <- lapply(drugClassesInTopHits, function(x){
  drugs <- intersect(rownames(drugClasses)[drugClasses$PossibleGDSC1000Class==x],ibx)
  rules <- lapply(drugs,function(drug){
    
    
    rule <- output_ensemble_CTRPv2_5CV3NS_LUNG[[drug]][[2]]
    rule <- unique(unlist(lapply(rule, function(y){
      y <- gsub("[&~|]","",y)
      r <- unlist(strsplit(y,split = " "))
      r <- r[r!=""]
      r_symbols <- genes_mappings[r,"Symbol"]
    })))
    
  })
  names(rules) <- drugs
  return(rules)
})

names(drugClassesInTopHits_associated_rules) <- drugClassesInTopHits



library(GSA) 
library(piano)

gsc1 <- loadGSC("~/Projects/general_functions/GeneSets_MsigDB_20180823/h.all.v6.2.symbols.gmt")
gsc1 <- loadGSC("~/Projects/general_functions/GeneSets_MsigDB_20180823/c5.bp.v6.2.symbols.gmt")
gsc1 <- loadGSC("~/Projects/general_functions/GeneSets_MsigDB_20180823/c2.cp.reactome.v6.2.symbols.gmt")
gsc1 <- loadGSC("~/Projects/general_functions/GeneSets_MsigDB_20180823/c2.cp.kegg.v6.2.symbols.gmt")


universe = unique(genes_mappings[protCodingGenes,"Symbol"])


library(EnrichmentBrowser)
drugClassesInTopHits_associated_rules_pathways <- lapply(drugClassesInTopHits_associated_rules, function(x){
  
  #  pathwaysPerDrug <- lapply(x, function(y){
  #  hsa.gs <- getGenesets("~/Projects/general_functions/GeneSets_MsigDB_20180823/c2.cp.reactome.v6.2.symbols.gmt")
  #  sbea.res <- sbea(method="ora", se=allSE, gs=hsa.gs, perm=0, alpha=0.1)
  #  gsRanking(sbea.res)$FDR <- p.adjust(gsRanking(sbea.res))
  #  gsRanking(sbea.res)
  
  gsea_out <- piano::runGSAhyper(genes = unique(unlist(x)), gsc=gsc1,universe = universe)
  gsea_out_final <- as.data.frame(gsea_out$resTab)
  gsea_out_final$fdr <- p.adjust(gsea_out_final$`p-value`,method = "fdr")
  gsea_out_final <- gsea_out_final[order(gsea_out_final$`p-value`),]
  #  })
  
  
})

capture.output(lapply(drugClassesInTopHits_associated_rules_pathways, function(x){
  head(x[,c("p-value","Adjusted p-value","fdr")])
}), file = "DrugClassesPathways.txt")

#drugClassesInTopHits_associated_rules_pathways_GO <- drugClassesInTopHits_associated_rules_pathways
#drugClassesInTopHits_associated_rules_pathways_Reactom <- drugClassesInTopHits_associated_rules_pathways

#############################
# GDSC1000

#load("~/Projects/PSets/PSets_updated/GDSCv2.RData",verbose = T)
#load("~/Projects/PSets/PSets/GDSC1000_kallisto_toil2.RData",verbose = T)
#GDSC1000 <- GDSC
#rm(GDSC)
#gc()

GDSC1000 <- readRDS("~/Projects/PSets/final_psets/GDSCv2.rds")
gdsc1000RNAseq_ALL <- t(exprs(summarizeMolecularProfiles(GDSC1000,mDataType = "rnaseq",fill.missing = F)))
gdsc1000Tissues <- GDSC1000@cell[rownames(gdsc1000RNAseq_ALL),"tissueid",drop=F]
ibx <- which(gdsc1000Tissues$tissueid=="lung")
gdsc1000Tissues_WithoutHAEMATO <- gdsc1000Tissues[ibx,,drop=F]
gdsc1000RNAseq <- gdsc1000RNAseq_ALL[rownames(gdsc1000Tissues_WithoutHAEMATO),]
GDSC1000_tissues <- data.frame(table(gdsc1000Tissues[rownames(gdsc1000RNAseq),]),stringsAsFactors = F)
gdsc1000RS <- gdsc1000RNAseq[,finalSetOfGenes]
gdsc1000RS_binary <- getBinaryValues(gdsc1000RS,cutoffs_ccle2[finalSetOfGenes])


AUC_GDSC1000 <- summarizeSensitivityProfiles(GDSC1000,sensitivity.measure = "aac_recomputed",fill.missing = F)


drugs_common <- intersect(rownames(AUC_GDSC1000),rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation))


ibx <- which(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[drugs_common,"mCI"]>mciTh)

final_drugs_common <- rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[drugs_common,])[ibx]





hist(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[final_drugs_common,"mCI.balanced"],breaks = 15)

results_final_common <- lapply(final_drugs_common, function(x,output,experssion,AUCs_Drugs){
  
  AUCs_Drug <- AUCs_Drugs[x,]
  AUCs_Drug <- AUCs_Drug[!is.na(AUCs_Drug)]
  
  common <- intersect(names(AUCs_Drug),rownames(experssion))
  
  
  final_results <- lapply(1:length(output[[x]][[3]]), function(y){
    print(y)
    SolMat <- output[[x]][[3]][[y]][[2]]
    K <- as.numeric(output[[x]][[3]][[y]][[3]])
    M <- as.numeric(output[[x]][[3]][[y]][[4]])
    print((match(colnames(SolMat),colnames(experssion))))
    results <- getPredictions_LOBICO_outer(SolMat,cbind("drug_AUC"=AUCs_Drug[common], experssion[common,colnames(SolMat)]),K,M)
  })
  
  final_results <- do.call(cbind,final_results)
  
  final_results <- cbind(final_results[,seq(1,length(output[[x]][[3]])*2,2)],final_results[,2])
  
  finalPredictions <- apply(final_results[,1:length(output[[x]][[3]]),drop=F], 1, function(x){
    round(mean(x))
  })
  
  return(cbind("Vote"=finalPredictions,"Obs"=final_results[,dim(final_results)[2]]))
},output=output_ensemble_CTRPv2_5CV3NS_LUNG_final,experssion=gdsc1000RS_binary[,finalSetOfGenes],AUCs_Drugs=AUC_GDSC1000)




names(results_final_common) <- final_drugs_common

library(wCI)
results_final_common_evaluation <- lapply(results_final_common, function(x){
  
  output <- x
  
  expressedCellLines <- rownames(output)[output[,"Vote"]==1]
  notExpressedCellLines <- rownames(output)[output[,"Vote"]==0]
  
  #integrCindex <- Hmisc::rcorr.cens(S=output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
  CI <- wCI::paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                      ,delta.pred = 0,delta.obs = 0,outx = F)
  
  mCI <- wCI::paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                       ,delta.pred = 0,delta.obs = 0.2,outx = F)
  
  results <- c("CI"=as.numeric(CI["cindex"]),"CI.pval"=as.numeric(CI["p.value"]),"mCI"=as.numeric(mCI["cindex"]),"mCI.pval"=as.numeric(mCI["p.value"]))
  return(results)
  
  
})


results_final_common_evaluation <- do.call(rbind,results_final_common_evaluation)



pdf("Plots/Validating_on_GDSC1000_LUNG.pdf")
plot(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"]
     ,results_final_common_evaluation[,"mCI"],xlab = "CTRPv2 - training [mCI]",ylab = "GDSC1000 - testing [mCI]",pch=19,main="Validating on GDSC1000"
     ,ylim=range(c(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"],results_final_common_evaluation[,"mCI"]),na.rm = T)
     ,xlim=range(c(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"],results_final_common_evaluation[,"mCI"]),na.rm = T))
#     ,ylim=c(min(results_final_common_evaluation[,"mCI"],na.rm = T),max(results_final_common_evaluation[,"mCI"],na.rm = T))
#     ,xlim=c(min(results_final_common_evaluation[,"mCI"],na.rm = T),max(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"],na.rm = T))
#     )


ibx <- output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"] >= mciTh
points(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation)[ibx],"mCI"]
       ,results_final_common_evaluation[ibx,"mCI"],col="red",pch=19)

iax <- output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"] >= mciTh & results_final_common_evaluation[,"mCI"] >= mciTh
iax[is.na(iax)] <- F
#legend("topleft",legend = c("mCI.balanced training >= 0.6"),fill = c("red"),bty = "n")
points(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation)[ibx & iax],"mCI"]
       ,results_final_common_evaluation[ibx & iax,"mCI"],col="blue",pch=19)
dev.off()

sum(iax)/sum(ibx)

idx <- unlist(lapply(rownames(results_final_common_evaluation), function(x){
  ((sum(AUC_GDSC1000[x,]>0.2,na.rm = T)/sum(AUC_GDSC1000[x,]>=0,na.rm = T)) <0.1)# |  ((sum(AUC_GDSC1000[x,]<0.2,na.rm = T)/sum(AUC_GDSC1000[x,]>=0,na.rm = T)) <0.1) 
}))


sum(!idx & iax)/sum(!idx & ibx)


points(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation)[idx],"mCI"]
       ,results_final_common_evaluation[idx,"mCI"],pch=16,col="skyblue")

text(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"]
     ,results_final_common_evaluation[,"mCI"],labels = rownames(results_final_common_evaluation),pos = 1)

legend("topleft",legend = c("CI.balanced training >= 0.6","AUC range is very skewed"),fill = c("red","skyblue"),bty = "n")



data <- data.frame("CTRPv2"=output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"]
                   ,"GDSC1000"=results_final_common_evaluation[,"mCI"],"drug"=rownames(results_final_common_evaluation),stringsAsFactors=F)

pdf("Plots/Validating_on_GDSC1000_LUNG.pdf",height = 6,width = 14)
dataF <- t(as.matrix(data[,c(1,2)]))
dataF <- dataF[,names(sort(dataF["CTRPv2",]))]
dataF <- dataF - 0.5
par(mai=c(1.5,1,1,1))
barplot(dataF,beside = T,xpd = F,col=rev(c("#e9a3c9","#a1d76a"))
        ,border = NA,ylab="mCI",axes=FALSE
        ,ylim = c(min(dataF,na.rm = T)-0.01,max(dataF,na.rm = T)+0.01),las=2)


axis(side=2,at=seq(0,1,0.1)-0.5,labels=(seq(0,1,0.1)))
abline(h=0)
legend("topleft",legend = c("CTRPv2","GDSC1000")
       ,fill = rev(c("#e9a3c9","#a1d76a")),bty = "n",border = F)

abline(h=0.1,col="red",lty=2)
dev.off()


pdf("Plots/topHits_GDSC1000_barPlots_LUNG.pdf",width = 12,height = 6)
lapply(rownames(results_final_common_evaluation)[ibx], function(x){
  results <- results_final_common[[x]]
  results <- results[order(results[,2],decreasing = T),]
  results_tmp <- results
  results_tmp[,2] <- results_tmp[,2]-0.2
  barplot(results_tmp[,2],col = ifelse(results_tmp[,1]==1,"#75d2ff","#f24720"),border = NA,axes=FALSE,names.arg = NA,ylab = "AAC")
  axis(side=2,at=seq(0,1,0.1)-0.2,labels=(seq(0,1,0.1)))
  legend("topright",legend = c("Sensitive","Resistant"),fill = c("#75d2ff","#f24720"),bty="n")
  title(main = paste(x,"[GDSC1000]\n","mCI:",sprintf("%.2f",as.numeric(results_final_common_evaluation[x,"mCI"]))))
})

dev.off()





#############################################
# validate on GCSI
library(Biobase)

gCSI <- readRDS("~/Projects/PSets/final_psets/gCSI_2018.rds")



gCSIRnaSeq_ALL <- t(exprs(summarizeMolecularProfiles(gCSI,mDataType = "rnaseq",fill.missing = F)))
gCSITissues <- gCSI@cell[rownames(gCSIRnaSeq_ALL),"tissueid",drop=F]
ibx <- which(gCSITissues$tissueid=="lung")
gCSITissues_WithoutHAEMATO <- gCSITissues[ibx,,drop=F]
gCSIRnaSeq <- gCSIRnaSeq_ALL[rownames(gCSITissues_WithoutHAEMATO),]

gCSI_tissues <- data.frame(table(gCSITissues[rownames(gCSIRnaSeq),]),stringsAsFactors = F)
gcsi <- gCSIRnaSeq[,finalSetOfGenes]
gcsi_binary <- getBinaryValues(gcsi,cutoffs_ccle2[finalSetOfGenes])


AUC_gCSI <- summarizeSensitivityProfiles(gCSI,sensitivity.measure = "aac_recomputed",fill.missing = F)


drugs_common <- intersect(rownames(AUC_gCSI),rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation))


ibx <- which(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[drugs_common,"mCI"]>mciTh)

final_drugs_common <- rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[drugs_common,])[ibx]






results_final_common <- lapply(final_drugs_common, function(x,output,experssion,AUCs_Drugs){
  
  AUCs_Drug <- AUCs_Drugs[x,]
  AUCs_Drug <- AUCs_Drug[!is.na(AUCs_Drug)]
  
  common <- intersect(names(AUCs_Drug),rownames(experssion))
  
  
  final_results <- lapply(1:length(output[[x]][[3]]), function(y){
    SolMat <- output[[x]][[3]][[y]][[2]]
    K <- as.numeric(output[[x]][[3]][[y]][[3]])
    M <- as.numeric(output[[x]][[3]][[y]][[4]])
    
    results <- getPredictions_LOBICO_outer(SolMat,cbind("drug_AUC"=AUCs_Drug[common], experssion[common,colnames(SolMat)]),K,M)
  })
  
  final_results <- do.call(cbind,final_results)
  
  final_results <- cbind(final_results[,seq(1,length(output[[x]][[3]])*2,2)],final_results[,2])
  
  finalPredictions <- apply(final_results[,1:length(output[[x]][[3]]),drop=F], 1, function(x){
    round(mean(x))
  })
  
  return(cbind("Vote"=finalPredictions,"Obs"=final_results[,dim(final_results)[2]]))
},output=output_ensemble_CTRPv2_5CV3NS_LUNG_final,experssion=gcsi_binary[,finalSetOfGenes],AUCs_Drugs=AUC_gCSI)




names(results_final_common) <- final_drugs_common


results_final_common_evaluation <- lapply(results_final_common, function(x){
  
  output <- x
  
  expressedCellLines <- rownames(output)[output[,"Vote"]==1]
  notExpressedCellLines <- rownames(output)[output[,"Vote"]==0]
  
  #integrCindex <- Hmisc::rcorr.cens(S=output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
  CI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                 ,delta.pred = 0,delta.obs = 0,outx = F)
  
  mCI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                  ,delta.pred = 0,delta.obs = 0.2,outx = F)
  
  results <- c("CI"=as.numeric(CI["cindex"]),"CI.pval"=as.numeric(CI["p.value"]),"mCI"=as.numeric(mCI["cindex"]),"mCI.pval"=as.numeric(mCI["p.value"]))
  return(results)
  
  
})


results_final_common_evaluation <- do.call(rbind,results_final_common_evaluation)


cor.test(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"]
         ,results_final_common_evaluation[,"mCI"])



pdf("Plots/Validating_on_gCSI_LUNG.pdf")
plot(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"]
     ,results_final_common_evaluation[,"mCI"],xlab = "CTRPv2 - training [mCI]",ylab = "gCSI - testing [mCI]",pch=19,main="Vaalidating on gCSI"
     ,ylim=range(c(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"],results_final_common_evaluation[,"mCI"]),na.rm = T)
     ,xlim=range(c(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"],results_final_common_evaluation[,"mCI"]),na.rm = T))


ibx <- output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"] >= mciTh
points(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation)[ibx],"mCI"]
       ,results_final_common_evaluation[ibx,"mCI"],col="red",pch=19)

text(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation)[ibx],"mCI"]
     ,results_final_common_evaluation[ibx,"mCI"],labels = rownames(results_final_common_evaluation)[ibx],pos = 2)



iax <- output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"] >= mciTh & results_final_common_evaluation[,"mCI"] >= mciTh
iax[is.na(iax)] <- F
#legend("topleft",legend = c("mCI.balanced training >= 0.6"),fill = c("red"),bty = "n")
points(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation)[ibx & iax],"mCI"]
       ,results_final_common_evaluation[ibx & iax,"mCI"],col="blue",pch=19)


idx <- unlist(lapply(rownames(results_final_common_evaluation), function(x){
  ((sum(AUC_gCSI[x,]>0.2,na.rm = T)/sum(AUC_gCSI[x,]>=0,na.rm = T)) <0.1)# |  ((sum(AUC_GDSC1000[x,]<0.2,na.rm = T)/sum(AUC_GDSC1000[x,]>=0,na.rm = T)) <0.1) 
}))


sum(!idx & iax)/sum(!idx & ibx)

sum(iax)/sum(ibx)
dev.off()


data <- data.frame("CTRPv2"=output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[rownames(results_final_common_evaluation),"mCI"]
                   ,"gCSI"=results_final_common_evaluation[,"mCI"],"drug"=rownames(results_final_common_evaluation),stringsAsFactors=F)


pdf("Plots/Validating_on_gCSI_LUNG.pdf",height = 6,width = 14)
dataF <- t(as.matrix(data[,c(1,2)]))
dataF <- dataF[,names(sort(dataF["CTRPv2",]))]
dataF <- dataF - 0.5
p <- barplot(dataF,beside = T,xpd = F,col=rev(c("#e9a3c9","#a1d76a"))
             ,border = NA,ylab="mCI",axes=FALSE#,las=2
             ,ylim = c(min(dataF,na.rm = T)-0.01,max(dataF,na.rm = T)+0.05))


axis(side=2,at=seq(0,1,0.1)-0.5,labels=(seq(0,1,0.1)))
abline(h=0)
#legend("topleft",legend = c("CTRPv2","gCSI","FDR < 0.05")
#       ,fill = rev(c(NA,"#e9a3c9","#a1d76a")),pch = c(NA,NA,8),bty = "n",border = F)
legend("topleft",legend = c("CTRPv2","gCSI")
       ,fill = rev(c("#e9a3c9","#a1d76a")),bty = "n",border = F)
abline(h=0.1,col="red",lty=2)


#points(colMeans(p),apply(dataF, 2, max)+0.015,pch=ifelse(p.adjust(results_final_common_evaluation[colnames(dataF),"mCI.pval"],method = "fdr")<0.05,8,NA))

dev.off()


pdf("Plots/topHits_gCSI_barPlots_LUNG.pdf",width = 12,height = 6)
lapply(rownames(results_final_common_evaluation)[ibx & iax], function(x){
  results <- results_final_common[[x]]
  results <- results[order(results[,2],decreasing = T),]
  results_tmp <- results
  results_tmp[,2] <- results_tmp[,2]-0.2
  barplot(results_tmp[,2],col = ifelse(results_tmp[,1]==1,"#75d2ff","#f24720"),border = NA,axes=FALSE,names.arg = NA,ylab = "AAC")
  axis(side=2,at=seq(0,1,0.1)-0.2,labels=(seq(0,1,0.1)))
  legend("topright",legend = c("Sensitive","Resistant"),fill = c("#75d2ff","#f24720"),bty="n")
  title(main = paste(x,"[GDSC1000]\n","mCI:",sprintf("%.2f",as.numeric(results_final_common_evaluation[x,"mCI"]))))
})

dev.off()


################################################
# validate on patients: NeoAlto

neoAlto <- readRDS("NeoALTTO_Exp_Mappping_Meta.rds")

normalized_neoAltoRS <- neoAlto$TPMMat_GeneLevel

library(preprocessCore)
CCLE_space <- normalize.quantiles.determine.target(t(ccle2)) # Getting Quantiles for CCLE expression matrix
normalized_neoAltoRS <- normalize.quantiles.use.target(normalized_neoAltoRS, target=CCLE_space) # Quantile normalization of BTSC expression matrix
colnames(normalized_neoAltoRS) <- colnames(neoAlto$TPMMat_GeneLevel)
rownames(normalized_neoAltoRS) <- rownames(neoAlto$TPMMat_GeneLevel)

neoAltoRS_binary <- getBinaryValues(t(normalized_neoAltoRS),cutoffs_ccle2)
neoAltoRS <- t(neoAlto$TPMMat_GeneLevel)

mappingRS <- neoAlto$Annotation_mapping$`Patient screening no.`[match(neoAlto$Annotation_mapping$`FASTQ files`,rownames(neoAltoRS_binary))]
names(mappingRS) <- neoAlto$Annotation_mapping$`FASTQ files`[match(neoAlto$Annotation_mapping$`FASTQ files`,rownames(neoAltoRS_binary))]
commonNeoAlto <- intersect(rownames(neoAltoRS_binary),names(mappingRS))

neoAltoRS_binary <- neoAltoRS_binary[commonNeoAlto,]
neoAltoRS <- neoAltoRS[commonNeoAlto,]
mappingRS <- mappingRS[commonNeoAlto]

patientsTreated <- unique(mappingRS)

neoAlto_response <- as.numeric(neoAlto$ClinicalInfo$pCR)
neoAlto_response <- neoAlto_response-1
neoAlto_response <- as.matrix(neoAlto_response)
rownames(neoAlto_response) <- neoAlto$ClinicalInfo$PatientID
neoAlto_response <- neoAlto_response[as.character(patientsTreated),,drop=F]


neoAlto_response_SRV <- as.numeric(neoAlto$ClinicalInfo$efsYears)
neoAlto_response_SRV <- as.matrix(neoAlto_response_SRV)
rownames(neoAlto_response_SRV) <- neoAlto$ClinicalInfo$PatientID
neoAlto_response_SRV <- neoAlto_response_SRV[as.character(patientsTreated),,drop=F]

neoAlto_response_SRV_ev <- as.numeric(neoAlto$ClinicalInfo$efs)
neoAlto_response_SRV_ev <- as.matrix(neoAlto_response_SRV_ev)
rownames(neoAlto_response_SRV_ev) <- neoAlto$ClinicalInfo$PatientID
neoAlto_response_SRV_ev <- neoAlto_response_SRV_ev[as.character(patientsTreated),,drop=F]



neoAlto_response_SRV_O <- as.numeric(neoAlto$ClinicalInfo$osYears)
neoAlto_response_SRV_O <- as.matrix(neoAlto_response_SRV_O)
rownames(neoAlto_response_SRV_O) <- neoAlto$ClinicalInfo$PatientID
neoAlto_response_SRV_O <- neoAlto_response_SRV_O[as.character(patientsTreated),,drop=F]

neoAlto_response_SRV_O_ev <- as.numeric(neoAlto$ClinicalInfo$os)
neoAlto_response_SRV_O_ev <- as.matrix(neoAlto_response_SRV_O_ev)
rownames(neoAlto_response_SRV_O_ev) <- neoAlto$ClinicalInfo$PatientID
neoAlto_response_SRV_O_ev <- neoAlto_response_SRV_O_ev[as.character(patientsTreated),,drop=F]


#ibx <- match(rownames(neoAlto_response),mappingRS)

#mappingRS_final <- mappingRS[ibx]
#neoAltoRS_binary_final <- neoAltoRS_binary[names(mappingRS_final),]

drug <- "LAPATINIB ALONE"

armDrug <- neoAlto$ClinicalInfo[neoAlto$ClinicalInfo$randarm == drug,c("PatientID","randarm")]
commonSamples <- intersect(armDrug$PatientID,rownames(neoAlto_response))
neoAlto_response_drug <- neoAlto_response[commonSamples,,drop=F]


ibx <- match(rownames(neoAlto_response_drug),mappingRS)

mappingRS_final <- mappingRS[ibx]
neoAltoRS_binary_final <- neoAltoRS_binary[names(mappingRS_final),]
neoAltoRS_final <- neoAltoRS[names(mappingRS_final),]

x <- "Lapatinib"

output_ensemble_CTRPv2_5CV3NS_LUNG_final <- output_ensemble_CTRPv2_5CV3NS_LUNG
#output_ensemble_CTRPv2_5CV3NS_LUNG_final <- output_CTRPv2_All_BimodalGenomeWide_lapatinib

final_results <- lapply(1:length(output_ensemble_CTRPv2_5CV3NS_LUNG_final[[x]][[3]]), function(y){
  SolMat <- output_ensemble_CTRPv2_5CV3NS_LUNG_final[[x]][[3]][[y]][[2]]
  K <- as.numeric(output_ensemble_CTRPv2_5CV3NS_LUNG_final[[x]][[3]][[y]][[3]])
  M <- as.numeric(output_ensemble_CTRPv2_5CV3NS_LUNG_final[[x]][[3]][[y]][[4]])
  
  results <- getPredictions_LOBICO_outer(SolMat,cbind("drug_AUC"=neoAlto_response_drug[,,drop=T], neoAltoRS_binary_final[,colnames(SolMat)]),K,M)
})

final_results <- do.call(cbind,final_results)

final_results <- cbind(final_results[,seq(1,length(output_ensemble_CTRPv2_5CV3NS_LUNG_final[[x]][[3]])*2,2)],final_results[,2])

finalPredictions <- apply(final_results[,1:length(output_ensemble_CTRPv2_5CV3NS_LUNG_final[[x]][[3]]),drop=F], 1, function(x){
  round(mean(x))
})

final_results_summary <- cbind("Vote"=finalPredictions,"Obs"=final_results[,dim(final_results)[2]])

output <- final_results_summary

expressedCellLines <- rownames(output)[output[,"Vote"]==1]
notExpressedCellLines <- rownames(output)[output[,"Vote"]==0]

library(wCI)
library(e1071)
#integrCindex <- Hmisc::rcorr.cens(S=output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
CI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),"Obs"],output[c(notExpressedCellLines,expressedCellLines),1]
                               ,delta.pred = 0,delta.obs = 0,outx = F)

mCI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),"Obs"],output[c(notExpressedCellLines,expressedCellLines),1]
                                ,delta.pred = 0,delta.obs = 0,outx = F)$cindex

boxplot(list(output[c(notExpressedCellLines),"Obs"],output[c(expressedCellLines),"Obs"]))

results <- c("CI"=as.numeric(CI["cindex"]),"CI.pval"=as.numeric(CI["p.value"]),"mCI"=as.numeric(mCI["cindex"]),"mCI.pval"=as.numeric(mCI["p.value"]))
mccr(pred = as.numeric(output[,"Vote"]),act = as.numeric(output[,"Obs"]))
output_caret <- output
output_caret[output_caret==0] <- 2
C <-caret::confusionMatrix(as.factor(output_caret[,"Obs"]), as.factor(output_caret[,"Vote"]))

fisher.test(y=as.numeric(output[,"Vote"]),x = as.numeric(output[,"Obs"]))


library(survcomp)
D.index(x=output[,"Vote"], surv.time=neoAlto_response_SRV[rownames(output),], surv.event=neoAlto_response_SRV_ev[rownames(output),],na.rm = TRUE)$p.val
D.index(x=output[,"Vote"], surv.time=neoAlto_response_SRV_O[rownames(output),], surv.event=neoAlto_response_SRV_O_ev[rownames(output),],na.rm = TRUE)$p.val

low <- neoAlto_response_SRV[neoAlto_response_SRV[,1]<=1.5,,drop=F]
high <- neoAlto_response_SRV[neoAlto_response_SRV[,1]>=6.5,,drop=F]

D.index(x=output[intersect(rownames(output),c(rownames(low),rownames(high))),"Vote"], surv.time=neoAlto_response_SRV[intersect(rownames(output),c(rownames(low),rownames(high))),], surv.event=neoAlto_response_SRV_ev[intersect(rownames(output),c(rownames(low),rownames(high))),],na.rm = TRUE)$p.val
D.index(x=output[intersect(rownames(output),c(rownames(low),rownames(high))),"Vote"], surv.time=neoAlto_response_SRV_O[intersect(rownames(output),c(rownames(low),rownames(high))),], surv.event=neoAlto_response_SRV_O_ev[intersect(rownames(output),c(rownames(low),rownames(high))),],na.rm = TRUE)$p.val

ERBB2_neoAlto <- neoAltoRS_final[,"ENSG00000141736.13"]

D.index(x=ERBB2_neoAlto, surv.time=neoAlto_response_SRV[rownames(output),], surv.event=neoAlto_response_SRV_ev[rownames(output),],na.rm = TRUE)$p.val


##########################
##########################

# comparison with global features
common <- intersect(rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation)
                    ,rownames(output_ensemble_CTRPv2_5CV3NS_evaluation))
tmp_data <- as.data.frame(cbind("Global"=output_ensemble_CTRPv2_5CV3NS_evaluation[common,"mCI"],"Lung_specific"=output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[common,"mCI"]))


tmp_data$type <- ifelse(tmp_data$Global > 0.6 & tmp_data$Lung_specific > 0.6, "Both"
                        ,ifelse(tmp_data$Global > 0.6,"Global"
                                ,ifelse(tmp_data$Lung_specific > 0.6,"Lung_specific","None")))


table(tmp_data$type)







