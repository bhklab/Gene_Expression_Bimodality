library(Biobase)
library(PharmacoGx)
library(ggplot2)
library(GSA) 
library(piano)
library(gsubfn)
library(wCI)


setwd("~/Projects/Bimodality_paper/")

source("src/getBiModalScore.R")
source("src/getPredictions_LOBICO.R")


# Flags to set whether to run long functions
perform_CL_bimodality <- F # flag for recomputing gene expression bimodality on cell lines.
perform_MCC_bimodalGenesSimilarity <- F # flag for recomputing bimodal genes similarity using MCC metric
perform_CTRPv2_training_BimodalGenes <- F # flag for retraining models based on CTRPv2 bimodal genes - RNAseq


######################################
######################################
# Data for cell lines bimodality

CCLE <- readRDS("PSets/CCLE.rds")
genes_mappings <- featureInfo(CCLE,"rnaseq")


############
# Processsing Data

ccleRnaSeq_ALL <- t(exprs(summarizeMolecularProfiles(CCLE,mDataType = "rnaseq",fill.missing = F)))
ccleTissues <- CCLE@cell[rownames(ccleRnaSeq_ALL),"tissueid",drop=F]
ibx <- which(ccleTissues$tissueid!="haematopoietic_and_lymphoid_tissue")
ccleTissues_WithoutHAEMATO <- ccleTissues[ibx,,drop=F]
ccleRnaSeq <- ccleRnaSeq_ALL[rownames(ccleTissues_WithoutHAEMATO),]
CCLE_tissues <- data.frame(table(ccleTissues[rownames(ccleRnaSeq),]),stringsAsFactors = F)




if(perform_CL_bimodality){

  library(parallel)
  cl <- makeCluster(3)
  
  BiModalScores_ccleRnaSeq_ALL <- parApply(cl = cl,FUN = function(x){source("getBiModalScore.R");getBiModalScore_Updated(x)},MARGIN = 2,X = ccleRnaSeq)
  BiModalScores_ccleRnaSeq <- unlist(lapply(BiModalScores_ccleRnaSeq_ALL, "[[",1))
  BiModalScores_ccleRnaSeq <- BiModalScores_ccleRnaSeq[order(BiModalScores_ccleRnaSeq,decreasing = T)]
  
  stopCluster(cl)
}else{
  BiModalScores_ccleRnaSeq <-readRDS("data/BiModalScores_ccleRnaSeq.rda")
  BiModalScores_ccleRnaSeq_ALL <-readRDS("data/BiModalScores_ccleRnaSeq_ALL.rda")
}

protCodingGenes <- genes_mappings[names(BiModalScores_ccleRnaSeq),"gene_type"]
names(protCodingGenes) <- names(BiModalScores_ccleRnaSeq)

protCodingGenes <- protCodingGenes[protCodingGenes=="protein_coding"]
protCodingGenes <- names(protCodingGenes)


percentile <- 0.80
ccle_RS_th <- quantile(BiModalScores_ccleRnaSeq[protCodingGenes],probs = percentile)


#####################
# Fig 1 - b
#####################
dens=density(BiModalScores_ccleRnaSeq[protCodingGenes])
pdf("plots/Fig_1_b_big.pdf")
hist(BiModalScores_ccleRnaSeq[protCodingGenes],main = "Dist. of bimodality score in CCLE",xlab = "BI score",breaks = 50,col="skyblue",border = NA)
dev.off()

# example of one of the top bimodal genes
ibx <- which(genes_mappings$Symbol == "RAB25")

pdf("plots/Fig_1_b_corner.pdf",height = 4,width = 6)
dens=density(ccleRnaSeq[,ibx])
hist(ccleRnaSeq[,ibx],main = featureInfo(CCLE,"rnaseq")[ibx,"Symbol"],xlab = "Expression",breaks = 20,col="skyblue",probability = T,border = NA)
lines(dens$x,length(ccleRnaSeq[,ibx])*dens$y/600,type="l",col="darkblue")
dev.off()


######################################
# Data for TCGA bimodality


#####################
# Fig 1 - c
#####################

TCGA_bimodal_genes <- readRDS("data/TCGA_bimodality_scores.rda")

TCGA_threshold <- quantile(TCGA_bimodal_genes[protCodingGenes],probs = 0.8)
ccle_Threshold <- quantile(BiModalScores_ccleRnaSeq[protCodingGenes],probs = 0.8)

ibx <- which(BiModalScores_ccleRnaSeq[names(TCGA_bimodal_genes[protCodingGenes])] >= ccle_Threshold & TCGA_bimodal_genes[protCodingGenes] >= TCGA_threshold  )
df1 <- data.frame(data.frame("TCGA"=TCGA_bimodal_genes[protCodingGenes],"CCLE"=BiModalScores_ccleRnaSeq[protCodingGenes],"group"=(TCGA_bimodal_genes[protCodingGenes]>TCGA_threshold & BiModalScores_ccleRnaSeq[protCodingGenes]>ccle_Threshold)))

finalSetOfGenes <- names(ibx)
names(finalSetOfGenes) <- genes_mappings[finalSetOfGenes,"Symbol"]


myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
myColor_scale_fill <- scale_fill_gradientn(colours = myColor)

pdf("plots/Fig_1_c.pdf",width = 11,height = 10)
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
dev.off()

cor.test(df1$TCGA,df1$CCLE,method = "p")

################################
################################
################################

# Binarizing the gene expression data based on the bimodaly features.
ccleRnaSeq_final <- ccleRnaSeq[,finalSetOfGenes]
colnames(ccleRnaSeq_final) <- finalSetOfGenes

cutoffs_ccle <- unlist(lapply(colnames(ccleRnaSeq_final), function(x){
  val <- BiModalScores_ccleRnaSeq_ALL[[x]][["mix"]][["m.step"]][["mu"]]
  return(mean(val))
}))
names(cutoffs_ccle) <- colnames(ccleRnaSeq_final)

ccleRnaSeq_final_binary <- getBinaryValues(ccleRnaSeq_final,cutoffs_ccle)



#####################
# Pathway enrichment analysis on the top bimodal genes

library(GSA) 
library(piano)

gsc1 <- loadGSC("data/c2.cp.reactome.v6.2.symbols.gmt")


finalSetOfGenes_tmp <- finalSetOfGenes[!duplicated(names(finalSetOfGenes))]

gsea_out <- piano::runGSAhyper(genes = names(finalSetOfGenes_tmp), gsc=gsc1,universe = unique(genes_mappings[protCodingGenes,"Symbol"]))
gsea_out_final <- as.data.frame(gsea_out$resTab)
gsea_out_final$fdr <- p.adjust(gsea_out_final$`p-value`,method = "fdr")
gsea_out_final <- gsea_out_final[order(gsea_out_final$`p-value`),]

gsea_out_final1 <- gsea_out_final[gsea_out_final$fdr<0.05,]

gsea_out_final1 <- gsea_out_final1[order(gsea_out_final1$fdr,decreasing = T),]

#####################
# Fig S1 - a
#####################

pdf("plots/Fig_S1_a.pdf",width = 11,height = 8)
par(mai=c(1,7,1,1))
#barplot(-log10(gsea_out_final1$fdr),horiz = T,names.arg = gsub("REACTOME_","",rownames(gsea_out_final1)),las=2,xlab = "-log10(FDR)")
barplot(-log10(gsea_out_final1$fdr),horiz = T,names.arg =gsub("REACTOME","",gsub("_"," ",rownames(gsea_out_final1))),las=2,xlab = "-log10(FDR)")
dev.off()
########
# performing the Matthew correlation coefficient to measure the similarity between the top bimodal genes

if(perform_MCC_bimodalGenesSimilarity){
mccCCLE <- matrix(nrow = dim(ccleRnaSeq_final_binary)[2],ncol = dim(ccleRnaSeq_final_binary)[2])
for (i in 1:(dim(ccleRnaSeq_final_binary)[2]-1)) {
  for (j in (i+1):dim(ccleRnaSeq_final_binary)[2]) {
    tmpEst <- mccr(as.factor(ccleRnaSeq_final_binary[,i]),as.factor(ccleRnaSeq_final_binary[,j]))
    mccCCLE[i,j] <- tmpEst
    mccCCLE[j,i] <- tmpEst
  }
}
}else{
  mccCCLE <- readRDS("data/mccCCLE.rda")
}

print(paste("IQR for similarity between the top bimodal genes (MCC):",sprintf("%0.2f",IQR(mccCCLE,na.rm = T))))


#####################
# Fig S1 - b
#####################
dens=density(mccCCLE[!is.na(mccCCLE)])
# Plot y-values scaled by number of observations against x values
pdf("plots/Fig_S1_b.pdf")
hist(mccCCLE[!is.na(mccCCLE)],col="skyblue",breaks = 500,main = "Bimodal genes pairwise-correlation",xlab = "MCC",border = NA)
lines(dens$x,length(mccCCLE[!is.na(mccCCLE)])*dens$y/200,type="l",col="darkblue")
dev.off()


#####################
#####################
# Training on CTRPv2 data set

CTRPv2 <- readRDS("PSets/CTRPv2_clean.rds")

AAC_CTRPv2 <- summarizeSensitivityProfiles(CTRPv2,sensitivity.measure = "aac_recomputed",fill.missing = F)

# filtering drugs to only those that have at least 10% of tested cell lines to "sensitive", i.e. AAC < 0.2
y <- 0.1
keptDrugs <- apply(AAC_CTRPv2, 1, function(x){
  x <- x[!is.na(x)]
  ((sum(x>0.2,na.rm = T)/sum(x>=0,na.rm = T)) >y) 
})




if(perform_CTRPv2_training_BimodalGenes){
  

featuresType <- "output_mCI_CCLE_TCGA_75_F"

listOfdrugsFiles <- dir(paste("AzureRun/",featuresType,"/",sep = ""))

output_ensemble_CTRPv2_5CV3NS <- lapply(listOfdrugsFiles, function(x){
  print(x)
  output <- readRDS(paste("AzureRun/",featuresType,"/",x,sep = ""))
  return(output)
})

drugsProcessed <- unlist(lapply(strsplit(listOfdrugsFiles,"_CTRPv2_clean"),"[[",1))
drugsProcessed <- gsub("%2C",",",drugsProcessed)
drugsProcessed <- gsub("%5B","[",drugsProcessed)
drugsProcessed <- gsub("%5D","]",drugsProcessed)
drugsProcessed <- gsub("%20"," ",drugsProcessed)
drugsProcessed <- gsub("%7B","{",drugsProcessed)
drugsProcessed <- gsub("%7D","}",drugsProcessed)

names(output_ensemble_CTRPv2_5CV3NS) <- drugsProcessed
}else{
  output_ensemble_CTRPv2_5CV3NS <- readRDS("data/output_ensemble_CTRPv2_5CV3NS.rda")
}

output_ensemble_CTRPv2_5CV3NS_final <- output_ensemble_CTRPv2_5CV3NS[intersect(names(keptDrugs)[keptDrugs],names(output_ensemble_CTRPv2_5CV3NS))]
mciTh = 0.6

output_ensemble_CTRPv2_5CV3NS_evaluation <- lapply(output_ensemble_CTRPv2_5CV3NS_final,function(x){
  return(x[[1]])
})

output_ensemble_CTRPv2_5CV3NS_CIs <- unlist(lapply(1:length(output_ensemble_CTRPv2_5CV3NS_evaluation), function(x){output_ensemble_CTRPv2_5CV3NS_evaluation[[x]][2]}))
output_ensemble_CTRPv2_5CV3NS_CIs <- as.numeric(output_ensemble_CTRPv2_5CV3NS_CIs)
names(output_ensemble_CTRPv2_5CV3NS_CIs) <- names(output_ensemble_CTRPv2_5CV3NS_evaluation)
output_ensemble_CTRPv2_5CV3NS_CIs <- output_ensemble_CTRPv2_5CV3NS_CIs[!is.na(output_ensemble_CTRPv2_5CV3NS_CIs)]

output_ensemble_CTRPv2_5CV3NS_final <- output_ensemble_CTRPv2_5CV3NS_final[names(output_ensemble_CTRPv2_5CV3NS_CIs)]

output_ensemble_CTRPv2_5CV3NS_evaluation <- output_ensemble_CTRPv2_5CV3NS_evaluation[names(output_ensemble_CTRPv2_5CV3NS_CIs)]
output_ensemble_CTRPv2_5CV3NS_evaluation <- do.call(rbind,output_ensemble_CTRPv2_5CV3NS_evaluation)
rownames(output_ensemble_CTRPv2_5CV3NS_evaluation) <- names(output_ensemble_CTRPv2_5CV3NS_CIs)
class(output_ensemble_CTRPv2_5CV3NS_evaluation) <- "numeric"

output_ensemble_CTRPv2_5CV3NS_evaluation <- output_ensemble_CTRPv2_5CV3NS_evaluation[order(output_ensemble_CTRPv2_5CV3NS_evaluation[,"mCI"],decreasing = T),,drop=F]
output_ensemble_CTRPv2_5CV3NS_evaluation <- cbind(output_ensemble_CTRPv2_5CV3NS_evaluation,"fdr"=p.adjust(output_ensemble_CTRPv2_5CV3NS_evaluation[,"Pvalue"]))
output_ensemble_CTRPv2_5CV3NS_evaluation[,"fdr"] <- p.adjust(output_ensemble_CTRPv2_5CV3NS_evaluation[,"Pvalue"],method = "bonferroni")

drugClasses <- readRDS("data/drugClasses_CTRPv2_GDSC_Nehme.rda")
output_ensemble_CTRPv2_5CV3NS_evaluation_drugClasses <- as.data.frame(output_ensemble_CTRPv2_5CV3NS_evaluation,stringsAsFactors=F)
output_ensemble_CTRPv2_5CV3NS_evaluation_drugClasses <- cbind(output_ensemble_CTRPv2_5CV3NS_evaluation_drugClasses
                                                              ,"Class"=drugClasses[rownames(output_ensemble_CTRPv2_5CV3NS_evaluation_drugClasses),"PossibleGDSC1000Class_final"],stringsAsFactors=F)

cols_DrugClasses <- c(` other` = "#004700", `ABL signaling` = "#004743", `apoptosis regulation` = "#b1ff64", 
                      `cell cycle` = "#f20040", chromatin = "#0185ff", cytoskeleton = "#707d00", 
                      `DNA replication` = "#b900a3", `EGFR signaling` = "#002562", 
                      `ERK MAPK signaling` = "#f8c500", `Genome integrity` = "#01e3a8", 
                      `IGFR signaling` = "#8c26d9", metabolism = "#01b631", mitosis = "#36000e", 
                      other = "#d7d3ff", `p53 pathway` = "#a1004d", `PI3K signaling` = "#483400", 
                      `RTK signaling` = "#00d9e5", `TOR signaling` = "#ffeea1")

output_ensemble_CTRPv2_5CV3NS_evaluation_drugClasses <- cbind(output_ensemble_CTRPv2_5CV3NS_evaluation_drugClasses,"ClassCol"=rep(NA,dim(output_ensemble_CTRPv2_5CV3NS_evaluation_drugClasses)[1]),stringsAsFactors=F)

for (i in 1:dim(output_ensemble_CTRPv2_5CV3NS_evaluation_drugClasses)[1]) {
  output_ensemble_CTRPv2_5CV3NS_evaluation_drugClasses[i,"ClassCol"] <- cols_DrugClasses[output_ensemble_CTRPv2_5CV3NS_evaluation_drugClasses[i,"Class"]]
}


sum(output_ensemble_CTRPv2_5CV3NS_evaluation[,"mCI"]>0.6 & output_ensemble_CTRPv2_5CV3NS_evaluation[,"Pvalue"]<0.05)/dim(output_ensemble_CTRPv2_5CV3NS_evaluation)[1]

#####################
# Fig 2 - a
#####################

pdf("plots/Fig_2_a.pdf",width = 12)
data1 <- sort(output_ensemble_CTRPv2_5CV3NS_evaluation[,"mCI"])
plotdata_tmp <- data.frame("CI"=data1,stringsAsFactors = F)
plotdata_tmp$drug = rownames(plotdata_tmp)
plotdata_tmp$drug <- factor(plotdata_tmp$drug, levels = plotdata_tmp$drug)
plotdata_tmp$CI[plotdata_tmp$CI<=0.5] = 0.499
ggplot(plotdata_tmp,aes(drug,CI,fill=CI)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_gradient(low = "gray",high = "darkblue") +
  coord_cartesian(ylim=c(0.49, 0.8)) +
  geom_hline(yintercept = 0.6,col="red",lty=2) +
  ylab("CI") + xlab("Drugs") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10)) +
  ggtitle("CTRPv2 - cross validation results [CI]") 

dev.off()





tmp <- cbind("drug"=names(data1)[data1>0.6],output_ensemble_CTRPv2_5CV3NS_evaluation_drugClasses[names(data1)[data1>0.6],c("mCI","Class","ClassCol")])
tmp <- rbind(tmp,cbind("drug"=names(data1)[data1<=0.6],output_ensemble_CTRPv2_5CV3NS_evaluation_drugClasses[names(data1)[data1<=0.6],c("mCI","Class","ClassCol")]))
tmp$type <- ifelse(tmp$mCI>0.6,0,1)
tmp$label <- ""
ibx <- which(tmp$Class==" other")
tmp$Class[ibx] <- "other"
drugClasses_levels <- c(
  "chromatin",
  "TOR signaling",
  "apoptosis regulation",
  "ABL signaling",
  "cell cycle",
  "cytoskeleton",
  "IGFR signaling",
  "EGFR signaling",
  "Genome integrity",
  "ERK MAPK signaling",
  "DNA replication" ,
  "mitosis",
  "metabolism",
  "PI3K signaling",
  "other" ,
  "p53 pathway",
  "RTK signaling"
)


source("PieDonut_mod.R")

#####################
# Fig 2 - b
#####################

pdf("plots/Fig_2_b.pdf",width = 15,height = 15)
PieDonut_mod(data = tmp,mapping = aes(pies=Class,donuts=type,label=label),piesLevels = drugClasses_levels,showlabelsDonut = F
             ,donutLabelSize = 5,pieLabelSize = 5.5,addDonutLabel = F,showRatioDonut = T
             ,showPieName = F,showDonutName = F,showRatioThreshold = 0,use.label = T)
dev.off()


###################################
###################################
###################################


topHits <- rownames(output_ensemble_CTRPv2_5CV3NS_evaluation)[output_ensemble_CTRPv2_5CV3NS_evaluation[,"mCI"] >= 0.6]

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
},output=output_ensemble_CTRPv2_5CV3NS_final,experssion=ccleRnaSeq_final_binary[,finalSetOfGenes],AUCs_Drugs=AAC_CTRPv2)


names(results_final_common) <- topHits

rules <- lapply(topHits,function(drug){
  
  
  rule <- unique(unlist(output_ensemble_CTRPv2_5CV3NS[[drug]][[2]]))
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

rules_symbols_final <- unlist(lapply(1:length(rules_symbols),function(x){
  paste("[",rules_symbols[[x]],"]",collapse = " | ")
}))

names(rules_symbols_final) <- names(rules_symbols)


#####################
# Supplementry File 1
#####################


pdf("plots/SupplementryFile1.pdf",width = 12,height = 18)

for (q in seq(1,length(topHits),3)) {
  
  par(mfrow=c(3,1))
  lapply(topHits[q:min(c(length(topHits),q+2))], function(x){
    print(x)
    results <- results_final_common[[x]]
    results <- results[order(results[,2],decreasing = T),]
    results_tmp <- results
    results_tmp[,2] <- results_tmp[,2]-0.2
    barplot(results_tmp[,2],col = ifelse(results_tmp[,1]==1,"#75d2ff","#f24720"),border = NA,axes=FALSE,names.arg = NA,ylab = "AAC")
    axis(side=2,at=seq(0,1,0.1)-0.2,labels=(seq(0,1,0.1)))
    legend("topright",legend = c("Sensitive","Resistant"),fill = c("#75d2ff","#f24720"),bty="n")
    title(main = paste(x,"[CTRPv2]\n","CI:",sprintf("%.2f",as.numeric(output_ensemble_CTRPv2_5CV3NS_final[[x]][[1]][,"mCI"])),
                       "\n",paste("[",rules_symbols[[x]],"]",collapse = " | ")))
    #  tmpTissues <- ccleTissues[rownames(results_tmp),1]
    #  ibx <- match(tmpTissues,names(colors))
    #  points(p,rep(-.25,dim(results_tmp)[1]),pch=16,col=colors[ibx],cex=0.5)
    #  legend("right",legend = names(colors),fill = colors,bty="n")
  })
}
dev.off()


#####################
# Fig 2 - c
#####################

pdf("plots/Fig_2_c.pdf",width = 12,height = 12)


  par(mfrow=c(2,1))
  lapply(topHits[c(1,7)], function(x){
    print(x)
    results <- results_final_common[[x]]
    results <- results[order(results[,2],decreasing = T),]
    results_tmp <- results
    results_tmp[,2] <- results_tmp[,2]-0.2
    barplot(results_tmp[,2],col = ifelse(results_tmp[,1]==1,"#75d2ff","#f24720"),border = NA,axes=FALSE,names.arg = NA,ylab = "AAC",cex.lab=1.3)
    axis(side=2,at=seq(0,1,0.1)-0.2,labels=(seq(0,1,0.1)),cex.axis=1.25)
    legend("topright",legend = c("Sensitive","Resistant"),fill = c("#75d2ff","#f24720"),bty="n",cex = 1.2)
    title(main = paste(x,"[CTRPv2]\n","CI:",sprintf("%.2f",as.numeric(output_ensemble_CTRPv2_5CV3NS_final[[x]][[1]][,"mCI"])),
                       "\n",paste("[",rules_symbols[[x]],"]",collapse = " | ")))
  })

dev.off()


# EGFR expresion correlation with Erlotinib Rule
tmp <- output_ensemble_CTRPv2_5CV3NS$Erlotinib[[4]][,c(4,5)]


EGFR =ccleRnaSeq[rownames(tmp),"ENSG00000146648.15"]
EGFR_bin = ifelse(EGFR>mean(EGFR),1,0)
tmp <- cbind(tmp,"EGFR"=EGFR,"EGFR_bin"=EGFR_bin)

a <- cor.test(tmp[,"Vote"],tmp[,"EGFR"],method = "p")

print(paste("EGFR expression correlation with Erlotinib rule => "
            ,"PCC:"
            ,sprintf("%0.2f",a$estimate)
            ,", P-val:"
            ,sprintf("%0.2E",a$p.value)
            ))



# ERBB2 expresion correlation with Erlotinib Rule
tmp <- output_ensemble_CTRPv2_5CV3NS$Lapatinib[[4]][,c(4,5)]


ERBB2 =ccleRnaSeq[rownames(tmp),"ENSG00000141736.13"]
ERBB2_bin = ifelse(ERBB2>mean(ERBB2),1,0)
tmp <- cbind(tmp,"ERBB2"=ERBB2,"ERBB2_bin"=ERBB2_bin)

a <- cor.test(tmp[,"Vote"],tmp[,"ERBB2"],method = "p")

print(paste("ERBB2 expression correlation with Erlotinib rule => "
            ,"PCC:"
            ,sprintf("%0.2f",a$estimate)
            ,", P-val:"
            ,sprintf("%0.2E",a$p.value)
))


#######################################
#######################################
# validating on gCSI and GDSC


gCSI <- readRDS("PSets/gCSI_2018.rds")
GDSC1000 <- readRDS("PSets/GDSCv2.rds")




AUC_gCSI <- summarizeSensitivityProfiles(gCSI,sensitivity.measure = "aac_recomputed",fill.missing = F)
AUC_GDSC1000 <- summarizeSensitivityProfiles(GDSC1000,sensitivity.measure = "aac_recomputed",fill.missing = F)

commonDrugsGDSC_gCSI <- Reduce(intersect,x = list(rownames(AUC_gCSI)
                                                  ,rownames(AUC_GDSC1000)
                                                  ,rownames(output_ensemble_CTRPv2_5CV3NS_evaluation)[output_ensemble_CTRPv2_5CV3NS_evaluation[,"mCI"]>=mciTh]))



#######################################
# validate on GCSI


gCSIRnaSeq_eSet <- summarizeMolecularProfiles(gCSI,mDataType = "rnaseq",fill.missing = F)
gCSIRnaSeq_ALL <- t(exprs(gCSIRnaSeq_eSet))
gCSITissues <- gCSI@cell[rownames(gCSIRnaSeq_ALL),"tissueid",drop=F]
ibx <- which(gCSITissues$tissueid!="haematopoietic_and_lymphoid_tissue")
gCSITissues_WithoutHAEMATO <- gCSITissues[ibx,,drop=F]
gCSIRnaSeq <- gCSIRnaSeq_ALL[rownames(gCSITissues_WithoutHAEMATO),]


gcsi <- gCSIRnaSeq[,finalSetOfGenes]
gcsi_binary <- getBinaryValues(gcsi,cutoffs_ccle[finalSetOfGenes])


drugs_common <- intersect(rownames(AUC_gCSI),rownames(output_ensemble_CTRPv2_5CV3NS_evaluation))
ibx <- which(output_ensemble_CTRPv2_5CV3NS_evaluation[drugs_common,"mCI"]>mciTh)
final_drugs_common <- rownames(output_ensemble_CTRPv2_5CV3NS_evaluation[drugs_common,])[ibx]


results_final_common_gCSI <- lapply(final_drugs_common, function(x,output,experssion,AUCs_Drugs){
  
  print(x)
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
},output=output_ensemble_CTRPv2_5CV3NS_final,experssion=gcsi_binary[,finalSetOfGenes],AUCs_Drugs=AUC_gCSI)




names(results_final_common_gCSI) <- final_drugs_common

results_final_common_evaluation_gCSI <- lapply(results_final_common_gCSI, function(x){
  
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


results_final_common_evaluation_gCSI <- do.call(rbind,results_final_common_evaluation_gCSI)


data <- data.frame("CTRPv2"=output_ensemble_CTRPv2_5CV3NS_evaluation[rownames(results_final_common_evaluation_gCSI),"mCI"]
                   ,"gCSI"=results_final_common_evaluation_gCSI[,"mCI"],"drug"=rownames(results_final_common_evaluation_gCSI),stringsAsFactors=F)


#####################
# Fig 3 - a
#####################

pdf("plots/Fig_3_a.pdf",height = 6,width = 14)
dataF <- t(as.matrix(data[,c(1,2)]))
dataF <- dataF[,names(sort(dataF["gCSI",],decreasing = T))]
dataF <- dataF - 0.5

par(mai=c(2,1,1,1))
x <- barplot(dataF,beside = T,xpd = F,col = c(rep(rev(c("#66beb0","#f1b09a")),12),rep(rev(c("#abc1ca","#f1b09a")),11)) #col=rev(c("#e9a3c9","#a1d76a"))
             ,border = NA,ylab="rCI",axes=FALSE#,las=2
             ,ylim = c(0,0.35), xaxt="n")

axis(side=2,at=seq(0,1,0.1)-0.5,labels=(seq(0,1,0.1)))
abline(h=0)
legend("topright",legend = c("CTRPv2","gCSI")
       ,fill = rev(c("#66beb0","#f1b09a")),bty = "n",border = F)
abline(h=0.1,col="red",lty=2)

text(cex=1, x=colMeans(x), y=-0.009, colnames(dataF), xpd=TRUE, srt=45,adj = 1,col = ifelse(colnames(dataF) %in% commonDrugsGDSC_gCSI,"red","black"))

dev.off()

#######################
# validate on GDSCv2

gdsc1000RNAseq_ALL <- t(exprs(summarizeMolecularProfiles(GDSC1000,mDataType = "rnaseq",fill.missing = F)))
gdsc1000Tissues <- GDSC1000@cell[rownames(gdsc1000RNAseq_ALL),"tissueid",drop=F]
ibx <- which(gdsc1000Tissues$tissueid!="haematopoietic_and_lymphoid_tissue")
gdsc1000Tissues_WithoutHAEMATO <- gdsc1000Tissues[ibx,,drop=F]
gdsc1000RNAseq <- gdsc1000RNAseq_ALL[rownames(gdsc1000Tissues_WithoutHAEMATO),]
gdsc1000RS <- gdsc1000RNAseq[,finalSetOfGenes]
gdsc1000RS_binary <- getBinaryValues(gdsc1000RS,cutoffs_ccle[finalSetOfGenes])

drugs_common <- intersect(rownames(AUC_GDSC1000),rownames(output_ensemble_CTRPv2_5CV3NS_evaluation))
ibx <- which(output_ensemble_CTRPv2_5CV3NS_evaluation[drugs_common,"mCI"]>mciTh)
final_drugs_common <- rownames(output_ensemble_CTRPv2_5CV3NS_evaluation[drugs_common,])[ibx]


results_final_common_gdsc <- lapply(final_drugs_common, function(x,output,experssion,AUCs_Drugs){
  
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
},output=output_ensemble_CTRPv2_5CV3NS_final,experssion=gdsc1000RS_binary[,finalSetOfGenes],AUCs_Drugs=AUC_GDSC1000)

names(results_final_common_gdsc) <- final_drugs_common

results_final_common_evaluation_gdsc <- lapply(results_final_common_gdsc, function(x){
  
  output <- x
  
  expressedCellLines <- rownames(output)[output[,"Vote"]==1]
  notExpressedCellLines <- rownames(output)[output[,"Vote"]==0]
  
  #integrCindex <- Hmisc::rcorr.cens(S=output[c(notExpressedCellLines,expressedCellLines),2], x = as.numeric(output[c(notExpressedCellLines,expressedCellLines),1]),outx = F )
  CI <- tryCatch(wCI::paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                      ,delta.pred = 0,delta.obs = 0,outx = F),error=function(e){
                                        return(list("cindex"=NA,"p.value"=NA))
                                      })
  
  mCI <- tryCatch(wCI::paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                       ,delta.pred = 0,delta.obs = 0.2,outx = F),error=function(e){
                                         return(list("cindex"=NA,"p.value"=NA))
                                       })
  
  results <- c("CI"=as.numeric(CI["cindex"]),"CI.pval"=as.numeric(CI["p.value"]),"mCI"=as.numeric(mCI["cindex"]),"mCI.pval"=as.numeric(mCI["p.value"]))
  return(results)
  
  
})

results_final_common_evaluation_gdsc <- do.call(rbind,results_final_common_evaluation_gdsc)


data <- data.frame("CTRPv2"=output_ensemble_CTRPv2_5CV3NS_evaluation[rownames(results_final_common_evaluation_gdsc),"mCI"]
                   ,"GDSCv2"=results_final_common_evaluation_gdsc[,"mCI"],"drug"=rownames(results_final_common_evaluation_gdsc),stringsAsFactors=F)

data[data$CTRPv2<0.5,"CTRPv2"] = 0.5
data[is.na(data$GDSCv2) | data$GDSCv2<0.5,"GDSCv2"] = rep(0.5,length(which(is.na(data$GDSCv2) | data$GDSCv2<0.5)))


#####################
# Fig 3 - b
#####################

pdf("plots/Fig_3_b.pdf",height = 6,width = 14)
dataF <- t(as.matrix(data[,c(1,2)]))
dataF <- dataF[,names(sort(dataF["GDSCv2",],decreasing = T))]
dataF <- dataF - 0.5
par(mai=c(2,1,1,1))
x <- barplot(dataF,beside = T,xpd = F, col = c(rep(rev(c("#66beb0","#f1b09a")),16),rep(rev(c("#abc1ca","#f1b09a")),12)) #  col=ifelse(dataF["GDSCv2",]>0.6,rev(c("#66beb0","#f1b09a")),rev(c("#a1c1c6","#f1b09a")))
             ,border = NA,ylab="CI",axes=FALSE
             ,ylim = c(0,0.35), xaxt="n")


axis(side=2,at=seq(0,1,0.1)-0.5,labels=(seq(0,1,0.1)))
abline(h=0)
legend("topright",legend = c("CTRPv2","GDSCv2")
       ,fill = rev(c("#66beb0","#f1b09a")),bty = "n",border = F)

abline(h=0.1,col="red",lty=2)
text(cex=1, x=colMeans(x)-.2, y=-0.009, colnames(dataF), xpd=TRUE, srt=45,adj = 1,col = ifelse(colnames(dataF) %in% commonDrugsGDSC_gCSI,"red","black"))

dev.off()


##################################
##################################
##################################
# Comparison with other data types (Mutations, CNVs, Tissues)

##################################
# Mutations

featuresType <- "output_mCI_CCLE_mutation"

listOfdrugsFiels <- dir(paste("AzureRun/",featuresType,"/",sep = ""))

output_ensemble_CTRPv2_5CV3NS_mutation <- lapply(listOfdrugsFiels, function(x){
  print(x)
  output <- tryCatch(readRDS(paste("AzureRun/",featuresType,"/",x,sep = "")),error = function(e) e)
  return(output)
})


drugsProcessed <- unlist(lapply(strsplit(listOfdrugsFiels,"_CTRPv2"),"[[",1))
drugsProcessed <- gsub("%2C",",",drugsProcessed)
drugsProcessed <- gsub("%5B","[",drugsProcessed)
drugsProcessed <- gsub("%5D","]",drugsProcessed)
drugsProcessed <- gsub("%20"," ",drugsProcessed)
drugsProcessed <- gsub("%7B","{",drugsProcessed)
drugsProcessed <- gsub("%7D","}",drugsProcessed)

names(output_ensemble_CTRPv2_5CV3NS_mutation) <- drugsProcessed

output_ensemble_CTRPv2_5CV3NS_mutation_final <- output_ensemble_CTRPv2_5CV3NS_mutation[intersect(names(keptDrugs)[keptDrugs],names(output_ensemble_CTRPv2_5CV3NS_mutation))]
mciTh = 0.6

output_ensemble_CTRPv2_5CV3NS_mutation_evaluation <- lapply(output_ensemble_CTRPv2_5CV3NS_mutation_final,function(x){
  return(x[[1]])
})


output_ensemble_CTRPv2_5CV3NS_mutation_CIs <- unlist(lapply(1:length(output_ensemble_CTRPv2_5CV3NS_mutation_evaluation), function(x){output_ensemble_CTRPv2_5CV3NS_mutation_evaluation[[x]][2]}))
output_ensemble_CTRPv2_5CV3NS_mutation_CIs <- as.numeric(output_ensemble_CTRPv2_5CV3NS_mutation_CIs)
names(output_ensemble_CTRPv2_5CV3NS_mutation_CIs) <- names(output_ensemble_CTRPv2_5CV3NS_mutation_evaluation)
output_ensemble_CTRPv2_5CV3NS_mutation_CIs <- output_ensemble_CTRPv2_5CV3NS_mutation_CIs[!is.na(output_ensemble_CTRPv2_5CV3NS_mutation_CIs)]

output_ensemble_CTRPv2_5CV3NS_mutation_final <- output_ensemble_CTRPv2_5CV3NS_mutation_final[names(output_ensemble_CTRPv2_5CV3NS_mutation_CIs)]

output_ensemble_CTRPv2_5CV3NS_mutation_evaluation <- output_ensemble_CTRPv2_5CV3NS_mutation_evaluation[names(output_ensemble_CTRPv2_5CV3NS_mutation_CIs)]
output_ensemble_CTRPv2_5CV3NS_mutation_evaluation <- do.call(rbind,output_ensemble_CTRPv2_5CV3NS_mutation_evaluation)
rownames(output_ensemble_CTRPv2_5CV3NS_mutation_evaluation) <- names(output_ensemble_CTRPv2_5CV3NS_mutation_CIs)
class(output_ensemble_CTRPv2_5CV3NS_mutation_evaluation) <- "numeric"

output_ensemble_CTRPv2_5CV3NS_mutation_evaluation <- output_ensemble_CTRPv2_5CV3NS_mutation_evaluation[order(output_ensemble_CTRPv2_5CV3NS_mutation_evaluation[,"mCI"],decreasing = T),,drop=F]

##################################
# CNVs

featuresType <- "output_mCI_CCLE_cnv_all"

listOfdrugsFiels <- dir(paste("AzureRun/",featuresType,"/",sep = ""))

output_ensemble_CTRPv2_5CV3NS_cnv_all <- lapply(listOfdrugsFiels, function(x){
  print(x)
  output <- tryCatch(readRDS(paste("AzureRun/",featuresType,"/",x,sep = "")),error = function(e) e)
  return(output)
})

drugsProcessed <- unlist(lapply(strsplit(listOfdrugsFiels,"_CTRPv2"),"[[",1))
drugsProcessed <- gsub("%2C",",",drugsProcessed)
drugsProcessed <- gsub("%5B","[",drugsProcessed)
drugsProcessed <- gsub("%5D","]",drugsProcessed)
drugsProcessed <- gsub("%20"," ",drugsProcessed)
drugsProcessed <- gsub("%7B","{",drugsProcessed)
drugsProcessed <- gsub("%7D","}",drugsProcessed)

names(output_ensemble_CTRPv2_5CV3NS_cnv_all) <- drugsProcessed

output_ensemble_CTRPv2_5CV3NS_cnv_all_final <- output_ensemble_CTRPv2_5CV3NS_cnv_all[intersect(names(keptDrugs)[keptDrugs],names(output_ensemble_CTRPv2_5CV3NS_cnv_all))]
mciTh = 0.6
output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation <- lapply(output_ensemble_CTRPv2_5CV3NS_cnv_all_final,function(x){
  return(x[[1]])
})

output_ensemble_CTRPv2_5CV3NS_cnv_all_CIs <- unlist(lapply(1:length(output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation), function(x){output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation[[x]][2]}))

output_ensemble_CTRPv2_5CV3NS_cnv_all_CIs <- as.numeric(output_ensemble_CTRPv2_5CV3NS_cnv_all_CIs)
names(output_ensemble_CTRPv2_5CV3NS_cnv_all_CIs) <- names(output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation)
output_ensemble_CTRPv2_5CV3NS_cnv_all_CIs <- output_ensemble_CTRPv2_5CV3NS_cnv_all_CIs[!is.na(output_ensemble_CTRPv2_5CV3NS_cnv_all_CIs)]

output_ensemble_CTRPv2_5CV3NS_cnv_all_final <- output_ensemble_CTRPv2_5CV3NS_cnv_all_final[names(output_ensemble_CTRPv2_5CV3NS_cnv_all_CIs)]

output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation <- output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation[names(output_ensemble_CTRPv2_5CV3NS_cnv_all_CIs)]
output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation <- do.call(rbind,output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation)
rownames(output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation) <- names(output_ensemble_CTRPv2_5CV3NS_cnv_all_CIs)
class(output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation) <- "numeric"

output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation <- output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation[order(output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation[,"mCI"],decreasing = T),,drop=F]

##################################
# Tissues

featuresType <- "output_mCI_CCLE_tissues"

listOfdrugsFiels <- dir(paste("AzureRun/",featuresType,"/",sep = ""))

output_ensemble_CTRPv2_5CV3NS_tissues <- lapply(listOfdrugsFiels, function(x){
  print(x)
  output <- tryCatch(readRDS(paste("AzureRun/",featuresType,"/",x,sep = "")),error = function(e) e)
  return(output)
})

drugsProcessed <- unlist(lapply(strsplit(listOfdrugsFiels,"_CTRPv2"),"[[",1))
drugsProcessed <- gsub("%2C",",",drugsProcessed)
drugsProcessed <- gsub("%5B","[",drugsProcessed)
drugsProcessed <- gsub("%5D","]",drugsProcessed)
drugsProcessed <- gsub("%20"," ",drugsProcessed)
drugsProcessed <- gsub("%7B","{",drugsProcessed)
drugsProcessed <- gsub("%7D","}",drugsProcessed)

names(output_ensemble_CTRPv2_5CV3NS_tissues) <- drugsProcessed

output_ensemble_CTRPv2_5CV3NS_tissues_final <- output_ensemble_CTRPv2_5CV3NS_tissues[intersect(names(keptDrugs)[keptDrugs],names(output_ensemble_CTRPv2_5CV3NS_tissues))]
mciTh = 0.6
output_ensemble_CTRPv2_5CV3NS_tissues_evaluation <- lapply(output_ensemble_CTRPv2_5CV3NS_tissues_final,function(x){
  return(x[[1]])
})

output_ensemble_CTRPv2_5CV3NS_tissues_CIs <- unlist(lapply(1:length(output_ensemble_CTRPv2_5CV3NS_tissues_evaluation), function(x){output_ensemble_CTRPv2_5CV3NS_tissues_evaluation[[x]][2]}))
output_ensemble_CTRPv2_5CV3NS_tissues_CIs <- as.numeric(output_ensemble_CTRPv2_5CV3NS_tissues_CIs)
names(output_ensemble_CTRPv2_5CV3NS_tissues_CIs) <- names(output_ensemble_CTRPv2_5CV3NS_tissues_evaluation)
output_ensemble_CTRPv2_5CV3NS_tissues_CIs <- output_ensemble_CTRPv2_5CV3NS_tissues_CIs[!is.na(output_ensemble_CTRPv2_5CV3NS_tissues_CIs)]

output_ensemble_CTRPv2_5CV3NS_tissues_final <- output_ensemble_CTRPv2_5CV3NS_tissues_final[names(output_ensemble_CTRPv2_5CV3NS_tissues_CIs)]

output_ensemble_CTRPv2_5CV3NS_tissues_evaluation <- output_ensemble_CTRPv2_5CV3NS_tissues_evaluation[names(output_ensemble_CTRPv2_5CV3NS_tissues_CIs)]
output_ensemble_CTRPv2_5CV3NS_tissues_evaluation <- do.call(rbind,output_ensemble_CTRPv2_5CV3NS_tissues_evaluation)
rownames(output_ensemble_CTRPv2_5CV3NS_tissues_evaluation) <- names(output_ensemble_CTRPv2_5CV3NS_tissues_CIs)
class(output_ensemble_CTRPv2_5CV3NS_tissues_evaluation) <- "numeric"

output_ensemble_CTRPv2_5CV3NS_tissues_evaluation <- output_ensemble_CTRPv2_5CV3NS_tissues_evaluation[order(output_ensemble_CTRPv2_5CV3NS_tissues_evaluation[,"mCI"],decreasing = T),,drop=F]

###############################
output_ensemble_CTRPv2_5CV3NS_mutation_evaluation <- cbind(output_ensemble_CTRPv2_5CV3NS_mutation_evaluation
                                                           ,"fdr"=p.adjust(output_ensemble_CTRPv2_5CV3NS_mutation_evaluation[,"Pvalue"],method = "fdr"))
output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation <- cbind(output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation
                                                          ,"fdr"=p.adjust(output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation[,"Pvalue"],method = "fdr"))
output_ensemble_CTRPv2_5CV3NS_tissues_evaluation <- cbind(output_ensemble_CTRPv2_5CV3NS_tissues_evaluation
                                                          ,"fdr"=p.adjust(output_ensemble_CTRPv2_5CV3NS_tissues_evaluation[,"Pvalue"],method = "fdr"))


output_ensemble_CTRPv2_5CV3NS_evaluation[,"fdr"] <- p.adjust(output_ensemble_CTRPv2_5CV3NS_evaluation[,"Pvalue"],method = "fdr")
output_ensemble_CTRPv2_5CV3NS_mutation_evaluation[,"fdr"] <- p.adjust(output_ensemble_CTRPv2_5CV3NS_mutation_evaluation[,"Pvalue"],method = "fdr")
output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation[,"fdr"] <- p.adjust(output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation[,"Pvalue"],method = "fdr")
output_ensemble_CTRPv2_5CV3NS_tissues_evaluation[,"fdr"] <- p.adjust(output_ensemble_CTRPv2_5CV3NS_tissues_evaluation[,"Pvalue"],method = "fdr")


commonDrugs <- Reduce(intersect,list(rownames(output_ensemble_CTRPv2_5CV3NS_evaluation)
                                     ,rownames(output_ensemble_CTRPv2_5CV3NS_mutation_evaluation)
                                     ,rownames(output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation)
                                     ,rownames(output_ensemble_CTRPv2_5CV3NS_tissues_evaluation)
))

df_comparison <- data.frame("RNAseq"=output_ensemble_CTRPv2_5CV3NS_evaluation[commonDrugs,"mCI"]
                            ,"Mutation"=output_ensemble_CTRPv2_5CV3NS_mutation_evaluation[commonDrugs,"mCI"]
                            ,"CNV"=output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation[commonDrugs,"mCI"]
                            ,"Tissues"=output_ensemble_CTRPv2_5CV3NS_tissues_evaluation[commonDrugs,"mCI"]
                            ,stringsAsFactors = F)

df_comparison_fdr <- data.frame("RNAseq"=output_ensemble_CTRPv2_5CV3NS_evaluation[commonDrugs,"fdr"]
                                ,"Mutation"=output_ensemble_CTRPv2_5CV3NS_mutation_evaluation[commonDrugs,"fdr"]
                                ,"CNV"=output_ensemble_CTRPv2_5CV3NS_cnv_all_evaluation[commonDrugs,"fdr"]
                                ,"Tissues"=output_ensemble_CTRPv2_5CV3NS_tissues_evaluation[commonDrugs,"fdr"]
                                ,stringsAsFactors = F)

df_comparison_tmp <- df_comparison[apply(df_comparison, 1, function(x){
  any(x>0.6)
}),]


df_comparison_best <- apply(df_comparison_tmp, 1, function(x){
  return(ifelse(x==max(abs(x),na.rm = T),1,0))
})

bestModels <- apply(df_comparison_best, 2, function(x){
  rownames(df_comparison_best)[which(x==1)]
})

df_comparison$Model <- "None"
df_comparison[names(bestModels),"Model"] <- bestModels

df_comparison$Model <- factor(df_comparison$Model,levels = c( "RNAseq",      "Tissues",   "Mutation", "CNV","None"))

df_comparison_best_final <- data.frame("Model"=rownames(df_comparison_best),"Freq"=apply(df_comparison_best, 1, sum),stringsAsFactors = F)

df_comparison_best_final$Model <- factor(df_comparison_best_final$Model,levels = c( "RNAseq",      "Tissues",   "Mutation", "CNV"))

df_comparison$drug <- rownames(df_comparison)

#####################
# Fig 4 - a
#####################

pdf("plots/Fig_4_a.pdf",width = 9,height = 8)
ggplot(df_comparison_best_final,aes(x="",y = Freq, fill = Model)) +
  geom_col(width = 1,position = 'stack') +
  geom_text(aes(label = paste(round(Freq / sum(Freq) * 100, 1), "%",sep = ""), x = 1),size=6,
            position = position_stack(vjust = 0.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),text = element_text(size = 20,face = "bold")) +
  labs(fill = "Model",
       x = NULL,
       y = NULL,
       title = "Distribution of best models across data types") + 
  coord_polar("y")
dev.off()


# tissue
#####################
# Fig 4 - b
#####################

pdf("plots/Fig_4_b.pdf",width = 8,height = 7)
ggplot(df_comparison,mapping = aes(x=Tissues,y = RNAseq,col=Model)) +
  geom_point()  +
  geom_text(data=subset(df_comparison, Tissues > 0.7),aes(Tissues,RNAseq,label=drug),  vjust = "inward", hjust = "inward",nudge_x = -0.005) +
  geom_abline(intercept =0 , slope = 1,col="lightgray",lty=2) +
  ylim(c(0.35,0.8)) + xlim(c(0.35,0.95)) +
  labs(fill = "Model",
       x = NULL,
       y = NULL,
       title = "Comparing tissues vs RNAseq models [CTRPv2]") + xlab("tissues basde models [rCI]") + ylab("RNAseq basde models [rCI]") +
  scale_color_manual( breaks = c("RNAseq",      "Tissues",   "Mutation", "CNV","None"),
                      values = c(scales::hue_pal()(4),"darkgray") ) +
  theme_bw() + theme(plot.title = element_text(hjust=0.5),text = element_text(size = 10,face = "bold")
                     ,axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     #                  panel.border = element_blank(),
                     panel.background = element_blank())
dev.off()



# mutation
#####################
# Fig 4 - c
#####################
pdf("plots/Fig_4_c.pdf",width = 8,height = 7)
ggplot(df_comparison,mapping = aes(x=Mutation,y = RNAseq,col=Model)) +
  geom_point()  +
  geom_text(data=subset(df_comparison, Mutation > 0.7),aes(Mutation,RNAseq,label=drug),  vjust = "inward", hjust = "inward",nudge_x = -0.005) +
  geom_abline(intercept =0 , slope = 1,col="lightgray",lty=2) +
  ylim(c(0.35,0.8)) + xlim(c(0.35,0.95)) +
  labs(fill = "Model",
       x = NULL,
       y = NULL,
       title = "Comparing Mutations vs RNAseq models [CTRPv2]") + xlab("Mutations basde models [rCI]") + ylab("RNAseq basde models [rCI]") +
  scale_color_manual( breaks = c("RNAseq",      "Tissues",   "Mutation", "CNV","None"),
                      values = c(scales::hue_pal()(4),"darkgray") ) +
  theme_bw() + theme(plot.title = element_text(hjust=0.5),text = element_text(size = 10,face = "bold")
                     ,axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     #                  panel.border = element_blank(),
                     panel.background = element_blank())
dev.off()



# CNV
#####################
# Fig 4 - d
#####################
pdf("plots/Fig_4_d.pdf",width = 8,height = 7)
ggplot(df_comparison,mapping = aes(x=CNV,y = RNAseq,col=Model)) +
  geom_point()  +
  geom_text(data=subset(df_comparison, CNV > 0.75),aes(CNV,RNAseq,label=drug),  vjust = "inward", hjust = "inward",nudge_x = -0.005) +
  geom_abline(intercept =0 , slope = 1,col="lightgray",lty=2) +
  ylim(c(0.35,0.8)) + xlim(c(0.35,0.95)) +
  labs(fill = "Model",
       x = NULL,
       y = NULL,
       title = "Comparing CNVs vs RNAseq models [CTRPv2]") + xlab("CNVs basde models [rCI]") + ylab("RNAseq basde models [rCI]") +
  scale_color_manual( breaks = c("RNAseq",      "Tissues",   "Mutation", "CNV","None"),
                      values = c(scales::hue_pal()(4),"darkgray") ) +
  theme_bw() + theme(plot.title = element_text(hjust=0.5),text = element_text(size = 10,face = "bold")
                     ,axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     #                  panel.border = element_blank(),
                     panel.background = element_blank())
dev.off()



library(reshape)
df_comparison_melt <- melt(df_comparison,id.vars = c("Model","drug"))
library(ggpubr)
my_comparisons = list(  c("RNAseq", "Tissues"),c("RNAseq", "Mutation"), c("RNAseq", "CNV") )

df_comparison_melt$variable <- factor(df_comparison_melt$variable
                                      ,levels = c("RNAseq","Tissues", "Mutation", "CNV"))

pdf("plots/Fig_4_b1.pdf",width = 5.5,height = 5)
ggplot(df_comparison_melt,aes(variable,value)) +
  geom_boxplot(aes(fill=variable), show.legend = FALSE) +
  xlab("Data type") + ylab("CI") +
  stat_compare_means(comparisons = my_comparisons,paired = T,method = "wilcox.test"
                     ,label.y = c(0.89,0.94,0.99))+
  stat_compare_means(label.y = 1.05) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),strip.text.x =element_text(size=14,face="bold"))
dev.off()


common <- intersect(rownames(output_ensemble_CTRPv2_5CV3NS_evaluation),rownames(output_ensemble_CTRPv2_5CV3NS_tissues_evaluation))
ibx <- which(output_ensemble_CTRPv2_5CV3NS_evaluation[common,"mCI"] > 0.6 | output_ensemble_CTRPv2_5CV3NS_tissues_evaluation[common,"mCI"] > 0.6)

tissuesVsRNAseq <- lapply(names(ibx), function(k){
  
  samples <- rownames(output_ensemble_CTRPv2_5CV3NS[[k]][[4]])
  unlist(mcc(x = as.factor(output_ensemble_CTRPv2_5CV3NS[[k]][[4]][samples,"Vote"]),y = as.factor(output_ensemble_CTRPv2_5CV3NS_tissues[[k]][[4]][samples,"Vote"])))
  
})

tissuesVsRNAseq_final <- do.call(rbind,tissuesVsRNAseq)

dens=density(tissuesVsRNAseq_final[,"estimate"])
pdf("plots/Fig_4_b2.pdf")
hist(tissuesVsRNAseq_final[,"estimate"],breaks = 20,main = "RNAseq Predictions vs Tissue Predictions",xlab = "MCC",xlim = c(-1,1),col = "skyblue",border = NA)
dev.off()

median(tissuesVsRNAseq_final[,"estimate"])
IQR(tissuesVsRNAseq_final[,"estimate"])
##################################
##################################
##################################
# Comparison with tissue specific models (Pan-cancer vs Lung)

ibx <- which(ccleTissues$tissueid=="lung")
ccleLung <- rownames(ccleTissues[ibx,,drop=F])

ccleLung_binary <- ccleRnaSeq_final_binary[ccleLung,]

drugs_common <- intersect(rownames(AAC_CTRPv2),rownames(output_ensemble_CTRPv2_5CV3NS_evaluation))
ibx <- which(output_ensemble_CTRPv2_5CV3NS_evaluation[drugs_common,"mCI"]>mciTh)
final_drugs_common <- rownames(output_ensemble_CTRPv2_5CV3NS_evaluation[drugs_common,])[ibx]


results_final_common_ccleLung <- lapply(final_drugs_common, function(x,output,experssion,AUCs_Drugs){
  
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
},output=output_ensemble_CTRPv2_5CV3NS_final,experssion=ccleLung_binary[,finalSetOfGenes],AUCs_Drugs=AAC_CTRPv2)




names(results_final_common_ccleLung) <- final_drugs_common

results_final_common_evaluation_ccleLung <- lapply(results_final_common_ccleLung, function(x){
  
  output <- x
  
  expressedCellLines <- rownames(output)[output[,"Vote"]==1]
  notExpressedCellLines <- rownames(output)[output[,"Vote"]==0]
  
  CI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                 ,delta.pred = 0,delta.obs = 0,outx = F)
  
  mCI <- paired.concordance.index(output[c(notExpressedCellLines,expressedCellLines),1],output[c(notExpressedCellLines,expressedCellLines),2]
                                  ,delta.pred = 0,delta.obs = 0.2,outx = F)
  
  results <- c("CI"=as.numeric(CI["cindex"]),"CI.pval"=as.numeric(CI["p.value"]),"mCI"=as.numeric(mCI["cindex"]),"mCI.pval"=as.numeric(mCI["p.value"]),"N1"=length(notExpressedCellLines),"N2"=length(expressedCellLines))
  return(results)
  
  
})


results_final_common_evaluation_ccleLung <- do.call(rbind,results_final_common_evaluation_ccleLung)

results_final_common_evaluation_ccleLung <- results_final_common_evaluation_ccleLung[order(results_final_common_evaluation_ccleLung[,"mCI"],decreasing = T),]


keptDrugs_Lung <- apply(AAC_CTRPv2, 1, function(x){
  x <- x[!is.na(x)]
  x <- x[intersect(names(x),ccleLung)]
  ((sum(x>0.2,na.rm = T)/sum(x>=0,na.rm = T)) >0.1)
})

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
mciTh = 0.6

output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation <- lapply(output_ensemble_CTRPv2_5CV3NS_LUNG_final,function(x){
  return(x[[1]])
})

output_ensemble_CTRPv2_5CV3NS_LUNG_CIs <- unlist(lapply(1:length(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation), function(x){output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[[x]][2]}))

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


common <- intersect(rownames(output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation),rownames(results_final_common_evaluation_ccleLung))

pdf("plots/PanCancerVsLung.pdf")
plot(x = output_ensemble_CTRPv2_5CV3NS_LUNG_evaluation[common,"mCI"],y = results_final_common_evaluation_ccleLung[common,"mCI"],ylim = c(0.35,0.95),xlim = c(0.35,0.95)
     ,main = "Comparing Lung-specific vs pan-cancer models [CTRPv2]",xlab = "lung-specific models [CI]",ylab = "pan-cancer models [CI]"
     ,pch=16,col=rgb(red = 0.001, green = 0.4, blue = 0.4, alpha = 0.8))
lines(x = c(0,1), y = c(0,1))
abline(h = 0.6,v = 0.6,lty=2,col="red")
rect(xleft = 0,0,0.6,0.6,col = "gray",border = NA,density = 30)
dev.off()


