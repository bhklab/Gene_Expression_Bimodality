

data <- readRDS("UHN_PDX_BR_log2TPM_0.001.rda")

data2 <- readRDS("PDX_expression_all.rds")


cbind(pData(data)["final_processedID"],rownames(pData(data)),colnames(exprs(data)))

table(pData(data)$SAMPLE.TYPE)

dim(data)


data3 <- readRDS("Cescon_PDXs.rda")

data3 <- data3$rnaseq

pData(data3)

colnames(data3)

saveRDS(data3,file = "Cescon_PDXs_ALL_TPM.rda")


ibx <- match(colnames(data3),pData(data)$rnaseqID)

notInData <- colnames(data3)[is.na(ibx)]

notInData <- notInData[! notInData%in% c("BPTO_51_S4","BPTO_63_S3","BXTO_64_S2","BXTO_81_S3","Sample_21046","Sample_23223","Sample_48342","Sample_Project_Cescon_RNA_39","Sample_REF006" )]


tmp_data <- data3[,notInData]
colnames(tmp_data)[7] <- "69693"


tmp_data <- tmp_data[,c("104987_P1_S5","73263-CX_Res_p8_B_S8","REF-S-032B_S1","REF-S-035_S2","69693")]


p_data_tmp <- data.frame("final_processedID"=colnames(tmp_data),row.names = colnames(tmp_data),stringsAsFactors = F)


sum(duplicated(c(sampleNames(tmp_data),sampleNames(data))))
