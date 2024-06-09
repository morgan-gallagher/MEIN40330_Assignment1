library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(stringr) 
library(DESeq2)
library(pheatmap)


#read in raw files
df_clinical_patient <- read.delim("C:/Users/morgs/OneDrive/Desktop/PM_Ass1/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", row.names=1, comment.char="#")
rnaseq <- read.delim("data_mrna_seq_v2_rsem.txt")


#remove duplicated rows from rnaseq
keep = !duplicated(rnaseq[,1])
rnaseq = rnaseq[keep,]
rnaseq <- rnaseq[-1, ]

#set rownames of rnaseq
rownames(rnaseq)  = rnaseq[,1]

#remove last 3 characters of patient ids in rnaseq
colnames_rnaseq <- colnames(rnaseq)[-c(1,2)]
colnames_remove_end <- lapply(colnames_rnaseq, function(x) {str_sub(x, 1, -4)})
colnames(rnaseq)[3:length(colnames(rnaseq))] <- colnames_remove_end

#change patient id's to be . instead of - from clinical data
rownames_clinical <- rownames(df_clinical_patient)
rownames_sub <- lapply(rownames_clinical, function(x) {gsub("-", "\\.", x)})
rownames(df_clinical_patient) <- rownames_sub


#subset dataset to female Black and White
df_clinical_patient <- df_clinical_patient |>
  filter(SEX == "Female") |>
  filter(RACE == "White" | RACE == "Black or African American")

#make sure same patients are included in clinical and rnaseq df's
rem_rows <- setdiff(rownames(df_clinical_patient), colnames(rnaseq)[- c(1, 2)])
rem_rows1 <- setdiff(colnames(rnaseq)[- c(1, 2)], rownames(df_clinical_patient))
rem_rows_all <- cbind(rem_rows, rem_rows1)
df_clinical_patient <- df_clinical_patient[!rownames(df_clinical_patient) %in% rem_rows_all, ]
rnaseq <- rnaseq[, !colnames(rnaseq) %in% rem_rows_all]


#get metadata on race
info <- data.frame(colnames(rnaseq)[- c(1, 2)], df_clinical_patient$RACE)

#select only black and white patients
colnames(info)[1] <- "PATIENT_ID"
colnames(info)[2] <- "RACE"

#round rnaseq data
rnaseq_count <- round(rnaseq[, 3:length(colnames(rnaseq))])

#replace negative and NA values with 0 for count data
rnaseq_count[is.na(rnaseq_count)] = 0  
rnaseq_count[rnaseq_count<0] = 0

#make metadata factor
info$RACE <- as.factor(info$RACE)


#write all files to csv
write.csv(info, 'race_meta.csv')
info <- read.csv('race_meta.csv', row.names = 1)

#DIFFERENTIAL GENE EXPRESSION

#make DESeq object from dataframes
dds <- DESeqDataSetFromMatrix(countData = rnaseq_count,
                              colData = info,
                              design = ~ RACE)

#set threshold for group (number of people with more than 10 counts)
smallestGroupSize <- 3

#only keep genes where 3 or more people have true values
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize

dds <- dds[keep, ]

#Run DESeq
ddsDE <- DESeq(dds)

#DESeq results
res <- results(ddsDE, alpha = 0.05)
summary(res)

#order DESeq results by padj value ascending
resOrdered <- res[order(res$padj),]

#export ordered DESeq results
write.csv(resOrdered, 'deSeq_order_output.csv')
ordered_output <- read.csv('deSeq_order_output.csv', row.names = 1)

#plotMA(ddsDE, ylim = c(-5, 5))

#label genes as significant or not
ordered_output$sig <- ifelse(ordered_output$padj <= 0.05, "yes", "no")

#only keep genes that are sig
signif <- subset(ordered_output, padj <= 0.05)

#variance stabilizing transformation 
vsd <- vst(ddsDE, blind = FALSE)
transformed_counts <- assay(vsd)

#pca
pca <- prcomp(t(transformed_counts))
df_pca <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
percentVar <- round(100*pca$sdev^2 / sum(pca$sdev^2), 1)
pca_patient_ids <- rownames(df_pca)
df_pca$patient_ids <- pca_patient_ids
df_clinical_patient_ids <- rownames(df_clinical_patient)
df_clinical_patient$patient_ids <- df_clinical_patient_ids

#combine pc coordintes with clinical data for plotting
df_pca_complete <- left_join(df_pca, df_clinical_patient, by = 'patient_ids')

#replace blank entries for subtype with NA
df_pca_complete <- df_pca_complete |>
  mutate(SUBTYPE = na_if(SUBTYPE, ""))

#plot pca colored by race
ggplot(df_pca_complete, aes(x = PC1, y = PC2, color = RACE)) +
  geom_point(size = 3) +
  xlab(paste0('PC1 (', percentVar[1], '% Variance)')) +
  ylab(paste0('PC2 (', percentVar[2], '% Variance)')) +
  ggtitle("PCA of Gene Expression - Race") +
  theme_classic()

#plot pca colored by subtype
ggplot(df_pca_complete, aes(x = PC1, y = PC2, color = SUBTYPE)) +
  geom_point(size = 3) +
  xlab(paste0('PC1 (', percentVar[1], '% Variance)')) +
  ylab(paste0('PC2 (', percentVar[2], '% Variance)')) +
  theme_classic() +
  ggtitle("PCA of Gene Expression - Subtype of Tumor")

#get top 10 upregulated and downregulated genes
top10_downreg <- signif[order(signif$log2FoldChange), ]
top10_downreg <- top10_downreg[1:10, ]
top10_upreg <- signif[order(signif$log2FoldChange, decreasing = TRUE), ]
top10_upreg <- top10_upreg[1:10, ]


#get names of up and downregulated genes
up_genes <- rownames(top10_upreg)
down_genes <- rownames(top10_downreg)

#get rnaseq data for these genes
all_genes <- cbind(up_genes, down_genes)
rnaseq_up_down <- rnaseq_count[all_genes, ]

#log transform genes
rnaseq_up_down_transform <- log2(rnaseq_up_down + 1)

#extract race info for annotation
race_info <- data.frame(info$RACE)
colnames(race_info) <- "Race"
row.names(race_info) <- info$PATIENT_ID

#plot heatmap
pheatmap(rnaseq_up_down_transform, scale = 'column', show_colnames = FALSE, annotation_col = race_info, main = "Gene Expression of Top 20 DEG")

