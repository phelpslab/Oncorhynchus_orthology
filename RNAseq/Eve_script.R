BiocManager::install("DESeq2")
install.packages("remotes")
remotes::install_github("Jmendo12/evemodel")
library(DESeq2)
library(evemodel)

Omykiss_RNA_loc <- read.csv("Omkyiss_RNA_loc.csv")
Omykiss_protein <- read.csv("Omykiss_protein_annotations.csv")
merge <- inner_join(Omykiss_RNA_loc,Omykiss_protein, by="Locus")
write.csv(merge, file = "Omykiss_annotations_txi.csv")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
library(tximport)

Trout_samples <- read.table("Trout_samples.txt", header = TRUE)
Trout_samples
Trout_files <- file.path(Trout_samples$Samples, "quant.sf")
names(Trout_files) <- paste0("sample", 1:42)
all(file.exists(Trout_files))
Trout_tx2gene <- read.csv("Omykiss_txt2gene.csv")
head(Trout_tx2gene)
Trout_txi <- tximport(Trout_files, type = 'salmon',tx2gene=Trout_tx2gene)
names(Trout_txi)
tail(Mcav_txi$counts)
Mcav_txi_sum <- summarizeToGene(Mcav_txi,Mcav_tx2gene)
all.equal(Mcav_txi$counts, Mcav_txi_sum$counts)

write.csv(Trout_txi, file = "Redband_counts.csv")

#Read files in to DESEQ2
countData <- as.matrix(read.csv("Redband_counts.csv", row.names="Protein"))
colnames(countData)
colData <- read.csv("Design.csv",header = TRUE)
colData
#Normalize by species ID
#Because Salmon gives abundances in terms of decimals so you need to round each transcript abundance
dds <- DESeqDataSetFromMatrix(countData = round(countData), colData = colData, design = ~Strain)
#Remove low abundnace transctipts
dds <- dds[ rowMeans(counts(dds)) > 10, ]
#Noramlize
vst <- vst(dds)
#Extract normalized counts
vvst <- assay(vst)
write.csv(vvst, file ="Normalized_counts.csv" )
#Read counts in
#I also provided the normalized coutns in case you didn't want to run everything
exprMat <- as.matrix(read.csv("Normalized_counts.csv", row.names = "Protein"))
#Same species tree we already used in Cafe
speciesTree <- read.tree("Redband_tree.tre")
#Before running EVE make sure that your samples have the species name with and underscore then the sample number and that the species name matches what is in the tree
plot(speciesTree)
speciesTree$tip.label
colnames(exprMat)
colSpecies <- sub("_.*$","",colnames(exprMat))
colSpecies
res <- betaSharedTest(tree = speciesTree, gene.data = exprMat, colSpecies = colSpecies)
#Eve takes a really long time to run ~1.5 hours on this size of data so I provide the output
results <- res$indivBetaRes$par
pval = pchisq(res$LRT,df = 1,lower.tail = F)
pvalue <- pval
Eve_results <- cbind(rownames(exprMat),results,pvalue)
write.csv(Eve_results, file = "Expression_shifts.csv")


Eve_redband <- read.csv("Expression_shifts.csv")
ggdensity(Eve_redband$beta)
log_beta <- log10(Eve_redband$beta)
Normalized_Eve_data <- cbind(Eve_redband,log_beta)
ggdensity(Normalized_Eve_data$beta)
Divergence <- filter(Eve_redband, beta < 57 & pvalue < 0.05)
#783
Redip <- read.csv("Omykiss_redip_gene_ids.csv")
Redip_merge <- inner_join(Eve_redband,Redip,by="Protein")
head(Redip_merge)
Redip_sig <- filter(Redip_merge, beta < 57 & pvalue < 0.05)
#169 put of 5734
mean(Redip_sig$beta)
#1.83
Ohnos <- read.csv("Mykiss_universal_ohnos.csv")
Ohno_merge <- inner_join(Eve_redband,Ohnos,by="Protein")
Ohno_sig <- filter(Ohno_merge, pvalue < 0.05)
#158 out of 6013
mean(Ohno_sig$beta)
#1.927905

Redipping <- read.csv("Omykiss_redipping_IDs.csv")


t.test(Ohno_sig$beta, subset$beta)
t.test(Redip_sig$beta, subset$beta)

subset <- Eve_redband[sample(nrow(Eve_redband), 6013), ]
random_sig <- filter(subset,  pvalue < 0.05)
#177



Ohno_divergence <- Eve_redband$Protein%in%Ohnos$Protein
Redip_divergence <- Eve_redband$Protein%in%Redip$Protein
Redipping_divergence <- Eve_redband$Protein%in%Redipping$Omykiss.y
Class_divergence <- cbind(Eve_redband,Ohno_divergence,Redip_divergence,Redipping_divergence)
write.csv(Class_divergence, file = "Class_divergence.csv")

Divergence_test <- read.csv("Class_divergence.csv")
Filter_Divergence_test <- filter(Divergence_test, Expression_divergence > 0.336)
Ohno_nummber <- filter(Filter_Divergence_test, Ohno_divergence == TRUE)
#2954
Sig_ohnos <- filter(Ohno_nummber,pvalue < 0.05)
#158
Redip_number <- filter(Filter_Divergence_test, Redip_divergence == TRUE)
#5734
Sig_redip <- filter(Redip_number, pvalue < 0.05)
#169
Other_number <- filter(Filter_Divergence_test, Redip_divergence == FALSE & Ohno_divergence == TRUE & Redipping_divergence == FALSE)
#2954
ggdensity(Filter_Divergence_test$Expression_divergence)
anova <- aov(Expression_divergence~Class, data = Filter_Divergence_test)
summary(anova)
TukeyHSD(anova)
subset <- Filter_Divergence_test[sample(nrow(Filter_Divergence_test), 2954), ]
random_sig <- filter(subset,  pvalue < 0.05)
#185



Expression_shifts_00019 <- matrix(c(158,2954,169,5734),
                                  nrow = 2,
                                  dimnames = list(c("Hits","Total"),
                                                  c("Hits","Total")))
fisher.test(Expression_shifts_00019, alternative = 'two.sided')
# odds ratio 1.814703 
# p-value = 2.108e-07

Filtered_EVE <- filter(Eve_redband, beta < 99)

Duplicates <- read.csv("Omykiss_ohno_gene_ids.csv")
Sig_duplicates <- inner_join(Filtered_EVE,Duplicates,by="Protein")
Sig_duplicate_number <- filter(Sig_duplicates, pvalue < 0.05)
test <- filter(Filtered_EVE, pvalue < 0.05)

Expression_shifts_00019 <- matrix(c(431,7436,169,5734),
                                  nrow = 2,
                                  dimnames = list(c("Hits","Total"),
                                                  c("Hits","Total")))
fisher.test(Expression_shifts_00019, alternative = 'two.sided')
# odds ratio 1.96645 
# p-value = 3.948e-14




Quants <- read.csv("Normalized_counts.csv")
Sig_ohno_counts <- inner_join(Quants,Sig_ohnos,by="Protein")
write.csv(Sig_ohno_counts,file = "Sig_ohno_counts.csv")
Ohnos_heatmap <- read.csv("Sig_ohno_counts.csv",row.names = 1)
Spruce <- colorRampPalette(c("#C9D4D4", "#A3BDB8", "#8DA992", "#394036", "#1E1D16"))
Spruce1 <- colorRampPalette(c("#1E1D16", "#394036", "#8DA992", "#A3BDB8", "#C9D4D4"))
Lichen <- colorRampPalette(c("#7BA2C7", "#606F64", "#0D1312", "#D35D28", "#FFB160"))

pheatmap(Ohnos_heatmap,
        show_rownames = FALSE,
        color   = Lichen(500),
        scale = "row"
        )

