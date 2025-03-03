#run DESeq2
library(DESeq2)
library(tximport)
library(jsonlite)
library(readr)
library(apeglm)
library(clusterProfiler)
library(biomaRt)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(EnhancedVolcano)
library(genefilter)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(pheatmap)
library(tibble)
#prepare salmon mapping data for deseq2 (per http://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#transcript-abundance-files-and-tximport-tximeta)
files <- c("naive1.sf", "naive2.sf", "naive3.sf", "naive4.sf", "naive5.sf", "infected6.sf", "infected7.sf", "infected8.sf", "infected9.sf", "infected10.sf")
samples <- read.table("2024.1.19_deseq_samples.txt", header=TRUE)
samples$condition <- factor(rep(c("naive", "infected"),each=5))
rownames(samples) <- samples$run
samples[,c("run","condition")]
tx2gene <- read.table("salmon_tx2gene.tsv", header = FALSE)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
#create deseq dataset
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
#remove low read count genes so that at least 10 reads are found in at least 2 samples
keep <- rowSums(counts(ddsTxi) >= 10) >= 2
ddsTxi_keep <- ddsTxi[keep,]
#refactor conditions to indicate the control to be compared to
ddsTxi_keep$condition <- relevel(ddsTxi_keep$condition, ref = "naive")
#Differential expression analysis
ddsTxi_keep_deseq <- DESeq(ddsTxi_keep)

##exploratory data analysis
plotDispEsts(ddsTxi_keep_deseq)
#pcaplot of rlog transformed dseq2 results
rld <- rlog(ddsTxi_keep_deseq)
plotPCA(rld, intgroup = c("run", "condition"))

#sampledistance matrix of rlog Euclidean distances of samples
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition, rld$run, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(sampleDists)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )

#heatmap and hierarchial clustering of the top 35 most variable genes
topVarGenes <- head(order(-rowVars(assay(rld))),35)
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$condition ]
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(rld$condition,"-",rld$run)
heatmap.2(mat, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row")

#plot pvalues
hist(resInfected_0.05$pvalue, breaks=20, col="grey50", border="white")

#adjust the alpha value from default 0.01 to 0.05 and then apply apeglm lfc shrinkage
resInfected_0.05 <- results(ddsTxi_keep_deseq, alpha=0.05, name="condition_infected_vs_naive")
resInfected_0.05 <- lfcShrink(ddsTxi_keep_deseq, coef=2, type="apeglm", res=resInfected_0.05)

#summary of how many results have p value below 0.05
sum(resInfected_0.05$padj < 0.05, na.rm=TRUE)
#total genes
sum(rowSums(counts(ddsTxi_keep_deseq) > 0) > 0)
#MA-plot
plotMA(resInfected_0.05, ylim=c(-2,2))
plotMA(resInfected_0.05, ylim=c(-2,2))

#####remove potential outliersf naive2 and infected10 and adjust lfcThreshold
#prepare salmon mapping data for deseq2 (per http://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#transcript-abundance-files-and-tximport-tximeta)
files <- c("naive1.sf", "naive3.sf", "naive4.sf", "naive5.sf", "infected6.sf", "infected7.sf", "infected8.sf", "infected9.sf")
samples <- read.table("2024.1.22_deseq_samples_no_naive2_or_infected10.txt", header=TRUE)
samples$condition <- factor(rep(c("naive", "infected"),each=4))
rownames(samples) <- samples$run
samples[,c("run","condition")]
tx2gene <- read.table("salmon_tx2gene.tsv", header = FALSE)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
#create deseq dataset
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
#remove low read count genes so that at least 10 reads are found in at least 4 samples
keep <- rowSums(counts(ddsTxi) >= 10) >= 4
ddsTxi_keep <- ddsTxi[keep,]
#refactor conditions to indicate the control to be compared to
ddsTxi_keep$condition <- relevel(ddsTxi_keep$condition, ref = "naive")
#Differential expression analysis
ddsTxi_keep_deseq <- DESeq(ddsTxi_keep)

##exploratory data analysis
plotDispEsts(ddsTxi_keep_deseq)
#pcaplot of rlog transformed dseq2 results
rld <- rlog(ddsTxi_keep_deseq)
plotPCA(rld, intgroup = c("run", "condition"))

#sampledistance matrix of rlog Euclidean distances of samples
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition, rld$run, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(sampleDists)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )

#heatmap and hierarchial clustering of the top 35 most variable genes
topVarGenes <- head(order(-rowVars(assay(rld))),35)
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$condition ]
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(rld$condition,"-",rld$run)
heatmap.2(mat, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row")

#plot pvalues
hist(resInfected_0.05$pvalue, breaks=20, col="grey50", border="white")

#adjust the alpha value from default 0.1 to 0.05 and then apply apeglm lfc shrinkage
resInfected_0.05 <- results(ddsTxi_keep_deseq, alpha=0.05, name="condition_infected_vs_naive", lfcThreshold = 0.585, altHypothesis = "greaterAbs")
resInfected_0.05_lfcshrink <- lfcShrink(ddsTxi_keep_deseq, coef=2, type="apeglm", res=resInfected_0.05)

write.csv(resInfected_0.05_lfcshrink, "2024.7.1_resInfected_0.05_lfcshrink.csv")

#summary of how many results have p value below 0.05
sum(resInfected_0.05$padj < 0.05, na.rm=TRUE)
sum(resInfected_0.05_lfcshrink$padj < 0.05, na.rm=TRUE)
summary(resInfected_0.05)
summary(resInfected_0.05_lfcshrink)
#total genes
sum(rowSums(counts(ddsTxi_keep_deseq) > 0) > 0)
#MA-plot
plotMA(resInfected_0.05, ylim=c(-10,10))
plotMA(resInfected_0.05, ylim=c(-10,10))

#convert ensembl id to entrezid
resInfected_0.05_lfcshrink$ENTREZID <- mapIds(EnsDb.Mmusculus.v79, keys = rownames(resInfected_0.05_lfcshrink), column = 'ENTREZID', keytype = 'GENEID', multiVals="first")
rownames(resInfected_0.05_lfcshrink) <- ifelse(is.na(resInfected_0.05_lfcshrink$ENTREZID), rownames(resInfected_0.05_lfcshrink), resInfected_0.05_lfcshrink$ENTREZID)
anyDuplicated(names(resInfected_0.05_lfcshrink))

# subset significant and more than 1log2 up/down
resInfected_0.05_lf2_up_entrez <- subset(resInfected_0.05_lfcshrink, padj < 0.05 & log2FoldChange >= 0.585)
resInfected_0.05_lf2_down_entrez <- subset(resInfected_0.05_lfcshrink, padj < 0.05 & log2FoldChange <= -0.585)

resInfected_0.05_lf2_sig_entrez <- subset(resInfected_0.05_lfcshrink, padj < 0.05)

#get list of genes for upregulated and downregulated genes
res_Infected_sig_entrez_names <- row.names(resInfected_0.05_lf2_sig_entrez)
res_Infected_all_genes_entrez_names <- row.names(resInfected_0.05_lfcshrink)

res_Infected_down_entrez_names <- row.names(resInfected_0.05_lf2_down_entrez)


#gene ontology GO analysis for upregulated genes
ego_resInfected_0.05_lf2_up <- enrichGO(gene = rownames(resInfected_0.05_lf2_up_entrez), 
                                        OrgDb         = org.Mm.eg.db,
                                        keyType ='ENTREZID',
                                        ont           = "ALL",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        qvalueCutoff  = 0.05,
                                        readable      = TRUE,
                                        universe = rownames(resInfected_0.05_lfcshrink))

#gene ontology GO analysis for downregulated genes
ego_resInfected_0.05_lf2_down <- enrichGO(gene = rownames(resInfected_0.05_lf2_down_entrez), 
                                          OrgDb         = org.Mm.eg.db,
                                          keyType ='ENTREZID',
                                          ont           = "ALL",
                                          pAdjustMethod = "BH",
                                          pvalueCutoff  = 0.05,
                                          qvalueCutoff  = 0.05,
                                          readable      = TRUE,
                                          universe = rownames(resInfected_0.05_lfcshrink))

#barplot of GO enrichment
barplot(ego_resInfected_0.05_lf2_up, showCategory = 40, minGSSize = 5, title= "Reinfected up lf2")
barplot(ego_resInfected_0.05_lf2_down, showCategory = 40, minGSSize = 5,title= "Reinfected down lf2")

#volcano plot with enhancedvolcano package (per https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html)

volcano_Infected <- EnhancedVolcano(resInfected_0.05_lfcshrink,lab = NA, x = 'log2FoldChange', y = 'pvalue', title = 'infected_v_naive', FCcutoff = 0.585, pCutoff = 0.05, pCutoffCol = 'padj')

print(volcano_Infected)

#convert ENSMBL to GENE SYMBOL and output results in order of log2foldchange
resInfected_0.05$SYMBOL <- mapIds(EnsDb.Mmusculus.v79, keys = rownames(resInfected_0.05), column = 'SYMBOL', keytype = 'GENEID', multiVals="first")
rownames(resInfected_0.05) <- ifelse(is.na(resInfected_0.05$SYMBOL), rownames(resInfected_0.05), resInfected_0.05$SYMBOL)

resInfected_0.05_reorder <- resInfected_0.05[order(resInfected_0.05$log2FoldChange),]
write.csv(as.data.frame(resInfected_0.05_reorder), file="reinfection_0.05_lfcthreshold_2_no_outliers_significant_results.csv")



write.csv(as.data.frame(ego_resInfected_0.05_lf2_down), file= "ego_resInfected_0.05_lf2_down_for_plotting.csv")
write.csv(as.data.frame(ego_resInfected_0.05_lf2_up), file= "ego_resInfected_0.05_lf2_up_for_plotting.csv")
write.csv(res_Infected_sig_entrez_names, file="reinfection_significant_genes.csv")
write.csv(res_Infected_all_genes_entrez_names, file="reinfection_all_genes.csv")


##plotting GO data
#read in csv files
ego_Infected_lfc_0.05_lf2_up_for_plotting <- read.csv("ego_resInfected_0.05_lf2_up_for_plotting.csv")
ego_Infected_lfc_0.05_lf2_down_for_plotting <- read.csv("ego_resInfected_0.05_lf2_down_for_plotting.csv")


# Order the data by p.adjust values
data_up <- ego_Infected_lfc_0.05_lf2_up_for_plotting[order(ego_Infected_lfc_0.05_lf2_up_for_plotting$p.adjust), ]
data_down <- ego_Infected_lfc_0.05_lf2_down_for_plotting[order(ego_Infected_lfc_0.05_lf2_down_for_plotting$p.adjust), ]

# Create a factor variable for category with the desired order
data_up$category <- factor(data_up$category, levels = rev(data_up$category))
data_down$category <- factor(data_down$category, levels = rev(data_down$category))

# Create the bar plot upregulated
ggplot(data_up, aes(x = count, y = category, fill = p.adjust)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_gradientn(colors = c("darkorchid4","#ff3edf"), guide = guide_colourbar(reverse = TRUE)) +
  labs(title = expression(bold("Upregulated GO Categories")),
       x = expression(bold("Gene Count")),
       y = expression(bold("Category"))) +
  theme_minimal()+
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 16),
    text = element_text(size = 14),
    axis.text= element_text(color="black")
  )

# Create the bar plot downregulated
ggplot(data_down, aes(x = count, y = category, fill = p.adjust)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_gradientn(colors = c("chartreuse4","chartreuse"), guide = guide_colourbar(reverse = TRUE)) +
  labs(title = expression(bold("Downregulated GO Categories")),
       x = expression(bold("Gene Count")),
       y = expression(bold("Category"))) +
  theme_minimal()+
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 16),
    text = element_text(size = 14),
    axis.text= element_text(color="black")
  )

##heatmap of tcell genes zscore normalized counts
#prep data
normalized_counts <- counts(ddsTxi_keep_deseq , normalized=TRUE)
tcellgenenames <- read.csv("mouse_tcell_gene_names.csv", header = T)
tcellgenes_sig_list <- read.csv("mouse_tcell_genes.csv", header = T)
tcellgenes_sig_list <-  as.character(tcellgenes_sig_list$ENSEMBL)
selected_counts <- as.data.frame(normalized_counts[tcellgenes_sig_list, ])
selected_counts <- rownames_to_column(selected_counts, var = "ENSEMBL")
merged_selected_counts <- merge(selected_counts, tcellgenenames, by= "ENSEMBL")
merged_selected_counts$ENSEMBL <- merged_selected_counts$GENE
merged_selected_counts <- merged_selected_counts[, -which(names(merged_selected_counts) == "GENE")]
merged_selected_counts_data <- column_to_rownames(merged_selected_counts, var = "ENSEMBL")


# Define a function to calculate z-scores
calculate_z_scores <- function(x) {
  (x - mean(x)) / sd(x)
}

# Apply the z-score function to each row (gene)
z_scores <- t(apply(merged_selected_counts_data, 1, calculate_z_scores))

#generate heatmap
pheatmap(z_scores,
         cluster_rows = TRUE,     # Cluster rows (genes)
         cluster_cols = FALSE,     # Cluster columns (samples)
         scale = "none",          # Data is already scaled
         annotation_col = NULL,   # Add column annotations if needed
         show_rownames = TRUE,    # Show row (gene) names
         show_colnames = TRUE     # Show column (sample) names
)

##heatmap of mucosal genes zscore normalized counts
#prep data
normalized_counts <- counts(ddsTxi_keep_deseq , normalized=TRUE)
mmggenenames <- read.csv("mouse_mucosal_gene_names.csv", header = T)
mmggenes_sig_list <- read.csv("mouse_mucosal_genes.csv", header = T)
mmggenes_sig_list <-  as.character(mmggenes_sig_list$ENSEMBL)
selected_counts <- as.data.frame(normalized_counts[mmggenes_sig_list, ])
selected_counts <- rownames_to_column(selected_counts, var = "ENSEMBL")
merged_selected_counts <- merge(selected_counts, mmggenenames, by= "ENSEMBL")
merged_selected_counts$ENSEMBL <- merged_selected_counts$GENE
merged_selected_counts <- merged_selected_counts[, -which(names(merged_selected_counts) == "GENE")]
merged_selected_counts_data <- column_to_rownames(merged_selected_counts, var = "ENSEMBL")

# Apply the z-score function to each row (gene)
z_scores <- t(apply(merged_selected_counts_data, 1, calculate_z_scores))

#generate heatmap
pheatmap(z_scores,
         cluster_rows = TRUE,
         cluster_cols = FALSE,   
         scale = "none",
         annotation_col = NULL,
         show_rownames = TRUE,
         show_colnames = TRUE
)


##heatmap of apoptosis genes zscore normalized counts
#prep data
normalized_counts <- counts(ddsTxi_keep_deseq , normalized=TRUE)
apoptosisgenenames <- read.csv("mouse_apoptosis_genes_names.csv", header = T)
apoptosisgenes_sig_list <- read.csv("mouse_apoptosis_genes.csv", header = T)
apoptosisgenes_sig_list <-  as.character(apoptosisgenes_sig_list$ENSEMBL)
selected_counts <- as.data.frame(normalized_counts[apoptosisgenes_sig_list, ])
selected_counts <- rownames_to_column(selected_counts, var = "ENSEMBL")
merged_selected_counts <- merge(selected_counts, apoptosisgenenames, by= "ENSEMBL")
merged_selected_counts$ENSEMBL <- merged_selected_counts$GENE
merged_selected_counts <- merged_selected_counts[, -which(names(merged_selected_counts) == "GENE")]
merged_selected_counts_data <- column_to_rownames(merged_selected_counts, var = "ENSEMBL")

# Apply the z-score function to each row (gene)
z_scores <- t(apply(merged_selected_counts_data, 1, calculate_z_scores))

#generate heatmap
pheatmap(z_scores,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         annotation_col = NULL,
         show_rownames = TRUE,
         show_colnames = TRUE
)

























