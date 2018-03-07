library('DESeq2')
library( "gplots" )
library("ggplot2")
library("RColorBrewer")
library("BiocParallel")
library(pheatmap)
library("IHW")
library( "biomaRt" )
register(MulticoreParam(32))
lnames <- load("GSF1635.Counts_DESeq2.RData")

res <- as.data.frame(res)
rn <- rownames(res)
rownames(res) <- gsub("\\.\\d+","",rownames(res))
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description"), filters = "ensembl_gene_id",values = rownames(res),mart = ensembl )
idx <- match( rownames(res), genemap$ensembl_gene_id )
res$entrez <- genemap$entrezgene[ idx ]
res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
res$description <- genemap$description[ idx ]
res$ENSEMBL_ID <- rn;

rlogMat2 <- as.data.frame(rlogMat)
rn <- rownames(rlogMat2)
rownames(rlogMat2) <- gsub("\\.\\d+","",rownames(rlogMat2))
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description"), filters = "ensembl_gene_id",values = rownames(rlogMat2),mart = ensembl )
idx <- match( rownames(rlogMat2), genemap$ensembl_gene_id )
rlogMat2$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
rlogMat2 <- rlogMat2[! is.na(rlogMat2$hgnc_symbol),]


# "DAB2IP Expressing" <- minus
# "Control" <- plus
pdf("heatmaps.pdf")
mat <- rlogMat2[rlogMat2$hgnc_symbol=="PROM1" | rlogMat2$hgnc_symbol=="TWIST1" | rlogMat2$hgnc_symbol=="ALDH1A1" | rlogMat2$hgnc_symbol=="LGR5" | rlogMat2$hgnc_symbol=="VIM" | rlogMat2$hgnc_symbol=="CD24",]
rownames(mat) <- mat$hgnc_symbol
mat <- mat[,-7]
mat <- mat - rowMeans(mat)

df <- as.data.frame(colData(rld)[,c("condition")])
colnames(df) <- "CellType"
rownames(df) <- colnames(mat)
df3 <- df
df3$CellType <- c("DAB2IP Expressing", "DAB2IP Expressing", "DAB2IP Expressing", "Control","Control","Control")
rownames(df3) <- c("DAB2IP Expressing 1", "DAB2IP Expressing 2", "DAB2IP Expressing 3", "Control 1","Control 2","Control 3")
df <- df3
colnames(mat) <- rownames(df)
mycols <- brewer.pal(8, "Dark2")[1:length(unique(df$CellType))]
ann_colors <- list(CellType=c(mycols))
names(ann_colors$CellType) <- unique(df$CellType)

title <- "Heatmap"
ph <- pheatmap(mat, color = colorRampPalette(c("navy", "white", "red"))(256), annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=TRUE, cluster_cols=FALSE, scale="row")
#ph <- pheatmap(mat, annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=4, fontsize_col=6, show_rownames=TRUE, cluster_cols=FALSE, scale="row")


mat <- rlogMat2[rlogMat2$hgnc_symbol=="ABCA2" | rlogMat2$hgnc_symbol=="ABCA3" | rlogMat2$hgnc_symbol=="ABCA5" | rlogMat2$hgnc_symbol=="ABCB11" | rlogMat2$hgnc_symbol=="ABCG1",]

savehistory("heatmaps.R")
