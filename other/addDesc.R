library( "biomaRt" )
args <- commandArgs(TRUE)
res <- read.table(file=args[1], sep="\t",row.names=1,header=TRUE)
rn <- rownames(res)
rownames(res) <- gsub("\\.\\d+","",rownames(res))

ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description"), filters = "ensembl_gene_id",values = rownames(res),mart = ensembl )

# ensembl = useMart( "ensembl", dataset = "mmusculus_gene_ensembl")
# genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "description"), filters = "ensembl_gene_id",values = rownames(res),mart = ensembl )
idx <- match( rownames(res), genemap$ensembl_gene_id )
res$entrez <- genemap$entrezgene[ idx ]
# res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
# res$mgi_symbol <- genemap$mgi_symbol[ idx ]
res$description <- genemap$description[ idx ]
res$ENSEMBL_ID <- rn;
# res <- res[,c("ENSEMBL_ID","Gene.Name","Gene.Type","description","Chr","Strand","Start","End","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","entrez","hgnc_symbol")]
res <- res[,c("ENSEMBL_ID","Gene.Name","Gene.Type","description","Chr","Strand","Start","End","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","entrez")]
write.table(res, file=args[1], quote=FALSE, sep="\t", row.names=FALSE)

