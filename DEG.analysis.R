# Analysis for differentially expressed genes

library(DESeq2)
library(tximport)
library(doParallel)
library(BiocParallel)

te.oper <- "./12.rsem.re/GNC/"
exp.files <- list.files(te.oper1,full.names=T)
file.nms <- do.call(rbind,strsplit(exp.files,"/"))
file.nms <- file.nms[,ncol(file.nms)]

exp.files <- paste(exp.files1,"/",file.nms,".result.genes.results",sep="")

exp.files <- c(exp.files)
te.cns <- gsub(".result.genes.results","",gsub(".*/","",exp.files))
tx.rsem <- tximport(exp.files, type = "rsem")
colnames(tx.rsem$"abundance") <- colnames(tx.rsem$"counts") <- colnames(tx.rsem$"length") <- te.cns

sample.in <- gsub(" ","",as.matrix(read.table("./sample.info",sep='\t')))
sample.in <- sample.in[sample.in[,2] == "GNC",]
sampleTable <- sample.in
rownames(sampleTable) <- sampleTable[,1]
sampleTable <- sampleTable[colnames(tx.rsem$"counts"),]
colnames(sampleTable) <- c("x","ct","condition")

dds <- DESeqDataSetFromTximport(tx.rsem, sampleTable, ~ condition)
dds.re <- DESeq(dds,BPPARAM=MulticoreParam(12))
save(dds.re,file="GNC.gene.DEG.re")


load("GNC.gene.DEG.re")
DEG.re <- results(dds.re)
sig.DEG.re <- DEG.re[which(as.double(DEG.re[,"padj"]) < 0.05 & abs(as.double(DEG.re[,"log2FoldChange"])) > log2(1.5)),]

library(biomaRt)
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

total.genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','ensembl_transcript_id'), 
    filters = 'ensembl_gene_id', values = rownames(tx.rsem[[1]]), mart = ensembl)
DEGs <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol','ensembl_transcript_id'), filters = 'ensembl_gene_id', values = rownames(sig.DEG.re), mart = ensembl)

DEGs <- unique(DEGs[,1:2])
rownames(DEGs) <- DEGs[,1]
DEG.list <- cbind(DEGs[rownames(sig.DEG.re),2],rownames(sig.DEG.re),sig.DEG.re[,"pvalue"],sig.DEG.re[,"log2FoldChange"])
write.table(DEG.list,"./results_v2/sig.sig.DEG.re",sep='\t',quote=F,row.names=F,col.names=F)

