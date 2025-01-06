library(GenomicFeatures)
library(VariantAnnotation)
library(doParallel)
library(BiocParallel)
library(VariantAnnotation)
library(data.table)

range.f <- function(sig.re,cn){
  min.v <- apply(sig.re[,cn],1,function(x)  min(as.integer(x)))
  max.v <- apply(sig.re[,cn],1,function(x)  max(as.integer(x)))
  g.ran <- GRanges(sig.re[,4],ranges=IRanges(min.v,max.v))
}

sQTLs.test <- function(test.mat,test.g.ran,alt.tp){
  over.vcf <- as.matrix(findOverlaps(test.g.ran,vcf.geno.ran))
  F.sqtl.re <- NULL
  for (i in 1:length(test.mat[,1])){
    Ay.exp <- rbind(as.double(strsplit(test.mat[i,"IncLevel1"],",")[[1]]))
    Eu.exp <- rbind(as.double(strsplit(test.mat[i,"IncLevel2"],",")[[1]]))
    colnames(Ay.exp) <- as.matrix(sample.info[sample.info[,1] == "Ay",3])
    colnames(Eu.exp) <- as.matrix(sample.info[sample.info[,1] == "Eu",3])
    ov.Ay <- intersect(colnames(Ay.exp),rownames(rna.sample.info))
    ov.Eu <- intersect(colnames(Eu.exp),rownames(rna.sample.info))
    Ay.exp <- Ay.exp[,ov.Ay]
    Eu.exp <- Eu.exp[,ov.Eu]
    names(Ay.exp) <- as.matrix(rna.sample.info[ov.Ay,2])
    names(Eu.exp) <- as.matrix(rna.sample.info[ov.Eu,2])
    total.exp <- cbind(names(c(Ay.exp,Eu.exp)),c(Ay.exp,Eu.exp))
    colnames(total.exp) <- c("sample","exp")
    ea.total.vcf <- rbind(final.geno.mat[over.vcf[over.vcf[,1] == i,2],])
    if (!length(ea.total.vcf))  next
    for (j in 1:length(ea.total.vcf[,1])){
      ea.geno <- cbind(colnames(ea.total.vcf)[grep("^J",colnames(ea.total.vcf))],ea.total.vcf[j,grep("^J",colnames(ea.total.vcf))])
      real.num <- length(which(ea.geno[,2] != "./."))
      eu.real.num <- length(which(ea.geno[names(Eu.exp),2] != "./."))
      if (real.num >= 7 & eu.real.num == 3){
        colnames(ea.geno) <- c("sample","geno")
        merged.data <- merge(total.exp,ea.geno,by.x="sample",by.y="sample")
        merged.data <- merged.data[merged.data[,3] != "./.",]
        ran.p <- NULL
        lm.re <- lm(as.double(as.matrix(merged.data[,2])) ~ as.matrix(merged.data[,3]))
        lm.p <- summary(lm.re)$coefficient[2,"Pr(>|t|)"]
	      lm.or <- summary(lm.re)$coefficient[2,"Estimate"]
        F.sqtl.re <- rbind(F.sqtl.re,c(test.mat[i,c(1,3)],alt.tp,lm.p,lm.or,ea.total.vcf[j,]))
      }
    }
  }
  return (F.sqtl.re)
}

# sQTLs analysis
A3SS.g.ran <- range.f(sig.GNC.A3SS.re,6:11)
A5SS.g.ran <- range.f(sig.GNC.A5SS.re[,4],6:11)
ES.g.ran <- range.f(sig.GNC.ES.re[,4],6:11)
IR.g.ran <- range.f(sig.GNC.IR.re[,4],6:11)
MXE.g.ran <- range.f(sig.GNC.MXE.re[,4],6:13)

sam.in <- read.csv("./5.rMATs.result/GNC/summary.txt",sep='\t',row.names=NULL)
sam.in <- gsub(".*GNC/|/Alig.*","",sam.in[grep("4.STAR.result",sam.in[,2]),2])

impute.vcf.geno <- read.table("./13.IMPUTE/result/vcffiles/total.vcf",sep='\t')
impute.vcf.mat <- gsub(" ","",as.matrix(impute.vcf.geno))

sample.info <- read.table("./rMATs.sample.info.GNC",sep='\t')
total.sample.info <- read.table("./GNC_WGS_RNA",sep='\t')
rna.sample.info <- total.sample.info
wgs.sample.info <- total.sample.info
rownames(rna.sample.info) <- rna.sample.info[,3]
rownames(wgs.sample.info) <- wgs.sample.info[,1]
ov.wgs.sample <- intersect(paste("11090X",1:19,sep=""),rownames(wgs.sample.info))

final.geno.mat <- impute.vcf.mat[,-c(6:9)]
colnames(final.geno.mat) <- c("chr","locus","id","ref","alt",paste("11090X",1:19,sep=""))
vcf.geno.ran <- GRanges(Rle(impute.vcf.mat[,1]),IRanges(as.integer(impute.vcf.mat[,2]),as.integer(impute.vcf.mat[,2])))


sQTLs.A3SS <- sQTLs.test(sig.GNC.A3SS.re,A3SS.g.ran,"A3SS")
sQTLs.A5SS <- sQTLs.test(sig.GNC.A5SS.re,A5SS.g.ran,"A5SS")
sQTLs.ES <- sQTLs.test(sig.GNC.ES.re,ES.g.ran,"ES")
sQTLs.IR <- sQTLs.test(sig.GNC.IR.re,IR.g.ran,"IR")
sQTLs.MXE <- sQTLs.test(sig.GNC.MXE.re,MXE.g.ran,"MXE")

total.sQTLs <- rbind(sQTLs.A3SS,sQTLs.A5SS,sQTLs.ES,sQTLs.IR,sQTLs.MXE)
sig.total.sQTLs <- total.sQTLs[which(as.double(total.sQTLs[,4]) < 0.05),]

# eQTLs analysis
txdb <- loadDb("~/txDB_GTF75")
fi.cns <- c("TXCHROM","TXNAME","GENEID","TXSTART","TXEND","TXSTRAND")
exon.re <- exonsBy(txdb,by="tx")
txTable <- gsub(" ","",as.matrix(select(txdb,keys=names(exon.re),columns=fi.cns,keytype="TXID")))
promoter.re <- promoters(txdb,upstream=1500,downstream=200)

ay.sam <- total.sample.info[is.element(total.sample.info[,3],sample.in[sample.in[,3] == "Ay",1]),1]
eu.sam <- total.sample.info[is.element(total.sample.info[,3],sample.in[sample.in[,3] == "Eu",1]),1]










