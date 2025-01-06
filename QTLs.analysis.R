library(GenomicFeatures)
library(VariantAnnotation)
library(doParallel)
library(BiocParallel)
library(VariantAnnotation)
library(data.table)

# sQTLs analysis
range.f <- function(sig.re,cn){
	min.v <- apply(sig.re[,cn],1,function(x){
		min(as.integer(x))})
	max.v <- apply(sig.re[,cn],1,function(x){
		max(as.integer(x))})
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
				lm.re <- lm(as.double(as.matrix(merged.data[,2])) ~ as.matrix(merged.data[,3]))
				lm.p <- summary(lm.re)$coefficient[2,"Pr(>|t|)"]
				lm.or <- summary(lm.re)$coefficient[2,"Estimate"]
				F.sqtl.re <- rbind(F.sqtl.re,c(test.mat[i,c(1,3)],alt.tp,lm.p,lm.or,ea.total.vcf[j,]))
				}
			}
		}
	return (F.sqtl.re)
	}

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
write.table(sig.total.sQTLs,"./results/sig.sQTLs.re",sep='\t',quote=F,row.names=F,col.names=F)

		 
# eQTLs analysis
txdb <- loadDb("~/txDB_GTF75")
fi.cns <- c("TXCHROM","TXNAME","GENEID","TXSTART","TXEND","TXSTRAND")
exon.re <- exonsBy(txdb,by="tx")
txTable <- gsub(" ","",as.matrix(select(txdb,keys=names(exon.re),columns=fi.cns,keytype="TXID")))
promoter.re <- promoters(txdb,upstream=1500,downstream=200)

sample.in <- gsub(" ","",as.matrix(read.table("./sample.info",sep='\t')))
sample.in <- sample.in[sample.in[,2] == "GNC",]
ay.sam <- total.sample.info[is.element(total.sample.info[,3],sample.in[sample.in[,3] == "Ay",1]),1]
eu.sam <- total.sample.info[is.element(total.sample.info[,3],sample.in[sample.in[,3] == "Eu",1]),1]
vcf.f <- "./13.IMPUTE/result/vcffiles/total.vcf.gz"

normalized_counts <- counts(dds.re, normalized=TRUE)
sig.exp <- normalized_counts[rownames(sig.DEG.re),]

f.eqtl.re <- NULL
for (i in 1:length(sig.exp[,1])){
	ea.exp <- sig.exp[i,]
	ea.pro.r <- promoter.re[rownames(sig.exp)[i],]
	ea.p <- ScanVcfParam(which=ea.pro.r)
	ea.vcf <- readVcf(TabixFile(vcf.f),"hg19",ea.p)
	ea.geno <- rbind(geno(ea.vcf)$GT)
	colnames(ea.geno) <- paste("11090X",1:19,sep="")
	row.al <- do.call(rbind,strsplit(gsub(".*_","",rownames(ea.geno)),"/"))
	ev1 <- is.element(row.al[,1],c("A","T","G","C")) & is.element(row.al[,2],c("A","T","G","C"))
	ev2 <- cbind(apply(ea.geno,1,function(x){
		c(length(which(x[ay.sam] != "./.")),length(which(x[eu.sam] != "./.")))
		}))
	ev2 <- ev2[1,] >= 7 & ev2[2,] >= 3
    	ev <- ev1 & ev2
    	if (!length(which(ev)))    next
    	te.geno <- rbind(ea.geno[ev,])
    	rownames(te.geno) <- rownames(ea.geno)[ev]
    	pre.re <- NULL
    	names(ea.exp) <- total.sample.info[names(ea.exp),1]
    	for (j in 1:length(te.geno[,1])){
		snp.pos <- strsplit(gsub("_.*","",rownames(te.geno)),":")[[j]]
        	alleles <- strsplit(gsub(".*_","",rownames(te.geno)[j]),"/")[[1]]
        	geno.m <- cbind(ids=colnames(te.geno),g=c(te.geno[j,]))
        	exp.m <- cbind(ids=names(ea.exp),e=ea.exp)
        	mer.mat <- merge(geno.m,exp.m,by.x="ids",by.y="ids")
        	ea.g <- do.call(rbind,strsplit(mer.mat[,2],"/"))
        	mer.mat[,2] <- as.integer(ea.g[,1]) + as.integer(ea.g[,2])
        	mer.mat[,3] <- as.double(mer.mat[,3])
        	mer.mat <- mer.mat[!is.na(mer.mat[,2]),]
        	lm.re <- try(lm(e ~ g,data=mer.mat),silent=T)
        	if (length(grep("Error",lm.re)))    next
		if (is.na(lm.re$coefficient[2]))    next
        	pv <- summary(lm.re)$coefficient[2,"Pr(>|t|)"]
		or <- summary(lm.re)$coefficient[2,"Estimate"]
        	pre.re <- rbind(pre.re,c(rownames(sig.exp)[i],snp.pos,pv,or,alleles))
		}
	f.eqtl.re <- rbind(f.eqtl.re,pre.re)
	}

sig.eqtl.re <- f.eqtl.re[which(as.double(f.eqtl.re[,4]) < 0.05),]
write.table(sig.eqtl.re,"./results/sig.eQTLs.re",sep='\t',quote=F,row.names=F)




