# Estimation for alternative splicing rates as PSI

comp <- commandArgs(TRUE)
comp <- comp[1]

oper.dir <- "4.STAR.result/"
gtf.path <- "/uufs/chpc.utah.edu/common/home/u0693520/2.reference.RNA.seq/ref/gtf/Homo_sapiens.GRCh37.75.gtf"
sample.info <- as.matrix(read.table("./Aymara/sample.info",sep='\t'))
out.dir <- paste("5.rMATs.result/",comp,"/",sep="")

if (comp == "GNC"){
    GNC.Ay <- sample.info[sample.info[,2] == "GNC" & sample.info[,3] == "Ay",]
    GNC.Eu <- sample.info[sample.info[,2] == "GNC" & sample.info[,3] == "Eu",]
    b1 <- paste(paste("4.STAR.result/GNC/",GNC.Ay[,1],"/Aligned.sortedByCoord.out.bam",sep=""),collapse=",")
    b2 <- paste(paste("4.STAR.result/GNC/",GNC.Eu[,1],"/Aligned.sortedByCoord.out.bam",sep=""),collapse=",")
}

rmats.oper <- paste("python ./Tools/rMATS.3.2.5/RNASeq-MATS.py -b1 ",b1," -b2 ",b2," -t paired -len 125 -gtf ",gtf.path," -novelSS 0 -analysis U -o ",out.dir,sep="")
system(rmats.oper)
