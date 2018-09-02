#############
# Fetch genotypes for HWE test --- Chromosome 22
#############
# 6/19/2017


# qsub -N hwe_Chr22 -m e -M wlhsu@uw.edu -q r420.q /projects/QCpipeline/runRscript.sh /projects/wlhsu/src/hwe/hwe_geno_Chr22.R

########################################

library(GWASTools);  sessionInfo()$otherPkgs$GWASTools$Version
library(QCpipeline); sessionInfo()$otherPkgs$QCpipeline$Version
date()


# read imputed sample table
machinePath <- "/projects"
filePath    <- "sample_snp_annot"
fileName    <- "sample.annot.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
sample      <- get(load(fileIn)); dim(sample)
head(pData(sample),n=3)

# chr 1 snp table - common to both dosage and best-guess genotype
machinePath <- "/projects"
filePath    <- "sample_snp_annot"
fileName    <- "snp.annot.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
snp.chr1    <- get(load(fileIn)); dim(snp.chr1)
head(pData(snp.chr1),n=3)

# chr 1 genotype GDS file
machinePath      <- "/projects"
filePath         <- "gds/samples"
fileName         <- "geno.gds"
gdsFile.geno.1  <- paste(machinePath,filePath,fileName,sep="/")


X <- snp.chr1[snp.chr1$chromosome == 22,]; dim(X)     
head(pData(X),n=3)


snpID <- c(X$snpID)
head(pData(snp.chr1)[snp.chr1$snpID %in% snpID,c("rsID","snpID","chromosome")])


GetGenotypeWrapper <- function(snpID,sample,snp,gdsFile){
    stopifnot(is.finite(snpID) & length(snpID)==1)
    gds      <- GdsGenotypeReader(gdsFile)
    genoData <- GenotypeData(gds,scanAnnot=sample,snpAnnot=snp)
    res      <- getGenotypeSelection(genoData, snpID=snpID, char=FALSE, sort=TRUE,drop=TRUE, use.names=TRUE)
    close(genoData)
    return(res)
}

# check the first one manually
geno <- GetGenotypeWrapper(snpID[1],sample,snp.chr1,gdsFile.geno.1)
all(names(geno)==sample$scanID)
table(geno,useNA="ifany")


# get all with a loop
geno <- matrix(NA,nrow=nrow(sample),ncol=length(snpID)); dim(geno)
colnames(geno) <- paste(snpID)
for(i in 1:length(snpID)) {
    res.geno <- GetGenotypeWrapper(snpID[i],sample,snp.chr1,gdsFile.geno.1)
    stopifnot(all(names(res.geno)==sample$scanID))
    if(i==1) {
        rownames(geno) <- names(res.geno)
    }
    geno[,i] <- res.geno
}
geno[1:6,1:6]


df <- cbind(pData(sample),geno); dim(df)
head(df, 2)

saveas(df, "Chr22.geno.v01.wlh.RData", "/projects/cidr/Peters_rep/wlhsu/results/hwe/HardyWeinberg/asian/Autosome")





