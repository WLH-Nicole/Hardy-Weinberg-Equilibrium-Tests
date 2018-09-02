#############
# Joint test for Hardy-Weinberg proportions in males and females --- Chromosome 22 _noEpsi
#############
# 6/20/2017

# qsub -N hwe_Chr22 -m e -M wlhsu@uw.edu -q bigmem.q /projects/QCpipeline/runRscript.sh /projects/wlhsu/src/hwe/jointHW/hwe_joint_Chr22_noEpsi.R

########################################

library(GWASTools);  sessionInfo()$otherPkgs$GWASTools$Version  # "1.19.2"
library(QCpipeline); sessionInfo()$otherPkgs$QCpipeline$Version # "0.10.6"
date()


# load genotype counts on X chromosome for males
machinePath <- "/projects"
filePath    <- "wlhsu/results/hwe/HardyWeinberg/asian/Autosome"
fileName    <- "Chr22.M.count.v02.wlh.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
cM <- get(load(fileIn)); dim(cM) #
head(cM)

# load genotype counts on X chromosome for females
machinePath <- "/projects"
filePath    <- "wlhsu/results/hwe/HardyWeinberg/asian/Autosome"
fileName    <- "Chr22.F.count.v02.wlh.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
cF <- get(load(fileIn)); dim(cF)
head(cF)

# load HWE asian results
machinePath <- "/projects"
filePath    <- "results/hwe2_asian"
fileName    <- "hwe.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
hwe.asian   <- get(load(fileIn)); dim(hwe.asian)
head(hwe.asian)

# remove monomorphics from data
asian <- hwe.asian[!hwe.asian$f %in% "NaN",]; dim(asian)
head(asian)

# seperate male and female for chromosome 22
asian.Ch <- asian[asian$chr == 22,]; dim(asian.Ch)
head(asian.Ch)

Ch.M <- subset(cM, select=c("snpID", "nAA", "nAB", "nBB")); head(Ch.M,3)
Ch.F <- subset(cF, select=c("snpID", "nAA", "nAB", "nBB")); head(Ch.F,3)

all(Ch.M$snpID == asian.Ch$snpID) # FALSE
all(Ch.F$snpID == asian.Ch$snpID) # FALSE
all(cM$snpID == cF$snpID)         # TRUE


# get snpIDs used for analysis
ID <- asian.Ch[c(1)]; head(ID)

# (1) male
ID.M <- merge(ID, Ch.M, all.x=TRUE); dim(ID.M)
IDc.M <- na.omit(ID.M); dim(IDc.M)
head(IDc.M)

# (2) female
ID.F <- merge(ID, Ch.F, all.x=TRUE); dim(ID.F)
IDc.F <- na.omit(ID.F); dim(IDc.F)
head(IDc.F)

all(IDc.M$snpID == IDc.F$snpID)         # TRUE

# (1) get counts only for male
M <- IDc.M[,2:ncol(IDc.M)]; dim(M)
head(M)

# (2) get counts only for female
F <- IDc.F[,2:ncol(IDc.F)]; dim(F)
head(F)

####### function codes from Jan Graffelman, June 20, 2017  ##########
# funcations for Joint test
outcomes <- function(x) {
    n <- sum(x)
    nm <- sum(x[1:3])
    nf <- n-nm
    nt.m <- 2*nm
    nt.f <- 2*nf
    nA <- 2*x[1]+x[2]+2*x[4]+x[5]
    nB <- 2*n - nA
    nA <- min(nA,nB)
    nAm.min <- max(0,nA-nt.f)
    nAm.max <- min(nA,nt.m)
    nam.seq <- seq(nAm.min,nAm.max)
    naf.seq <- nA-nam.seq
    nbm.seq <- nt.m-nam.seq
    nbf.seq <- nt.f-naf.seq
    m.het.seq <- pmin(nam.seq,nbm.seq)
    f.het.seq <- pmin(naf.seq,nbf.seq)
    nam.out <- floor(m.het.seq/2)+1
    naf.out <- floor(f.het.seq/2)+1
    tot <- sum(nam.out*naf.out)
    return(tot)
}

generate.outcomes <- function(x) {
    nm <- sum(x[1:3])
    nt.m <- 2*nm
    nf <- sum(x[4:6])
    nt.f <- 2*nf
    nA <- 2*x[1]+x[2]+2*x[4]+x[5]
    nout <- outcomes(x)
    out <- matrix(NA,nrow=nout,ncol=6)
    nam.min <- max(0,nA-2*nf)
    nam.max <- min(nA,2*nm)
    naf.min <- max(0,nA-2*nm)
    naf.max <- min(nA,2*nf)
    i <- 1
    for(nam in nam.max:nam.min) {
        
        nbm <- nt.m-nam
        het.max.m <- min(nam,nbm)
        
        naf <- nA-nam
        nbf <- nt.f-naf
        het.max.f <- min(naf,nbf)
        
        for(het.i in seq(het.max.m,0,-2)) {
            maa <- (nam-het.i)/2
            mbb <- nm-(het.i+maa)
            
            for(j in seq(het.max.f,0,-2)) {
                faa <- (naf-j)/2
                fbb <- nf-(j+faa)
                out[i,] <- c(maa,het.i,mbb,faa,j,fbb)
                i <- i+1
            }
        }
    }
    rownames(out) <- 1:nrow(out)
    return(out)
}

calcprob <- function(x,num,n) { # format x maa; mab; mbb; faa; fab; fbb;
    nab <- x[2]+x[5]
    den <- sum(lgamma(x+1)) + lgamma(2*n + 1)
    quo <- num - den + nab*log(2)
    prob <- exp(quo)
    return(prob)
}

pvalue.joint.autosomal.exact <- function(x,verbose=TRUE,midpvalue=TRUE) {
    #    epsilon <- 2^(-44)
    epsilon <- 0
    nm <- sum(x[1:3])
    nf <- sum(x[4:6])
    n <- nm+nf
    na <- 2*x[1]+x[2]+2*x[4]+x[5]
    nb <- 2*n - na
    num <- lgamma(na+1) + lgamma(nb+1) + lgamma(nf+1) + lgamma(nm+1)
    psamp <- calcprob(x,num,n)
    X <- generate.outcomes(x)
    nab <- X[,2]+X[,5]
    lg2n <- lgamma(2*n+1)
    den <- rowSums(lgamma(X+1)) + lg2n # vector
    quo <- num - den + nab*log(2)
    pvector <- exp(quo)
    psum <- sum(pvector)
    plessoe <- sum(pvector[pvector <= psamp + epsilon])
    midp <- plessoe - 0.5*psamp
    #    if(verbose) {
    #    cat("Autosomal exact test for Hardy-Weinberg proportions and equality of allele frequencies in the sexes.\n")
    #    if(midpvalue) {
    #        cat("Mid p-value =",midp,"\n")
    #    } else {
    #        cat("Standard p-value =",plessoe,"\n")
    #    }
    #}
    
    return(list(plessoe=plessoe,midp=midp,psum=psum))
}

# need to change datafram to matrix for joint test
MM <- as.matrix(M); head(MM)
FF <- as.matrix(F); head(FF)

# modify Jan's code to get p values from all SNPs from both exact test and midp test at the same time
out <- data.frame()
pval <- data.frame()

for(i in 1:nrow(MM)) {
    z <- c(MM[i,],FF[i,])
    out <- pvalue.joint.autosomal.exact(z,midpvalue = TRUE)
    pval <- rbind(pval,data.frame(out))
}

rownames(pval) <- c(IDc.M$snpID)
dim(pval)
head(pval)


#idpval <- data.frame(snpID, pval); dim(idpval)
#pval$snpID <- c(IDc.M$snpID)

# check p-value
summary(pval$plessoe)  # exact test


summary(pval$midp)     # midp test


asian.Ch$exact.pv.asian.jmf.ch22 <- pval$plessoe
asian.Ch$midp.pv.asian.jmf.ch22 <- pval$midp

dim(asian.Ch)
head(asian.Ch)

saveas(asian.Ch, "Chr22.jmf.noepsi.pv.wlh.RData", "/projects/wlhsu/results/hwe/HardyWeinberg/asian/Autosome")


all(hwe.asian$snpID == asian.Ch$snpID) # FALSE

asi<-asian.Ch[,c("snpID","exact.pv.asian.jmf.ch22","midp.pv.asian.jmf.ch22")]
new.asian<-merge(hwe.asian,asi, all.x=TRUE); dim(new.asian)  # 513419     11
head(new.asian)

# read in snp table with project data from fhcrc
machinePath <- "/projects"
filePath    <- "wlhsu/results/hwe/HardyWeinberg/asian"
fileName    <- "snp.annot.v03.wlh.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
snp         <- get(load(fileIn)); dim(snp) # 514948        32
head(pData(snp))
varLabels(snp)

# HWE contains chromosomes 1:23
all(snp$snpID[snp$chromosome %in% 1:23] == new.asian$snpID) # TRUE

# add HWE p-values into snp table
snp$exact.pv.asian.jmf.ch22[snp$chromosome %in% 1:23] <- new.asian$exact.pv.asian.jmf.ch22
snp$midp.pv.asian.jmf.ch22[snp$chromosome %in% 1:23] <- new.asian$midp.pv.asian.jmf.ch22

# update annotations
meta <- varMetadata(snp); dim(meta) # 34 1
head(meta,n=3)
idx <- which(is.na(meta$labelDescription)); length(idx) # 2
rownames(meta)[idx]    
meta["exact.pv.asian.jmf.ch22", "labelDescription"] <- "p-value from joint exact test without epsilon for Hardy-Weinberg proportions in males and females on chromosome 22 in set of 352 hwe.asian.mcd=TRUE samples"
meta["midp.pv.asian.jmf.ch22", "labelDescription"] <- "p-value from joint mid-p test without epsilon for Hardy-Weinberg proportions in males and females on chromosome 22 in set of 352 hwe.asian.mcd=TRUE samples"
varMetadata(snp) <- meta

# save snp data frame
dim(snp) # 514948        34
all(sort(getVariable(snp,"snpID"))==getVariable(snp,"snpID")) #  TRUE
head(pData(snp))
saveas(snp, "snp.annot.v04.noepsi.ch22.wlh.RData", "/projects/wlhsu/results/hwe/HardyWeinberg/asian/Autosome")


rm(list=objects())


