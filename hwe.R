# CONTENTS
# HWExactStats (hwe2_asian)                                    # snp.annot.v01.wlh.RData
# hwe analysis (hwe2_asian)
# HWExactStats (hwe2_asian) - X-chromosome                     # snp.annot.v02.wlh.RData
# hwe analysis (hwe2_asian) - X-chromosome

# Fetch genotypes for HWE test --- Chromosome 23               # XChr.geno.v01.wlh.RData
# Make Counts                                                  # XChr.geno.v02.wlh.RData
                                                               # XChr.M.count.v02.wlh.RData
                                                               # XChr.F.count.v02.wlh.RData
# HWExactStats (hwe2_asian) - X-chromosome - both gender       # snp.annot.v03.wlh.RData
# hwe analysis (hwe2_asian) - X-chromosome - both gender
# p-values in different tests

# Fetch genotypes for HWE test                                 # dfgeno.v01.wlh.RData

#################################################################
# HWExactStats (hwe2_asian)
#################################################################
library(GWASTools);  sessionInfo()$otherPkgs$GWASTools$Version  # "1.20.0"
library(QCpipeline); sessionInfo()$otherPkgs$QCpipeline$Version # "0.10.7"
library("HardyWeinberg"); sessionInfo()$otherPkgs$HardyWeinberg$Version  # "1.5.8"
date() # "Thu Jun  8 15:18:05 2017"


# load HWE asian results
machinePath <- "/projects"
filePath    <- "results/hwe2_asian"
fileName    <- "hwe.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
hwe.asian   <- get(load(fileIn)); dim(hwe.asian) # 513419      9
class(hwe.asian)    # "data.frame"
head(hwe.asian)
tail(hwe.asian)


table(hwe.asian$pval == "NA")
## FALSE
## 362884

# check if any monomorphics
table(hwe.asian$f %in% "NaN", useNA="ifany")
## FALSE   TRUE
## 362884 150535

# remove monomorphics from data
asian <- hwe.asian[!hwe.asian$f %in% "NaN",]; dim(asian)  # 362884      9
head(asian)

# make a matrix including counts of AA, AB and BB for each SNP
asian.nAA <- c(asian$nAA)
asian.nAB <- c(asian$nAB)
asian.nBB <- c(asian$nBB)
asian.count <- cbind(asian.nAA,asian.nAB,asian.nBB)
colnames(asian.count) <- c("nAA", "nAB", "nBB")
rownames(asian.count) <- asian$snpID
class(asian.count) # "matrix"
dim(asian.count)   # 362884      3
head(asian.count)

# use "HWExactStats" function from Jan's package:HardyWeinberg to compute Exact test or Mid p-values test
Exact.pvalues.asian <- HWExactStats(asian.count,x.linked=FALSE);  length(Exact.pvalues.asian) # 362884
MidP.pvalues.asian <- HWExactStats(asian.count,x.linked=FALSE, midp = TRUE); length(MidP.pvalues.asian) # 362884

class(Exact.pvalues.asian) # "numeric"
class(MidP.pvalues.asian)  # "numeric"

asian$exact.pval.asian <- Exact.pvalues.asian
asian$midp.pval.asian <- MidP.pvalues.asian

dim(asian)   # 362884     11
head(asian)

all(hwe.asian$snpID == asian$snpID) # FALSE


asi<-asian[,c("snpID","exact.pval.asian","midp.pval.asian")]
new.asian<-merge(hwe.asian,asi, all.x=TRUE); dim(new.asian)  # 513419     11
head(new.asian)

# read in snp table with project data
machinePath <- "/projects"
filePath    <- "sample_snp_annot"
fileName    <- "snp.annot.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
snp         <- get(load(fileIn)); dim(snp) # 514948        26
head(pData(snp))


# HWE contains chromosomes 1:23
all(snp$snpID[snp$chromosome %in% 1:23] == new.asian$snpID) # TRUE

# add HWE p-values into snp table
snp$exact.pval.asian[snp$chromosome %in% 1:23] <- new.asian$exact.pval.asian
snp$midp.pval.asian[snp$chromosome %in% 1:23] <- new.asian$midp.pval.asian


# update annotations
meta <- varMetadata(snp); dim(meta) # 28 1
head(meta,n=3)
idx <- which(is.na(meta$labelDescription)); length(idx) # 2
rownames(meta)[idx]    # [1] "exact.pval.asian" "midp.pval.asian"
meta["exact.pval.asian", "labelDescription"] <- "p-value from HardyWeinberg-HWExactStats exact test of Hardy-Weinberg Equilibrium in set of 352 hwe.asian.mcd=TRUE samples"
meta["midp.pval.asian", "labelDescription"] <- "p-value from HardyWeinberg-HWExactStats mid p-value test of Hardy-Weinberg Equilibrium in set of 352 hwe.asian.mcd=TRUE samples"
varMetadata(snp) <- meta

# save snp data frame
dim(snp) # 514948        28
all(sort(getVariable(snp,"snpID"))==getVariable(snp,"snpID")) #  TRUE
head(pData(snp))
saveas(snp, "snp.annot.v01.wlh.RData", "/projects/wlhsu/results/hwe/HardyWeinberg/asian")

rm(list=objects())








###################################################
# hwe analysis (hwe2_asian)
###################################################
library(GWASTools);  sessionInfo()$otherPkgs$GWASTools$Version  # "1.20.0"
library(QCpipeline); sessionInfo()$otherPkgs$QCpipeline$Version # "0.10.7"
date() # "Fri Jun  9 09:52:14 2017"
options(digits=7)

# sample table
machinePath <- "/projects"
filePath    <- "sample_snp_annot"
fileName    <- "sample.annot.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
sample      <- getobj(fileIn); dim(sample) # 30752    79
head(pData(sample),n=3)

# read in snp table with project data
machinePath <- "/projects"
filePath    <- "wlhsu/results/hwe/HardyWeinberg/asian"
fileName    <- "snp.annot.v01.wlh.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
snp         <- get(load(fileIn)); dim(snp) # 514948        28
head(pData(snp))

# look at p-values
summary(abs(snp$hwe.asian.pval-snp$exact.pval.asian))

summary(abs(snp$hwe.asian.pval-snp$midp.pval.asian))

summary(abs(snp$hwe.asian.pval-snp$exact.pval.asian.X))

summary(abs(snp$hwe.asian.pval-snp$midp.pval.asian.X))

# look at p-values -- typically we use cutoff of p < 1e-4,
summary(snp$hwe.asian.pval[snp$chromosome<=23])

summary(snp$exact.pval.asian[snp$chromosome<=23])

summary(snp$midp.pval.asian[snp$chromosome<=23])

# get count of HWE filters by MAF bin
res <- matrix(NA,nrow=4,ncol=4); dim(res)
colnames(res) <- c("pval","hwe.asian","exact.asian","midp.asian")
res[,1] <- c(1e-6,1e-5,1e-4,1e-3)

chk <- pData(snp)[!is.na(snp$hwe.asian.pval),]; dim(chk) # 362884     26
res[,2] <- c(sum(chk$hwe.asian.pval < 1e-6),
             sum(chk$hwe.asian.pval < 1e-5),
             sum(chk$hwe.asian.pval < 1e-4),
             sum(chk$hwe.asian.pval < 1e-3))

chk <- pData(snp)[!is.na(snp$exact.pval.asian),]; dim(chk) # 362884     28
res[,3] <- c(sum(chk$exact.pval.asian < 1e-6),
             sum(chk$exact.pval.asian < 1e-5),
             sum(chk$exact.pval.asian < 1e-4),
             sum(chk$exact.pval.asian < 1e-3))

chk <- pData(snp)[!is.na(snp$midp.pval.asian),]; dim(chk)  # 362884     28
res[,4] <- c(sum(chk$midp.pval.asian < 1e-6),
             sum(chk$midp.pval.asian < 1e-5),
             sum(chk$midp.pval.asian < 1e-4),
             sum(chk$midp.pval.asian < 1e-3))

res


# look at custom content
idx <- which(substr(snp$rsID,1,3) %in% "gco"); length(idx) # 15802
table(snp$chromosome[idx])

sum(snp$missing.n1[idx] %in% 1) # 2480
sum(snp$missing.n1[idx] %in% 1)/sum(snp$missing.n1 %in% 1) # 0.1799187
sum(substr(snp$rsID,1,3) %in% "gco")/nrow(snp) # 0.03068659
#
# so custom content is 3% of the SNPs, but 17% of the technical failures

custom <- matrix(NA,nrow=4,ncol=4); dim(custom)
colnames(custom) <- c("pval","hwe.asian","exact.asian","midp.asian")
custom[,1] <- c(1e-6,1e-5,1e-4,1e-3)

chk <- pData(snp)[!is.na(snp$hwe.asian.pval) & substr(snp$rsID,1,3) %in% "gco",]; dim(chk) # 3662   28
custom[,2] <- c(sum(chk$hwe.asian.pval < 1e-6),
                sum(chk$hwe.asian.pval < 1e-5),
                sum(chk$hwe.asian.pval < 1e-4),
                sum(chk$hwe.asian.pval < 1e-3))

chk <- pData(snp)[!is.na(snp$exact.pval.asian) & substr(snp$rsID,1,3) %in% "gco",]; dim(chk) # 3662   28
custom[,3] <- c(sum(chk$exact.pval.asian < 1e-6),
                sum(chk$exact.pval.asian < 1e-5),
                sum(chk$exact.pval.asian < 1e-4),
                sum(chk$exact.pval.asian < 1e-3))

chk <- pData(snp)[!is.na(snp$midp.pval.asian) & substr(snp$rsID,1,3) %in% "gco",]; dim(chk) # 3662   28
custom[,3] <- c(sum(chk$midp.pval.asian < 1e-6),
                sum(chk$midp.pval.asian < 1e-5),
                sum(chk$midp.pval.asian < 1e-4),
                sum(chk$midp.pval.asian < 1e-3))
custom


# qq plots
# (1) hwe.asian.pval
# check if any X chrom p-values are valid (in case of all-male study)
qq_outfile<-"/projects/wlhsu/plots/hwe/HardyWeinberg/asian/qq_plot_hwe.png"
plotX <- sum(!is.na(snp$hwe.asian.pval[snp$chromosome == 23])) > 0

if (plotX) nrow <- 2 else nrow <- 1
png(qq_outfile, width=720, height=(360*nrow))
par(mfrow=c(nrow,2), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
qqPlot(snp$hwe.asian.pval[snp$chromosome < 23],  trunc=FALSE, main="Autosomes, all")
qqPlot(snp$hwe.asian.pval[snp$chromosome < 23], trunc=TRUE, main="Autosomes, truncated")

if (plotX)
{
    qqPlot(snp$hwe.asian.pval[snp$chromosome == 23], trunc=FALSE, main="X chr, all")
    qqPlot(snp$hwe.asian.pval[snp$chromosome == 23], trunc=TRUE, main="X chromosome, truncated")
}
dev.off()


# (2) exact.pval.asian
# check if any X chrom p-values are valid (in case of all-male study)
qq_outfile<-"/projects/wlhsu/plots/hwe/HardyWeinberg/asian/qq_plot_exact.png"
plotX <- sum(!is.na(snp$exact.pval.asian[snp$chromosome == 23])) > 0

if (plotX) nrow <- 2 else nrow <- 1
png(qq_outfile, width=720, height=(360*nrow))
par(mfrow=c(nrow,2), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
qqPlot(snp$exact.pval.asian[snp$chromosome < 23],  trunc=FALSE, main="Autosomes, all")
qqPlot(snp$exact.pval.asian[snp$chromosome < 23], trunc=TRUE, main="Autosomes, truncated")

if (plotX)
{
    qqPlot(snp$exact.pval.asian[snp$chromosome == 23], trunc=FALSE, main="X chr, all")
    qqPlot(snp$exact.pval.asian[snp$chromosome == 23], trunc=TRUE, main="X chromosome, truncated")
}
dev.off()


# (3) midp.pval.asian
# check if any X chrom p-values are valid (in case of all-male study)
qq_outfile<-"/projects/wlhsu/plots/hwe/HardyWeinberg/asian/qq_plot_midp.png"
plotX <- sum(!is.na(snp$midp.pval.asian[snp$chromosome == 23])) > 0

if (plotX) nrow <- 2 else nrow <- 1
png(qq_outfile, width=720, height=(360*nrow))
par(mfrow=c(nrow,2), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
qqPlot(snp$midp.pval.asian[snp$chromosome < 23],  trunc=FALSE, main="Autosomes, all")
qqPlot(snp$midp.pval.asian[snp$chromosome < 23], trunc=TRUE, main="Autosomes, truncated")

if (plotX)
{
    qqPlot(snp$midp.pval.asian[snp$chromosome == 23], trunc=FALSE, main="X chr, all")
    qqPlot(snp$midp.pval.asian[snp$chromosome == 23], trunc=TRUE, main="X chromosome, truncated")
}
dev.off()






##########################################################################################
# remove monomorphics ? --- test

snpNew <- snp[1:5,]
head(pData(snpNew))

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asian/hwe_vs_exact.test.png",width=720, height=720)
plot(-log10(snpNew$hwe.asian.pval), -log10(snpNew$exact.pval.asian),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
     main="P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
dev.off()








##########################################################################################

# check bias
# (1) all chromosome
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asian/hwe_vs_exact.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval), -log10(snp$exact.pval.asian),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
     main="P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asian/hwe_vs_midp.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval), -log10(snp$midp.pval.asian),
     xlab="-log10(p-values) GWASTools HWE Exact Test", ylab="-log10(p-values) HardyWeinberg Mid-P Test",
     main="P-Values from GWASTools HWE Exact Test vs. HardyWeinberg Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()

# (2) Autosomes (chromosomes 1 ~ 22)
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asian/hwe_vs_exact_Autosomes.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome < 23]), -log10(snp$exact.pval.asian[snp$chromosome < 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
     main="For Autosomes, P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asian/hwe_vs_midp_Autosomes.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome < 23]), -log10(snp$midp.pval.asian[snp$chromosome < 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Mid-P Test",
     main="For Autosomes, P-Values from GWASTools HWE Exact Test vs. HardyWeinberg Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()

# (3) X chromosome --- only female [Y chromosome == 25]
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asian/hwe_vs_exact_XChr.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome == 23]), -log10(snp$exact.pval.asian[snp$chromosome == 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
     main="X chromosome, P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asian/hwe_vs_midp_XChr.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome == 23]), -log10(snp$midp.pval.asian[snp$chromosome == 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Mid-P Test",
     main="X chromosome, P-Values from GWASTools HWE Exact Test vs. HardyWeinberg Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()

# (4) Autosomes (chromosomes 1 ~ 22) --- truncated
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asian/hwe_vs_exact_Autosomes_truncated.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome < 23]), -log10(snp$exact.pval.asian[snp$chromosome < 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
     xlim=c(0,4), ylim=c(0,4),
     main="Autosomes - truncated, P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asian/hwe_vs_midp_Autosomes_truncated.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome < 23]), -log10(snp$midp.pval.asian[snp$chromosome < 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Mid-P Test",
     xlim=c(0,4), ylim=c(0,4),
main="Autosomes - truncated, P-Values from GWASTools HWE Exact Test vs. HardyWeinberg Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()


rm(list=objects())




#################################################################
# HWExactStats (hwe2_asian) - X-chromosome
#################################################################
library(GWASTools);  sessionInfo()$otherPkgs$GWASTools$Version  # "1.20.0"
library(QCpipeline); sessionInfo()$otherPkgs$QCpipeline$Version # "0.10.7"
library("HardyWeinberg"); sessionInfo()$otherPkgs$HardyWeinberg$Version  # "1.5.8"
date() # "Thu Jun  8 15:18:05 2017"


# load HWE asian results
machinePath <- "/projects"
filePath    <- "results/hwe2_asian"
fileName    <- "hwe.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
hwe.asian   <- get(load(fileIn)); dim(hwe.asian) # 513419      9
class(hwe.asian)    # "data.frame"
head(hwe.asian)
tail(hwe.asian)

table(hwe.asian$pval == "NA")
## FALSE
## 362884

# check if any monomorphics
table(hwe.asian$f %in% "NaN", useNA="ifany")
## FALSE   TRUE
## 362884 150535

# remove monomorphics from data
asian <- hwe.asian[!hwe.asian$f %in% "NaN",]; dim(asian)  # 362884      9
head(asian)


# seperate autosome and X chromosome
# (1) Autosome (chr 1- 22)
asian.A <- asian[asian$chr < 23,]; dim(asian.A) # 355193      9
tail(asian.A)


# make a matrix including counts of AA, AB and BB for each SNP
asian.A.nAA <- c(asian.A$nAA)
asian.A.nAB <- c(asian.A$nAB)
asian.A.nBB <- c(asian.A$nBB)
asian.A.count <- cbind(asian.A.nAA,asian.A.nAB,asian.A.nBB)
colnames(asian.A.count) <- c("nAA", "nAB", "nBB")
rownames(asian.A.count) <- asian.A$snpID
class(asian.A.count) # "matrix"
dim(asian.A.count)   # 355193     3
head(asian.A.count)

# use "HWExactStats" function from Jan's package:HardyWeinberg to compute Exact test or Mid p-values test
Exact.pvalues.asian.A <- HWExactStats(asian.A.count,x.linked=FALSE);  length(Exact.pvalues.asian.A) # 355193
MidP.pvalues.asian.A <- HWExactStats(asian.A.count,x.linked=FALSE, midp = TRUE); length(MidP.pvalues.asian.A) # 355193

# (2) X chromosome (chr 23)
asian.X <- asian[asian$chr == 23,]; dim(asian.X) # 7691    9
head(asian.X)


# make a matrix including counts of A, B, AA, AB and BB for each SNP
asian.X.nA <- matrix(0,nrow(asian.X),1); dim(asian.nA)  # 7691    1
asian.X.nB <- matrix(0,nrow(asian.X),1); dim(asian.nB)  # 7691    1
asian.X.nAA <- c(asian.X$nAA)
asian.X.nAB <- c(asian.X$nAB)
asian.X.nBB <- c(asian.X$nBB)
asian.X.count <- cbind(asian.X.nA, asian.X.nB, asian.X.nAA,asian.X.nAB,asian.X.nBB);
colnames(asian.X.count) <- c("A", "B", "nAA", "nAB", "nBB")
rownames(asian.X.count) <- asian.X$snpID
class(asian.X.count) # "matrix"
dim(asian.X.count)   # 7691    5
head(asian.X.count)

# use "HWExactStats" function from Jan's package:HardyWeinberg to compute Exact test or Mid p-values test
Exact.pvalues.asian.X <- HWExactStats(asian.X.count,x.linked=TRUE);  length(Exact.pvalues.asian.X) # 7691
MidP.pvalues.asian.X <- HWExactStats(asian.X.count,x.linked=TRUE, midp = TRUE); length(MidP.pvalues.asian.X) # 7691

# Combine P values of autosome and X chromosome
Exact.pvalues.asian <- c(Exact.pvalues.asian.A, Exact.pvalues.asian.X); length(Exact.pvalues.asian)# 362884
MidP.pvalues.asian <- c(MidP.pvalues.asian.A, MidP.pvalues.asian.X); length(MidP.pvalues.asian)    # 362884

asian$exact.pval.asian.X <- Exact.pvalues.asian
asian$midp.pval.asian.X <- MidP.pvalues.asian

summary(Exact.pvalues.asian.A)

summary(MidP.pvalues.asian.A)

summary(Exact.pvalues.asian.X)

summary(MidP.pvalues.asian.X)

summary(Exact.pvalues.asian)

summary(MidP.pvalues.asian)

dim(asian)   # 362884     11
head(asian)


all(hwe.asian$snpID == asian$snpID) # FALSE

asi<-asian[,c("snpID","exact.pval.asian.X","midp.pval.asian.X")]
new.asian<-merge(hwe.asian,asi, all.x=TRUE); dim(new.asian)  # 513419     11
head(new.asian)

# read in snp table with project data
machinePath <- "/projects"
filePath    <- "wlhsu/results/hwe/HardyWeinberg/asian"
fileName    <- "snp.annot.v01.wlh.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
snp         <- get(load(fileIn)); dim(snp) # 514948        28
head(pData(snp))


# HWE contains chromosomes 1:23
all(snp$snpID[snp$chromosome %in% 1:23] == new.asian$snpID) # TRUE

# add HWE p-values into snp table
snp$exact.pval.asian.X[snp$chromosome %in% 1:23] <- new.asian$exact.pval.asian.X
snp$midp.pval.asian.X[snp$chromosome %in% 1:23] <- new.asian$midp.pval.asian.X


# update annotations
meta <- varMetadata(snp); dim(meta) # 30 1
head(meta,n=3)
idx <- which(is.na(meta$labelDescription)); length(idx) # 2
rownames(meta)[idx]    # [1] "exact.pval.asian.X" "midp.pval.asian.X"
meta["exact.pval.asian.X", "labelDescription"] <- "p-value from HardyWeinberg-HWExactStats exact test of Hardy-Weinberg Equilibrium on the X-chromosome in set of 352 hwe.asian.mcd=TRUE samples"
meta["midp.pval.asian.X", "labelDescription"] <- "p-value from HardyWeinberg-HWExactStats mid p-value test of Hardy-Weinberg Equilibrium on the X-chromosome in set of 352 hwe.asian.mcd=TRUE samples"
varMetadata(snp) <- meta

# save snp data frame
dim(snp) # 514948        30
all(sort(getVariable(snp,"snpID"))==getVariable(snp,"snpID")) #  TRUE
head(pData(snp))
saveas(snp, "snp.annot.v02.wlh.RData", "/projects/wlhsu/results/hwe/HardyWeinberg/asian")

rm(list=objects())






###################################################
# hwe analysis (hwe2_asian) - X-chromosome
###################################################
library(GWASTools);  sessionInfo()$otherPkgs$GWASTools$Version  # "1.20.0"
library(QCpipeline); sessionInfo()$otherPkgs$QCpipeline$Version # "0.10.7"
date() # "Mon Jun 12 15:36:53 2017"


# sample table
machinePath <- "/projects"
filePath    <- "sample_snp_annot"
fileName    <- "sample.annot.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
sample      <- getobj(fileIn); dim(sample) # 30752    79
head(pData(sample),n=3)

# read in snp table with project data
machinePath <- "/projects"
filePath    <- "wlhsu/results/hwe/HardyWeinberg/asian"
fileName    <- "snp.annot.v02.wlh.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
snp         <- get(load(fileIn)); dim(snp) # 514948        30
head(pData(snp))


# look at p-values -- typically we use cutoff of p < 1e-4,
summary(snp$hwe.asian.pval[snp$chromosome<=23])

summary(snp$exact.pval.asian.X[snp$chromosome<=23])

summary(snp$midp.pval.asian.X[snp$chromosome<=23])

# get count of HWE filters by MAF bin
res <- matrix(NA,nrow=4,ncol=4); dim(res)
colnames(res) <- c("pval","hwe.asian","exact.asian","midp.asian")
res[,1] <- c(1e-6,1e-5,1e-4,1e-3)

chk <- pData(snp)[!is.na(snp$hwe.asian.pval),]; dim(chk) # 362884     26
res[,2] <- c(sum(chk$hwe.asian.pval < 1e-6),
sum(chk$hwe.asian.pval < 1e-5),
sum(chk$hwe.asian.pval < 1e-4),
sum(chk$hwe.asian.pval < 1e-3))

chk <- pData(snp)[!is.na(snp$exact.pval.asian.X),]; dim(chk) # 362884     28
res[,3] <- c(sum(chk$exact.pval.asian.X < 1e-6),
sum(chk$exact.pval.asian.X < 1e-5),
sum(chk$exact.pval.asian.X < 1e-4),
sum(chk$exact.pval.asian.X < 1e-3))

chk <- pData(snp)[!is.na(snp$midp.pval.asian.X),]; dim(chk)  # 362884     28
res[,4] <- c(sum(chk$midp.pval.asian.X < 1e-6),
sum(chk$midp.pval.asian.X < 1e-5),
sum(chk$midp.pval.asian.X < 1e-4),
sum(chk$midp.pval.asian.X < 1e-3))

res



# look at custom content
idx <- which(substr(snp$rsID,1,3) %in% "gco"); length(idx) # 15802
table(snp$chromosome[idx])

sum(snp$missing.n1[idx] %in% 1) # 2480
sum(snp$missing.n1[idx] %in% 1)/sum(snp$missing.n1 %in% 1) # 0.1799187
sum(substr(snp$rsID,1,3) %in% "gco")/nrow(snp) # 0.03068659


custom <- matrix(NA,nrow=4,ncol=4); dim(custom)
colnames(custom) <- c("pval","hwe.asian","exact.asian","midp.asian")
custom[,1] <- c(1e-6,1e-5,1e-4,1e-3)

chk <- pData(snp)[!is.na(snp$hwe.asian.pval) & substr(snp$rsID,1,3) %in% "gco",]; dim(chk) # 3662   28
custom[,2] <- c(sum(chk$hwe.asian.pval < 1e-6),
sum(chk$hwe.asian.pval < 1e-5),
sum(chk$hwe.asian.pval < 1e-4),
sum(chk$hwe.asian.pval < 1e-3))

chk <- pData(snp)[!is.na(snp$exact.pval.asian.X) & substr(snp$rsID,1,3) %in% "gco",]; dim(chk) # 3662   28
custom[,3] <- c(sum(chk$exact.pval.asian.X < 1e-6),
sum(chk$exact.pval.asian.X < 1e-5),
sum(chk$exact.pval.asian.X < 1e-4),
sum(chk$exact.pval.asian.X < 1e-3))

chk <- pData(snp)[!is.na(snp$midp.pval.asian.X) & substr(snp$rsID,1,3) %in% "gco",]; dim(chk) # 3662   28
custom[,3] <- c(sum(chk$midp.pval.asian.X < 1e-6),
sum(chk$midp.pval.asian.X < 1e-5),
sum(chk$midp.pval.asian.X < 1e-4),
sum(chk$midp.pval.asian.X < 1e-3))
custom




# qq plots
# (1) hwe.asian.pval
# check if any X chrom p-values are valid (in case of all-male study)
qq_outfile<-"/projects/wlhsu/plots/hwe/HardyWeinberg/asianX/qq_plot_hwe.png"
plotX <- sum(!is.na(snp$hwe.asian.pval[snp$chromosome == 23])) > 0

if (plotX) nrow <- 2 else nrow <- 1
png(qq_outfile, width=720, height=(360*nrow))
par(mfrow=c(nrow,2), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
qqPlot(snp$hwe.asian.pval[snp$chromosome < 23],  trunc=FALSE, main="Autosomes, all")
qqPlot(snp$hwe.asian.pval[snp$chromosome < 23], trunc=TRUE, main="Autosomes, truncated")

if (plotX)
{
    qqPlot(snp$hwe.asian.pval[snp$chromosome == 23], trunc=FALSE, main="X chr, all")
    qqPlot(snp$hwe.asian.pval[snp$chromosome == 23], trunc=TRUE, main="X chromosome, truncated")
}
dev.off()


# (2) exact.pval.asian.X
# check if any X chrom p-values are valid (in case of all-male study)
qq_outfile<-"/projects/wlhsu/plots/hwe/HardyWeinberg/asianX/qq_plot_exact.png"
plotX <- sum(!is.na(snp$exact.pval.asian.X[snp$chromosome == 23])) > 0

if (plotX) nrow <- 2 else nrow <- 1
png(qq_outfile, width=720, height=(360*nrow))
par(mfrow=c(nrow,2), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
qqPlot(snp$exact.pval.asian.X[snp$chromosome < 23],  trunc=FALSE, main="Autosomes, all")
qqPlot(snp$exact.pval.asian.X[snp$chromosome < 23], trunc=TRUE, main="Autosomes, truncated")

if (plotX)
{
    qqPlot(snp$exact.pval.asian.X[snp$chromosome == 23], trunc=FALSE, main="X chr, all")
    qqPlot(snp$exact.pval.asian.X[snp$chromosome == 23], trunc=TRUE, main="X chromosome, truncated")
}
dev.off()


# (3) midp.pval.asian.X
# check if any X chrom p-values are valid (in case of all-male study)
qq_outfile<-"/projects/wlhsu/plots/hwe/HardyWeinberg/asianX/qq_plot_midp.png"
plotX <- sum(!is.na(snp$midp.pval.asian.X[snp$chromosome == 23])) > 0

if (plotX) nrow <- 2 else nrow <- 1
png(qq_outfile, width=720, height=(360*nrow))
par(mfrow=c(nrow,2), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
qqPlot(snp$midp.pval.asian.X[snp$chromosome < 23],  trunc=FALSE, main="Autosomes, all")
qqPlot(snp$midp.pval.asian.X[snp$chromosome < 23], trunc=TRUE, main="Autosomes, truncated")

if (plotX)
{
    qqPlot(snp$midp.pval.asian.X[snp$chromosome == 23], trunc=FALSE, main="X chr, all")
    qqPlot(snp$midp.pval.asian.X[snp$chromosome == 23], trunc=TRUE, main="X chromosome, truncated")
}
dev.off()


# check bias
# (1) all chromosome
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianX/hwe_vs_exact.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval), -log10(snp$exact.pval.asian.X),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
     main="P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianX/hwe_vs_midp.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval), -log10(snp$midp.pval.asian.X),
     xlab="-log10(p-values) GWASTools HWE Exact Test", ylab="-log10(p-values) HardyWeinberg Mid-P Test",
     main="P-Values from GWASTools HWE Exact Test vs. HardyWeinberg Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()

# (2) Autosomes (chromosomes 1 ~ 22)
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianX/hwe_vs_exact_Autosomes.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome < 23]), -log10(snp$exact.pval.asian.X[snp$chromosome < 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
     main="For Autosomes, P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianX/hwe_vs_midp_Autosomes.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome < 23]), -log10(snp$midp.pval.asian.X[snp$chromosome < 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Mid-P Test",
     main="For Autosomes, P-Values from GWASTools HWE Exact Test vs. HardyWeinberg Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()

# (3) X chromosome --- only female [Y chromosome == 25]
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianX/hwe_vs_exact_XChr.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome == 23]), -log10(snp$exact.pval.asian.X[snp$chromosome == 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
     main="X chromosome, P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianX/hwe_vs_midp_XChr.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome == 23]), -log10(snp$midp.pval.asian.X[snp$chromosome == 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Mid-P Test",
     main="X chromosome, P-Values from GWASTools HWE Exact Test vs. HardyWeinberg Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()

# (4) Autosomes (chromosomes 1 ~ 22) --- truncated
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianX/hwe_vs_exact_Autosomes_truncated.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome < 23]), -log10(snp$exact.pval.asian.X[snp$chromosome < 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
     xlim=c(0,4), ylim=c(0,4),
     main="Autosomes - truncated, P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianX/hwe_vs_midp_Autosomes_truncated.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome < 23]), -log10(snp$midp.pval.asian.X[snp$chromosome < 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Mid-P Test",
     xlim=c(0,4), ylim=c(0,4),
     main="Autosomes - truncated, P-Values from GWASTools HWE Exact Test vs. HardyWeinberg Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()



rm(list=objects())








####### runRscript: /projects/wlhsu/src/hwe/hwe_geno.R   ###############
###################################################
# Fetch genotypes for HWE test --- Chromosome 23
###################################################
library(GWASTools);  sessionInfo()$otherPkgs$GWASTools$Version  # "1.19.2"
library(QCpipeline); sessionInfo()$otherPkgs$QCpipeline$Version # "0.10.6"
date()


# read imputed sample table
machinePath <- "/projects"
filePath    <- "sample_snp_annot"
fileName    <- "sample.annot.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
sample      <- get(load(fileIn)); dim(sample) # 30752        79
head(pData(sample),n=3)

# chr 23 snp table - common to both dosage and best-guess genotype
machinePath <- "/projects"
filePath    <- "sample_snp_annot"
fileName    <- "snp.annot.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
snp.chr1    <- get(load(fileIn)); dim(snp.chr1) # 514948        26
head(pData(snp.chr1),n=3)

# chr 23 genotype GDS file
machinePath      <- "/projects"
filePath         <- "gds/samples"
fileName         <- "geno.gds"
gdsFile.geno.1  <- paste(machinePath,filePath,fileName,sep="/")


X <- snp.chr1[snp.chr1$chromosome == 23,]; dim(X)  # 14210        26
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
all(names(geno)==sample$scanID) # TRUE
table(geno,useNA="ifany")

# get all with a loop
geno <- matrix(NA,nrow=nrow(sample),ncol=length(snpID)); dim(geno)   # 30752 14210
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


df <- cbind(pData(sample),geno); dim(df) # 30752    82
head(df)

saveas(df, "XChr.geno.v01.wlh.RData", "/projects/wlhsu/results/hwe/HardyWeinberg")


rm(list=objects())







###################################################
# Make Counts
###################################################
library(GWASTools);  sessionInfo()$otherPkgs$GWASTools$Version
library(QCpipeline); sessionInfo()$otherPkgs$QCpipeline$Version # "0.10.6"
date() # "Fri Jun 16 14:07:30 2017"


# read imputed sample table
machinePath <- "/projects"
filePath    <- "wlhsu/results/hwe/HardyWeinberg"
fileName    <- "XChr.geno.v01.wlh.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
snp         <- get(load(fileIn)); dim(snp) # 30752 14289
head(snp,n=3)
names(snp)

table(snp$hwe.asian.mcd)
## FALSE  TRUE
## 30400   352

# get those 352 individuals in "hwe.asian.mcd"  only
asian <- snp[snp$hwe.asian.mcd,]; dim(asian)  # 352 14289

# remove SNP with all NA with variables "scanID" and "sex"
rNA.asi <- cbind(asian$scanID,asian$sex,Filter(function(x) !all(is.na(x)), asian[,80:ncol(asian)])); dim(rNA.asi)  # 352 13037
head(rNA.asi,2)

colnames(rNA.asi)[colnames(rNA.asi)=="asian$scanID"] <- "scanID"
colnames(rNA.asi)[colnames(rNA.asi)=="asian$sex"] <- "sex"
head(rNA.asi)

# save data with "scanID", "sex" and SNP without all NA
saveas(rNA.asi, "XChr.geno.v02.wlh.RData", "/projects/wlhsu/results/hwe/HardyWeinberg")

# check how many males and femals
table(rNA.asi$sex)
##   F   M
## 131 221

# (1) get those 221 males only
asian.M <- rNA.asi[rNA.asi$sex %in% "M",]; dim(asian.M)  # 221 13037
M <- asian.M[,3:ncol(asian.M)]; dim(M) # 221 13035
head(M,2)

# (2) get those 131 females only
asian.F <- rNA.asi[rNA.asi$sex %in% "F",]; dim(asian.F)  # 131 13037
F <- asian.F[,3:ncol(asian.F)]; dim(F) # 131 13035
head(F,2)


# Make Counts
# (1) males
M <- as.matrix(M) # coerce data frames to matrices.
n <- nrow(M)
p <- ncol(M)
C <- matrix(NA,nrow=p,ncol=4)
for (j in 1:p) {
    snp <- M[, j]
    nB  <- sum(snp == 0, na.rm = TRUE)
    nAB <- sum(snp == 1, na.rm = TRUE)
    nA  <- sum(snp == 2, na.rm = TRUE)
    nNA <- sum(is.na(snp))
    tot <- nA + nAB + nB + nNA
    if (tot != n) {
        cat(j, "\n")
        stop("genotypes and missings do not sum n")
    }
    C[j,] <- c(nA, nAB, nB, nNA)
}

head(C)

snpID = colnames(M)
colnames(C) <- c("nA", "nAB", "nB", "nNA")
cM <- data.frame(snpID, C); dim(cM)   #  13035     5
head(cM)


saveas(cM, "XChr.M.count.v02.wlh.RData", "/projects/wlhsu/results/hwe/HardyWeinberg")

# (2) females
F <- as.matrix(F) # coerce data frames to matrices.
n <- nrow(F)
p <- ncol(F)
C <- matrix(NA,nrow=p,ncol=4)
for (j in 1:p) {
    snp <- F[, j]
    nBB <- sum(snp == 0, na.rm = TRUE)
    nAB <- sum(snp == 1, na.rm = TRUE)
    nAA <- sum(snp == 2, na.rm = TRUE)
    nNA <- sum(is.na(snp))
    tot <- nAA + nAB + nBB + nNA
    if (tot != n) {
        cat(j, "\n")
        stop("genotypes and missings do not sum n")
    }
    C[j,] <- c(nAA, nAB, nBB, nNA)
}

head(C)

snpID = colnames(F)
colnames(C) <- c("nAA", "nAB", "nBB", "nNA")
cF <- data.frame(snpID, C); dim(cF)   #  7 5
head(cF)

saveas(cF, "XChr.F.count.v02.wlh.RData", "/projects/wlhsu/results/hwe/HardyWeinberg")


rm(list=objects())








#################################################################
# HWExactStats (hwe2_asian) - X-chromosome - both gender
#################################################################
library(GWASTools);  sessionInfo()$otherPkgs$GWASTools$Version  # "1.20.0"
library(QCpipeline); sessionInfo()$otherPkgs$QCpipeline$Version # "0.10.7"
library("HardyWeinberg"); sessionInfo()$otherPkgs$HardyWeinberg$Version  # "1.5.8"
date()  # "Fri Jun 16 14:14:32 2017"

# load genotype counts on X chromosome for males
machinePath <- "/projects"
filePath    <- "wlhsu/results/hwe/HardyWeinberg"
fileName    <- "XChr.M.count.v02.wlh.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
ChX.M.asian <- get(load(fileIn)); dim(ChX.M.asian) #  13035     5
head(ChX.M.asian)

# load genotype counts on X chromosome for females
machinePath <- "/projects"
filePath    <- "wlhsu/results/hwe/HardyWeinberg"
fileName    <- "XChr.F.count.v02.wlh.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
ChX.F.asian <- get(load(fileIn)); dim(ChX.F.asian) # 13035     5
head(ChX.F.asian)

# load HWE asian results
machinePath <- "/projects"
filePath    <- "results/hwe2_asian"
fileName    <- "hwe.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
hwe.asian   <- get(load(fileIn)); dim(hwe.asian) # 513419      9
head(hwe.asian)
tail(hwe.asian)


table(hwe.asian$pval == "NA")
## FALSE
## 362884

# check if any monomorphics
table(hwe.asian$f %in% "NaN", useNA="ifany")
## FALSE   TRUE
## 362884 150535

# remove monomorphics from data
asian <- hwe.asian[!hwe.asian$f %in% "NaN",]; dim(asian)  # 362884      9
head(asian)

# seperate autosome and X chromosome
# (1) Autosome (chr 1- 22)
asian.A <- asian[asian$chr < 23,]; dim(asian.A) # 355193      9
tail(asian.A)


# make a matrix including counts of AA, AB and BB for each SNP
asian.A.nAA <- c(asian.A$nAA)
asian.A.nAB <- c(asian.A$nAB)
asian.A.nBB <- c(asian.A$nBB)
asian.A.count <- cbind(asian.A.nAA,asian.A.nAB,asian.A.nBB)
colnames(asian.A.count) <- c("nAA", "nAB", "nBB")
rownames(asian.A.count) <- asian.A$snpID
dim(asian.A.count)   # 355193     3
head(asian.A.count)

# use "HWExactStats" function from Jan's package:HardyWeinberg to compute Exact test or Mid p-values test
Exact.pvalues.asian.A <- HWExactStats(asian.A.count,x.linked=FALSE);  length(Exact.pvalues.asian.A) # 355193
MidP.pvalues.asian.A <- HWExactStats(asian.A.count,x.linked=FALSE, midp = TRUE); length(MidP.pvalues.asian.A) # 355193


# (2) X chromosome (chr 23)
asian.X <- asian[asian$chr == 23,]; dim(asian.X) # 7691    9
head(asian.X)


all(ChX.M.asian$snpID == ChX.F.asian$snpID) # TRUE

ChX.M <- subset(ChX.M.asian, select=c("snpID", "nA", "nB")); head(ChX.M,3)
ChX.F <- subset(ChX.F.asian, select=c("snpID", "nAA", "nAB", "nBB")); head(ChX.F,3)

ChX   <- merge(ChX.M,ChX.F, by = "snpID"); dim(ChX) # 13035     6
head(ChX,3)

all(ChX$snpID == asian.X$snpID)    # FALSE

XID <- asian.X[c(1)]; head(XID)

X <- merge(XID, ChX, all.x=TRUE); dim(X)  # 7691    6
Xc <- na.omit(X); dim(Xc)    # 7691    6
head(Xc)

asian.X.count <- Xc[,2:ncol(Xc)]; dim(asian.X.count)   # 7691    5
head(asian.X.count)

# use "HWExactStats" function from Jan's package:HardyWeinberg to compute Exact test or Mid p-values test
Exact.pvalues.asian.X <- HWExactStats(asian.X.count,x.linked=TRUE);  length(Exact.pvalues.asian.X) # 7691
MidP.pvalues.asian.X <- HWExactStats(asian.X.count,x.linked=TRUE, midp = TRUE); length(MidP.pvalues.asian.X) # 7691

# Combine P values of autosome and X chromosome
Exact.pvalues.asian <- c(Exact.pvalues.asian.A, Exact.pvalues.asian.X); length(Exact.pvalues.asian)# 362884
MidP.pvalues.asian <- c(MidP.pvalues.asian.A, MidP.pvalues.asian.X); length(MidP.pvalues.asian)    # 362884

# check p-value
summary(Exact.pvalues.asian.A)
summary(MidP.pvalues.asian.A)
summary(Exact.pvalues.asian.X)
summary(MidP.pvalues.asian.X)
summary(Exact.pvalues.asian)
summary(MidP.pvalues.asian)

asian$exact.pval.asian.X.MF <- Exact.pvalues.asian
asian$midp.pval.asian.X.MF <- MidP.pvalues.asian

dim(asian)   # 362884     11
head(asian)

all(hwe.asian$snpID == asian$snpID) # FALSE

asi<-asian[,c("snpID","exact.pval.asian.X.MF","midp.pval.asian.X.MF")]
new.asian<-merge(hwe.asian,asi, all.x=TRUE); dim(new.asian)  # 513419     11
head(new.asian)

# read in snp table with project data
machinePath <- "/projects"
filePath    <- "wlhsu/results/hwe/HardyWeinberg/asian"
fileName    <- "snp.annot.v02.wlh.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
snp         <- get(load(fileIn)); dim(snp) # 514948        30
head(pData(snp))
varLabels(snp)

# HWE contains chromosomes 1:23
all(snp$snpID[snp$chromosome %in% 1:23] == new.asian$snpID) # TRUE

# add HWE p-values into snp table
snp$exact.pval.asian.X.MF[snp$chromosome %in% 1:23] <- new.asian$exact.pval.asian.X.MF
snp$midp.pval.asian.X.MF[snp$chromosome %in% 1:23] <- new.asian$midp.pval.asian.X.MF

# update annotations
meta <- varMetadata(snp); dim(meta) # 32 1
head(meta,n=3)
idx <- which(is.na(meta$labelDescription)); length(idx) # 2
rownames(meta)[idx]    # [1] "exact.pval.asian.X.MF" "midp.pval.asian.X.MF"
meta["exact.pval.asian.X.MF", "labelDescription"] <- "p-value from HardyWeinberg-HWExactStats exact test of Hardy-Weinberg Equilibrium on the X-chromosome with both male and female in set of 352 hwe.asian.mcd=TRUE samples"
meta["midp.pval.asian.X.MF", "labelDescription"] <- "p-value from HardyWeinberg-HWExactStats mid p-value test of Hardy-Weinberg Equilibrium on the X-chromosome with both male and female in set of 352 hwe.asian.mcd=TRUE samples"
varMetadata(snp) <- meta

# save snp data frame
dim(snp) # 514948        32
all(sort(getVariable(snp,"snpID"))==getVariable(snp,"snpID")) #  TRUE
head(pData(snp))
saveas(snp, "snp.annot.v03.wlh.RData", "/projects/wlhsu/results/hwe/HardyWeinberg/asian")

rm(list=objects())



##########################################################
# hwe analysis (hwe2_asian) - X-chromosome - both gender
##########################################################
library(GWASTools);  sessionInfo()$otherPkgs$GWASTools$Version  # "1.20.0"
library(QCpipeline); sessionInfo()$otherPkgs$QCpipeline$Version # "0.10.7"
date() #  "Fri Jun 16 16:19:04 2017"


# sample table
machinePath <- "/projects"
filePath    <- "sample_snp_annot"
fileName    <- "sample.annot.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
sample      <- getobj(fileIn); dim(sample) # 30752    79
head(pData(sample),n=3)

# read in snp table with project data
machinePath <- "/projects"
filePath    <- "wlhsu/results/hwe/HardyWeinberg/asian"
fileName    <- "snp.annot.v03.wlh.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
snp         <- get(load(fileIn)); dim(snp) # 514948        30
head(pData(snp))


# look at p-values -- typically we use cutoff of p < 1e-4,
summary(snp$hwe.asian.pval[snp$chromosome<=23])
summary(snp$exact.pval.asian.X.MF[snp$chromosome<=23])
summary(snp$midp.pval.asian.X.MF[snp$chromosome<=23])

# get count of HWE filters by MAF bin
res <- matrix(NA,nrow=4,ncol=4); dim(res)
colnames(res) <- c("pval","hwe.asian","exact.asian","midp.asian")
res[,1] <- c(1e-6,1e-5,1e-4,1e-3)

chk <- pData(snp)[!is.na(snp$hwe.asian.pval),]; dim(chk) # 362884     26
res[,2] <- c(sum(chk$hwe.asian.pval < 1e-6),
sum(chk$hwe.asian.pval < 1e-5),
sum(chk$hwe.asian.pval < 1e-4),
sum(chk$hwe.asian.pval < 1e-3))

chk <- pData(snp)[!is.na(snp$exact.pval.asian.X.MF),]; dim(chk) # 362884     28
res[,3] <- c(sum(chk$exact.pval.asian.X.MF < 1e-6),
sum(chk$exact.pval.asian.X.MF < 1e-5),
sum(chk$exact.pval.asian.X.MF < 1e-4),
sum(chk$exact.pval.asian.X.MF < 1e-3))

chk <- pData(snp)[!is.na(snp$midp.pval.asian.X.MF),]; dim(chk)  # 362884     28
res[,4] <- c(sum(chk$midp.pval.asian.X.MF < 1e-6),
sum(chk$midp.pval.asian.X.MF < 1e-5),
sum(chk$midp.pval.asian.X.MF < 1e-4),
sum(chk$midp.pval.asian.X.MF < 1e-3))

res

# look at custom content
idx <- which(substr(snp$rsID,1,3) %in% "gco"); length(idx) # 15802
table(snp$chromosome[idx])

sum(snp$missing.n1[idx] %in% 1) # 2480
sum(snp$missing.n1[idx] %in% 1)/sum(snp$missing.n1 %in% 1) # 0.1799187
sum(substr(snp$rsID,1,3) %in% "gco")/nrow(snp) # 0.03068659
#
# so custom content is 3% of the SNPs, but 17% of the technical failures

custom <- matrix(NA,nrow=4,ncol=4); dim(custom)
colnames(custom) <- c("pval","hwe.asian","exact.asian","midp.asian")
custom[,1] <- c(1e-6,1e-5,1e-4,1e-3)

chk <- pData(snp)[!is.na(snp$hwe.asian.pval) & substr(snp$rsID,1,3) %in% "gco",]; dim(chk) # 3662   28
custom[,2] <- c(sum(chk$hwe.asian.pval < 1e-6),
sum(chk$hwe.asian.pval < 1e-5),
sum(chk$hwe.asian.pval < 1e-4),
sum(chk$hwe.asian.pval < 1e-3))

chk <- pData(snp)[!is.na(snp$exact.pval.asian.X.MF) & substr(snp$rsID,1,3) %in% "gco",]; dim(chk) # 3662   28
custom[,3] <- c(sum(chk$exact.pval.asian.X.MF < 1e-6),
sum(chk$exact.pval.asian.X.MF < 1e-5),
sum(chk$exact.pval.asian.X.MF < 1e-4),
sum(chk$exact.pval.asian.X.MF < 1e-3))

chk <- pData(snp)[!is.na(snp$midp.pval.asian.X.MF) & substr(snp$rsID,1,3) %in% "gco",]; dim(chk) # 3662   28
custom[,3] <- c(sum(chk$midp.pval.asian.X.MF < 1e-6),
sum(chk$midp.pval.asian.X.MF < 1e-5),
sum(chk$midp.pval.asian.X.MF < 1e-4),
sum(chk$midp.pval.asian.X.MF < 1e-3))
custom


# qq plots
# (1) hwe.asian.pval
# check if any X chrom p-values are valid (in case of all-male study)
qq_outfile<-"/projects/wlhsu/plots/hwe/HardyWeinberg/asianXMF/qq_plot_hwe.png"
plotX <- sum(!is.na(snp$hwe.asian.pval[snp$chromosome == 23])) > 0

if (plotX) nrow <- 2 else nrow <- 1
png(qq_outfile, width=720, height=(360*nrow))
par(mfrow=c(nrow,2), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
qqPlot(snp$hwe.asian.pval[snp$chromosome < 23],  trunc=FALSE, main="Autosomes, all")
qqPlot(snp$hwe.asian.pval[snp$chromosome < 23], trunc=TRUE, main="Autosomes, truncated")

if (plotX)
{
    qqPlot(snp$hwe.asian.pval[snp$chromosome == 23], trunc=FALSE, main="X chr, all")
    qqPlot(snp$hwe.asian.pval[snp$chromosome == 23], trunc=TRUE, main="X chromosome, truncated")
}
dev.off()


# (2) exact.pval.asian.X
# check if any X chrom p-values are valid (in case of all-male study)
qq_outfile<-"/projects/wlhsu/plots/hwe/HardyWeinberg/asianXMF/qq_plot_exact.png"
plotX <- sum(!is.na(snp$exact.pval.asian.X.MF[snp$chromosome == 23])) > 0

if (plotX) nrow <- 2 else nrow <- 1
png(qq_outfile, width=720, height=(360*nrow))
par(mfrow=c(nrow,2), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
qqPlot(snp$exact.pval.asian.X.MF[snp$chromosome < 23],  trunc=FALSE, main="Autosomes, all")
qqPlot(snp$exact.pval.asian.X.MF[snp$chromosome < 23], trunc=TRUE, main="Autosomes, truncated")

if (plotX)
{
    qqPlot(snp$exact.pval.asian.X.MF[snp$chromosome == 23], trunc=FALSE, main="X chr, all")
    qqPlot(snp$exact.pval.asian.X.MF[snp$chromosome == 23], trunc=TRUE, main="X chromosome, truncated")
}
dev.off()


# (3) midp.pval.asian.X
# check if any X chrom p-values are valid (in case of all-male study)
qq_outfile<-"/projects/wlhsu/plots/hwe/HardyWeinberg/asianXMF/qq_plot_midp.png"
plotX <- sum(!is.na(snp$midp.pval.asian.X.MF[snp$chromosome == 23])) > 0

if (plotX) nrow <- 2 else nrow <- 1
png(qq_outfile, width=720, height=(360*nrow))
par(mfrow=c(nrow,2), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
qqPlot(snp$midp.pval.asian.X.MF[snp$chromosome < 23],  trunc=FALSE, main="Autosomes, all")
qqPlot(snp$midp.pval.asian.X.MF[snp$chromosome < 23], trunc=TRUE, main="Autosomes, truncated")

if (plotX)
{
    qqPlot(snp$midp.pval.asian.X.MF[snp$chromosome == 23], trunc=FALSE, main="X chr, all")
    qqPlot(snp$midp.pval.asian.X.MF[snp$chromosome == 23], trunc=TRUE, main="X chromosome, truncated")
}
dev.off()


# check bias
# (1) all chromosome
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianXMF/hwe_vs_exact.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval), -log10(snp$exact.pval.asian.X.MF),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
     main="P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianXMF/hwe_vs_midp.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval), -log10(snp$midp.pval.asian.X.MF),
     xlab="-log10(p-values) GWASTools HWE Exact Test", ylab="-log10(p-values) HardyWeinberg Mid-P Test",
     main="P-Values from GWASTools HWE Exact Test vs. HardyWeinberg Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()

# (2) Autosomes (chromosomes 1 ~ 22)
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianXMF/hwe_vs_exact_Autosomes.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome < 23]), -log10(snp$exact.pval.asian.X.MF[snp$chromosome < 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
     main="For Autosomes, P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianXMF/hwe_vs_midp_Autosomes.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome < 23]), -log10(snp$midp.pval.asian.X.MF[snp$chromosome < 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Mid-P Test",
     main="For Autosomes, P-Values from GWASTools HWE Exact Test vs. HardyWeinberg Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()

# (3) X chromosome --- both males and females
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianXMF/hwe_vs_exact_XChr_MF.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome == 23]), -log10(snp$exact.pval.asian.X.MF[snp$chromosome == 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
     main="X chromosome, P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianXMF/hwe_vs_midp_XChr_MF.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome == 23]), -log10(snp$midp.pval.asian.X.MF[snp$chromosome == 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Mid-P Test",
     main="X chromosome, P-Values from GWASTools HWE Exact Test vs. HardyWeinberg Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()


# (4) Autosomes (chromosomes 1 ~ 22) --- truncated
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianXMF/hwe_vs_exact_Autosomes_truncated.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome < 23]), -log10(snp$exact.pval.asian.X.MF[snp$chromosome < 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
     xlim=c(0,4), ylim=c(0,4),
     main="Autosomes - truncated, P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianXMF/hwe_vs_midp_Autosomes_truncated.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome < 23]), -log10(snp$midp.pval.asian.X.MF[snp$chromosome < 23]),
     xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Mid-P Test",
     xlim=c(0,4), ylim=c(0,4),
     main="Autosomes - truncated, P-Values from GWASTools HWE Exact Test vs. HardyWeinberg Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()



rm(list=objects())

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianXMF/test.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome == 23]), -log10(snp$exact.pval.asian.X.MF[snp$chromosome == 23]),
xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Exact Test",
main="X chromosome, P-Values from HWE Exact Test: GWASTools vs. HardyWeinberg")
abline (a=0, b=1, col="red")
abline(h = 0:4, v = 0:4, col = "lightgray", lty=3)
abline(h = 1.5, v = 1.5, col = "blue", lty=4)
dev.off()











#############################################
# p-values in different tests
#############################################
# read in snp table with project data
machinePath <- "/projects"
filePath    <- "wlhsu/results/hwe/HardyWeinberg/asian"
fileName    <- "snp.annot.v03.wlh.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
snp         <- get(load(fileIn)); dim(snp) # 514948        30
head(pData(snp))


# look at p-values
# (1) chromosome 1-23
summary(snp$hwe.asian.pval[snp$chromosome<=23])

summary(snp$exact.pval.asian.X[snp$chromosome<=23])

summary(snp$midp.pval.asian.X[snp$chromosome<=23])

summary(snp$exact.pval.asian.X.MF[snp$chromosome<=23])

summary(snp$midp.pval.asian.X.MF[snp$chromosome<=23])


# (2)chromosome 1-22
summary(snp$hwe.asian.pval[snp$chromosome < 23])

summary(snp$exact.pval.asian[snp$chromosome < 23])

summary(snp$midp.pval.asian[snp$chromosome < 23])


# (3)chromosome 23 = X chromosome
summary(snp$hwe.asian.pval[snp$chromosome == 23])

summary(snp$exact.pval.asian.X[snp$chromosome == 23])

summary(snp$midp.pval.asian.X[snp$chromosome == 23])

summary(snp$exact.pval.asian.X.MF[snp$chromosome == 23])

summary(snp$midp.pval.asian.X.MF[snp$chromosome == 23])


# Different of p-values between tests
# (1)chromosome 1-22
summary(abs(snp$hwe.asian.pval[snp$chromosome < 23]-snp$exact.pval.asian[snp$chromosome < 23]))

summary(abs(snp$hwe.asian.pval[snp$chromosome < 23]-snp$midp.pval.asian[snp$chromosome < 23]))


# (2)chromosome 23 = X chromosome
summary(abs(snp$hwe.asian.pval[snp$chromosome == 23]-snp$exact.pval.asian[snp$chromosome == 23]))

summary(abs(snp$hwe.asian.pval[snp$chromosome == 23]-snp$midp.pval.asian[snp$chromosome == 23]))

summary(abs(snp$hwe.asian.pval[snp$chromosome == 23]-snp$exact.pval.asian.X[snp$chromosome == 23]))

summary(abs(snp$hwe.asian.pval[snp$chromosome == 23]-snp$midp.pval.asian.X[snp$chromosome == 23]))

summary(abs(snp$hwe.asian.pval[snp$chromosome == 23]-snp$exact.pval.asian.X.MF[snp$chromosome == 23]))

summary(abs(snp$hwe.asian.pval[snp$chromosome == 23]-snp$midp.pval.asian.X.MF[snp$chromosome == 23]))




##  Chromosome	    SNP
##  1		     35,796
##  2			 42,593
##  3			 34,493
##  4			 27,737
##  5			 31,580
##  6			 39,880
##  7			 25,527
##  8			 26,826
##  9			 21,483
##  10			 25,641
##  11			 24,051
##  12			 25,683
##  13			 16,703
##  14			 15,918
##  15			 14,791
##  16			 15,599
##  17			 17,537
##  18			 14,535
##  19			 13,056
##  20			 14,133
##  21			  6,207
##  22			  9,440
##  23			 14,210
##  total	    513,419














####### runRscript: /projects/wlhsu/src/hwe/hwe_geno.R   ###############
###################################################
# Fetch genotypes for HWE test
###################################################
library(GWASTools);  sessionInfo()$otherPkgs$GWASTools$Version  # "1.19.2"
library(QCpipeline); sessionInfo()$otherPkgs$QCpipeline$Version # "0.10.6"
date() # "Wed Aug 17 17:13:58 2016"


# read imputed sample table
machinePath <- "/projects"
filePath    <- "sample_snp_annot"
fileName    <- "sample.annot.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
sample      <- get(load(fileIn)); dim(sample) # 30752        79
head(pData(sample),n=3)

# chr 1 snp table - common to both dosage and best-guess genotype
machinePath <- "/projects"
filePath    <- "sample_snp_annot"
fileName    <- "snp.annot.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
snp.chr1    <- get(load(fileIn)); dim(snp.chr1) # 514948        26
head(pData(snp.chr1),n=3)

# chr 1 genotype GDS file
machinePath      <- "/projects"
filePath         <- "gds/samples"
fileName         <- "geno.gds"
gdsFile.geno.1  <- paste(machinePath,filePath,fileName,sep="/")


snpID <- c(snp.chr1$snpID)
head(pData(snp.chr1)[snp.chr1$snpID %in% snpID,c("rsID","snpID")])

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
all(names(geno)==sample$scanID) # TRUE
table(geno,useNA="ifany")

# get all with a loop
geno <- matrix(NA,nrow=nrow(sample),ncol=length(snpID)); dim(geno)
colnames(geno) <- paste(snpID)
#colnames(geno) <- paste("geno",snpID,sep=".")
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
head(df)

saveas(df, "dfgeno.v01.wlh.RData", "/projects/wlhsu/results/hwe/HardyWeinberg")


rm(list=objects())



