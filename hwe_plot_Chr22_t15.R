#############
# Joint test for Hardy-Weinberg proportions in males and females --- Chromosome 22 - t15
#############
# 6/20/2017

# qsub -N hwe_plot_Chr22_t15 -m e -M wlhsu@uw.edu -q bigmem.q /projects/QCpipeline/runRscript.sh /projects/wlhsu/src/hwe/jointHW/Chplot/hwe_plot_Chr22_t15.R

########################################

library(GWASTools);  sessionInfo()$otherPkgs$GWASTools$Version  # "1.19.2"
library(QCpipeline); sessionInfo()$otherPkgs$QCpipeline$Version # "0.10.6"
date()

# read in snp table with project data from fhcrc
machinePath <- "/projects/cidr"
filePath    <- "Peters_rep/wlhsu/results/hwe/HardyWeinberg/asian/Autosome"
fileName    <- "Peters_rep.snp.annot.v04.t15.ch22.wlh.RData"
fileIn      <- paste(machinePath,filePath,fileName,sep="/")
snp         <- get(load(fileIn)); dim(snp) # 514948        34
varLabels(snp)


# look at p-values -- typically we use cutoff of p < 1e-4,
summary(snp$hwe.asian.pval[snp$chromosome == 22])
summary(snp$exact.pv.asian.jmf.ch22[snp$chromosome == 22])
summary(snp$midp.pv.asian.jmf.ch22[snp$chromosome == 22])


# qq plots

qq_outfile<-"/projects/wlhsu/plots/hwe/HardyWeinberg/asianAMF/ch22/qq_plot_hwe_Ch22_joint_t15.png"
png(qq_outfile, width=720, height=(360*3))
par(mfrow=c(3,2), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
qqPlot(snp$hwe.asian.pval[snp$chromosome == 22],  trunc=FALSE, main="GWASTools Exact Test: Chr 22, all")
qqPlot(snp$hwe.asian.pval[snp$chromosome == 22], trunc=TRUE, main="GWASTools Exact Test: Chr 22, truncated")
qqPlot(snp$exact.pv.asian.jmf.ch22[snp$chromosome == 22], ylim=c(0,40), trunc=FALSE, main="HardyWeinberg Joint Exact Test: Chr 22, all")
qqPlot(snp$exact.pv.asian.jmf.ch22[snp$chromosome == 22], trunc=TRUE, main="HardyWeinberg Joint Exact Test: Chr 22, truncated")
qqPlot(snp$midp.pv.asian.jmf.ch22[snp$chromosome == 22], ylim=c(0,40), trunc=FALSE, main="HardyWeinberg Joint Mid-P Test: Chr 22, all")
qqPlot(snp$midp.pv.asian.jmf.ch22[snp$chromosome == 22], trunc=TRUE, main="HardyWeinberg Joint Mid-P Test: Chr 22, truncated")
dev.off()


# check bias

# (1) chromosome 22 - all
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianAMF/ch22/hwe_vs_exact_Ch22_joint_t15.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome == 22]), -log10(snp$exact.pv.asian.jmf.ch22[snp$chromosome == 22]),
xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Joint Exact Test",
xlim=c(0,40), ylim=c(0,40),
main="Chr 22, P-Values from GWASTools Exact Test vs. HardyWeinberg Joint Exact Test")
abline (a=0, b=1, col="red")
#abline(h = 4, v = 4, col = "blue", lty=4)
segments(0, 4, 4, 4, col = "blue", lty=2)
segments(4, 0, 4, 4, col = "blue", lty=2)
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianAMF/ch22/hwe_vs_midp_Ch22_joint_t15.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome == 22]), -log10(snp$midp.pv.asian.jmf.ch22[snp$chromosome == 22]),
xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Joint Mid-P Test",
xlim=c(0,40), ylim=c(0,40),
main="Chr 22, P-Values from GWASTools Exact Test vs. HardyWeinberg Joint Mid-P Test")
abline (a=0, b=1, col="red")
#abline(h = 4, v = 4, col = "blue", lty=4)
segments(0, 4, 4, 4, col = "blue", lty=2)
segments(4, 0, 4, 4, col = "blue", lty=2)
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianAMF/ch22/exact_vs_midp_Ch22_joint_t15.png",width=720, height=720)
plot(-log10(snp$exact.pv.asian.jmf.ch22[snp$chromosome == 22]), -log10(snp$midp.pv.asian.jmf.ch22[snp$chromosome == 22]),
xlab="-log10(p-values) HardyWeinberg Joint Exact Test", ylab="-log10(p-values) HardyWeinberg Joint Mid-P Test",
xlim=c(0,40), ylim=c(0,40),
main="Chr 22, P-Values from HardyWeinberg: Joint Exact Test vs. Joint Mid-P Test")
abline (a=0, b=1, col="red")
#abline(h = 4, v = 4, col = "blue", lty=4)
segments(0, 4, 4, 4, col = "blue", lty=2)
segments(4, 0, 4, 4, col = "blue", lty=2)
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianAMF/ch22/midp_vs_midp_Ch22_joint_t15.png",width=720, height=720)
plot(-log10(snp$midp.pval.asian[snp$chromosome == 22]), -log10(snp$midp.pv.asian.jmf.ch22[snp$chromosome == 22]),
xlab="-log10(p-values) HardyWeinberg Mid-P Test", ylab="-log10(p-values) HardyWeinberg Joint Mid-P Test",
xlim=c(0,40), ylim=c(0,40),
main="Chr 22, P-Values from HardyWeinberg: Mid-P Test vs. Joint Mid-P Test")
abline (a=0, b=1, col="red")
#abline(h = 4, v = 4, col = "blue", lty=4)
segments(0, 4, 4, 4, col = "blue", lty=2)
segments(4, 0, 4, 4, col = "blue", lty=2)
dev.off()


# (4) chromosome 22 --- truncated
png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianAMF/ch22/hwe_vs_exact_Ch22_joint_truncated_t15.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome == 22]), -log10(snp$exact.pv.asian.jmf.ch22[snp$chromosome == 22]),
xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Joint Exact Test",
xlim=c(0,4), ylim=c(0,4),
main="Chr 22 - truncated, P-Values from GWASTools Exact Test vs. HardyWeinberg Joint Exact Test")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianAMF/ch22/hwe_vs_midp_Ch22_joint_truncated_t15.png",width=720, height=720)
plot(-log10(snp$hwe.asian.pval[snp$chromosome == 22]), -log10(snp$midp.pv.asian.jmf.ch22[snp$chromosome == 22]),
xlab="-log10(p-values) GWASTools Exact Test", ylab="-log10(p-values) HardyWeinberg Joint Mid-P Test",
xlim=c(0,4), ylim=c(0,4),
main="Chr 22 - truncated, P-Values from GWASTools Exact Test vs. HardyWeinberg Joint Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianAMF/ch22/exact_vs_midp_Ch22_joint_truncated_t15.png",width=720, height=720)
plot(-log10(snp$exact.pv.asian.jmf.ch22[snp$chromosome == 22]), -log10(snp$midp.pv.asian.jmf.ch22[snp$chromosome == 22]),
xlab="-log10(p-values) HardyWeinberg Joint Exact Test", ylab="-log10(p-values) HardyWeinberg Joint Mid-P Test",
xlim=c(0,4), ylim=c(0,4),
main="Chr 22 - truncated, P-Values from HardyWeinberg: Joint Exact Test vs. Joint Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()

png("/projects/wlhsu/plots/hwe/HardyWeinberg/asianAMF/ch22/midp_vs_midp_Ch22_joint_truncated_t15.png",width=720, height=720)
plot(-log10(snp$midp.pval.asian[snp$chromosome == 22]), -log10(snp$midp.pv.asian.jmf.ch22[snp$chromosome == 22]),
xlab="-log10(p-values) HardyWeinberg Mid-P Test", ylab="-log10(p-values) HardyWeinberg Joint Mid-P Test",
xlim=c(0,4), ylim=c(0,4),
main="Chr 22 - truncated, P-Values from HardyWeinberg: Mid-P Test vs. Joint Mid-P Test")
abline (a=0, b=1, col="red")
dev.off()





rm(list=objects())


