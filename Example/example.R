#Fit RBMR using CAD-CAD study
library(MR.LDP)
library(RBMR)
filescreen<-'heart attack myocardial infarction.txt'
fileexposure<-'exposure c4d.txt'
fileoutcome<-'outcome cardiogram.txt'
block_file<-'fourier_ls-all.bed'
stringname3<-'all_chr_1000G'
pva_cutoff<-1e-4;
scrres <-matchscreen(filescreen, fileexposure, fileoutcome,stringname3,pva_cutoff)
bh1 <- as.numeric(scrres$bh1)
bh2 <- as.numeric(scrres$bh2)
s12 <- as.numeric(scrres$s12)
s22 <- as.numeric(scrres$s22)
chr <- as.numeric(scrres$chr)
bp <- scrres$bp
rsname <- scrres$rsname
avbIndex <- scrres$idxin
idx4panel<-scrres$idx4panel
QCresult <- summaryQC(mhcstart, mhcend, bh1, bh2, s12, s22, bp, chr, rsname, avbIndex,idx4panel, Inf, Inf)
bh1new <- QCresult$bh1new
bh2new <- QCresult$bh2new
s12new <- QCresult$s12new
s22new <- QCresult$s22new
bpnew <- QCresult$bpnew
chrnew <- QCresult$chrnew
avbIndexnew <- QCresult$avbIndexnew
idx4panelnew<-QCresult$idx4panelnew
rsnamenew <- QCresult$rsnamenew
pmhc <- QCresult$pmhc
px <- QCresult$px
py <- QCresult$py
p <- length(avbIndexnew);
coreNum <- 1
lam <- 0.15
Rblockres<-Cal_blockR(bpnew,chrnew,avbIndexnew-1,idx4panelnew,block_file, stringname3,coreNum,lam)
F4Rblock<-Rblockres$F4Rblock
block_inf<-Rblockres$block_inf
nblocks<-Rblockres$nblocks
bh1<-bh1new
bh2<-bh2new
se1<-s12new
se2<-s22new
mu<-rep(0.01,p)
muA<-rep(0.01,p)
sgga2<-0.01
sgal2<-0.01
beta0<-0
IterMax<-10000
epsStopLogLik<-1e-7
alphag<-8
betag<-4
RBMR_Ha<- RBMR_func_block(F4Rblock, block_inf, nblocks, bh1, bh2, se1, se2,
                          gamma, alpha, sgga2, sgal2, beta0, constr=0, epsStopLogLik, IterMax,alphag,betag)

RBMR_H0<- RBMR_func_block(F4Rblock, block_inf, nblocks, bh1, bh2, se1, se2,
                          gamma, alpha, sgga2, sgal2, beta0, constr=1, epsStopLogLik, IterMax,alphag,betag)
#calculate the p-value
tstat <- 2*(RBMR_Ha$tstat-RBMR_H0$tstat)
pvalue<-pchisq(tstat,1,lower.tail=F)
#causal effect
beta_hat <- RBMR_Ha$beta0
#standard error of causal effect
se_hat <- abs(beta_hat/sqrt(tstat))
