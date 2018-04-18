source("https://bioconductor.org/biocLite.R")
biocLite("bsseq") 
biocLite("DSS")

library(bsseq)
library(DSS)

#####Load in Bismark Data####

getwd()
# setwd("/Users/laurenschmitz/Desktop/03-methylation_call/")

bsseq_combined <- read.bismark(files = c("T101d_Ctrl_M_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T102c_Lead_F_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T103c_Lead_F_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T103d_Lead_M_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T105d_Lead_M_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T106c_Lead_F_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T106d_Lead_M_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T107d_Lead_M_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T108d_Ctrl_M_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T109c_Ctrl_F_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T109d_Ctrl_M_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T112c_Lead_F_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T112d_Lead_M_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T113c_Ctrl_F_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T116c_Lead_F_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T116d_Lead_M_liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T118c_Ctrl_F_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T118d_Ctrl_M_liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T118h_Ctrl_M_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T119e_Ctrl_F_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "T119g_Ctrl_F_Liver_5mo_R1_trimmed.fq.gz_bismark_bt2.bismark.cov"),
                               sampleNames = c("101d_CM","102c_TF","103c_TF","103d_TM",
                                               "105d_TM", "106c_TF","106d_TM","107d_TM","108d_CM","109c_CF",
                                               "109d_CM","112c_TF","112d_TM","113c_CF","116c_TF", "116d_TM",
                                               "118c_CF","118d_CM","118h_CM","119e_CF","119g_CF"),
                               rmZeroCov = TRUE,
                               strandCollapse = FALSE,
                               fileType = "cov",
                               verbose = TRUE)

#### Process data ####

# populate pData object
pData(bsseq_combined) <- data.frame(pData(bsseq_combined), type = exposure, sex = sex)

# get average coverage across all samples
# based on results below, we will have to remove sample 112c_TF from analysis as it has only 4x coverage
round(colMeans(getCoverage(bsseq_combined)), 1)

length(bsseq_combined)

# number of CpGs covered by at least one read across all samples
sum(rowSums(getCoverage(bsseq_combined) >= 1) == 21)

# remove CpGs with little or no coverage, that is
# CpG site in atleast 5 samples per group with more than 10 reads
# DICUSS: What is the proper threshold?

bsseq_cov <- getCoverage(bsseq_combined)
keep_loci_bsseq <- which(rowSums(bsseq_cov[, bsseq_combined$type == "Control"] >= 10) >= 5 &
                         rowSums(bsseq_cov[, bsseq_combined$type == "Lead"] >= 10) >= 5)
length(keep_loci_bsseq) 

# keep only those sites
bsseq_filtered <- bsseq_combined[keep_loci_bsseq, ]


#####Design Covariate Matrix####

# Must establish design prior to running multifactor function.
# all samples combined in single matrix 
n <- c(1:21)

exposure <- c("Control","Lead", "Lead","Lead", "Lead","Lead","Lead", "Lead","Control", "Control",
              "Control", "Lead","Lead", "Control","Lead", "Lead","Control", "Control","Control", 
              "Control", "Control")

sex <- c("M", "F", "F", "M","M", "F", "M", "M", "M", "F",
         "M", "F", "M", "F", "F", "M","F", "M", "M", "F", "F")

design_combined <- data.frame(n,exposure,sex)


####CALCULATE PERCENT DIFFERENCE METHYLATION IN MALES AND FEMALES####

DMLfit_m = DMLtest(bsseq_filtered, group1=c("103d_TM","105d_TM", "106d_TM","107d_TM","112d_TM","116d_TM"),
                                   group2=c("101d_CM","108d_CM","109d_CM","118d_CM","118h_CM"))

DMLfit_f = DMLtest(bsseq_filtered, group1=c("102c_TF","103c_TF","106c_TF","112c_TF","116c_TF"),
                                   group2=c("109c_CF","113c_CF","118c_CF","119e_CF","119g_CF"))

DMLfit_fold = DMLtest(bsseq_combined, group1=c("102c_TF","103c_TF","106c_TF","112c_TF","116c_TF"),
                                   group2=c("109c_CF","113c_CF","118c_CF","119e_CF","119g_CF"))

# volcano-plots 
par(mfrow = c(1,2))
plot(DMLfit_m$diff*100, -log10(DMLfit_m$pval), pch = 20, cex = .5, col = ifelse(DMLfit_m$fdr < 0.05, 'red', 'grey'),
     xlab = "%-diff methylation (Trt-Crt), Male", ylab = "-log10(p-val)")
abline(h = -log10(0.05)) # nominal-pvalue cutoff line

plot(DMLfit_f$diff*100, -log10(DMLfit_f$pval), pch = 20, cex = .5, col = ifelse(DMLfit_f$fdr < 0.05, 'red', 'grey'),
     xlab = "%-diff methylation (Trt-Crt), Female", ylab = "-log10(p-val)")
abline(h = -log10(0.05)) # nominal-pvalue cutoff line

# further stats
head(DMLfit_m)
dim(DMLfit_m) ##1,484,830 total CpG sites tested in males

head(DMLfit_f)
dim(DMLfit_f) ##1,484,291 total CpG sites tested in females

# save files
write.csv(DMLfit_m, file = "DMLtest_Male.csv")
write.csv(DMLfit_f, file = "DMLtest_Female.csv")

# subset data frames based on FDR cutoff of 0.05
DMLfit_m_05 <- subset(DMLfit_m, fdr < 0.05)
dim(DMLfit_m_05) ##1215 CpG sites

DMLfit_f_05 <- subset(DMLfit_f, fdr < 0.05)
dim(DMLfit_f_05) ##994 sites

## Histogram of p-values and test-statistics
# Males
par(mfrow=c(1,2))
hist(DMLfit_f$stat, 50, main="Test statistics, Males", xlab="")
hist(DMLfit_f$pval, 50, main="P-values, Males", xlab="")

#Females
par(mfrow=c(1,2))
hist(DMLfit_fold$stat, 50, main="Test statistics, Females", xlab="")
hist(DMLfit_fold$pval, 50, main="P-values, Females", xlab="")

####CALCULATE DMRs####
#p-value cut-off of 0.001 (~= to FDR=0.05 in males and females)

#Males
dmr_m_001 <- callDMR(DMLfit_m, p.threshold=0.001)
dim(dmr_m_001) ##622 DMRs

#Females
dmr_f_001 <- callDMR(DMLfit_f, p.threshold=0.001)
dim(dmr_f_001) ##307 DMRs

write.csv(dmr_m_001, file = "DMR_Male_001.csv")
write.csv(dmr_f_001, file = "DMR_Female_001.csv")

##p-value cut-off of 0.0001 (~= to FDR=0.01 in males and females)
dmr_m_0001 <- callDMR(DMLfit_m, p.threshold=0.0001)
dim(dmr_m_0001) ##538 DMRs

#Females
dmr_f_0001 <- callDMR(DMLfit_f, p.threshold=0.0001)
dim(dmr_f_0001) ##222 DMRs

write.csv(dmr_m_0001, file = "DMR_Male_0001.csv")
write.csv(dmr_f_0001, file = "DMR_Female_0001.csv")

####DSS Package - multifactor modeling ####

# Multifactorial models for examining sex and exposure effects in mice while 
# controlling for mouse and sex (NO INTERACTION):
DMLfit_combined = DMLfit.multiFactor(bsseq_filtered, design=design_combined,  formula=~exposure+sex)
head(DMLfit_combined)

#FIND SIGNIFICANT DMLs FOR EFFECT OF EXPOSURE
#Note: Results are data frames with chromosome number, CpG site position, test statistics, p-values (from
#normal distribution), and FDR.
DMLtest_combined.exp = DMLtest.multiFactor(DMLfit_combined, coef=2)
head(DMLtest_combined.exp)
dim(DMLtest_combined.exp) ##7,806,583 sites tested total


#Sort by exposure term
#sort the data frames (and remove NAs) using following codes:
DMLtest_combined.exp_sort1 <- DMLtest_combined.exp[order(DMLtest_combined.exp$pvals, na.last=NA), ]
head(DMLtest_combined.exp_sort1)
#Subset data frame without NAs
sum(is.na(DMLtest_combined.exp$stat)) ##5,553,056 sites are could not be tested (NAs)

#Subset data frame based on FDR cutoff of 0.05
DMLtest_combined.exp_sort2 <- subset(DMLtest_combined.exp_sort1, fdrs < 0.05)
head(DMLtest_combined.exp_sort2) 
dim(DMLtest_combined.exp_sort1) ##2253527 CpG sites tested total
dim(DMLtest_combined.exp_sort2) ##61225 CpG sites in data frame

#Write CSVs for exposure term p-values from models for exposure 
write.csv(DMLtest_combined.exp_sort2, file = "DSS_DMLtest_Exp.csv")

#Histogram of p-values and test-statistics
par(mfrow=c(1,2))
hist(DMLtest_combined.exp$stat, 50, main="test statistics", xlab="")
hist(DMLtest_combined.exp$pvals, 50, main="P values", xlab="")

####FIND SIGNIFICANT DMLs FOR EFFECT OF SEX####
DMLtest_combined.sex = DMLtest.multiFactor(DMLfit_combined, coef=3)
head(DMLtest_combined.sex)

DMLtest_combined.sex_sort1 <- DMLtest_combined.sex[order(DMLtest_combined.sex$pvals, na.last=NA), ]
head(DMLtest_combined.sex_sort1)

#Subset data frame based on FDR cutoff of 0.05
DMLtest_combined.sex_sort2 <- subset(DMLtest_combined.sex_sort1, fdrs < 0.05)
head(DMLtest_combined.sex_sort2) 
dim(DMLtest_combined.sex_sort1) ##2253527 CpG sites tested total
dim(DMLtest_combined.sex_sort2) ##89455 CpG sites in data frame

#Write CSVs for exposure term p-values from models for sex
write.csv(DMLtest_combined.sex_sort2, file = "DSS_DMLtest_Sex.csv")

#Histogram of p-values and test-statistics
par(mfrow=c(1,2))
hist(DMLtest_combined.sex$stat, 50, main="test statistics", xlab="")
hist(DMLtest_combined.sex$pvals, 50, main="P values", xlab="")

####DMR ANALYSIS--MODEL WITHOUT INTERACTION TERM####

#DMR ANALYSIS FOR EXPOSURE
#P-Value threshold=0.05
dmrscombined <- callDMR(DMLtest_combined.exp, p.threshold=0.05)
dim(dmrscombined) ##2,419 DMRs found

#DMR ANALYSIS for SEX

#P-Value threshold=0.05
dmr_sex_05 <- callDMR(DMLtest_combined.sex, p.threshold=0.05)
dim(dmr_sex_05) ##3,785 DMRs found
write.csv(dmr_sex_05, file = "DMRtest_Sex_pval_05.csv")

####EXTRA CODE####
DML_all <- DMLtest(bsseq_filtered, c("101d_CM", "108d_CM","109c_CF","109d_CM","113c_CF",
                                             "118c_CF","118d_CM","118h_CM","119e_CF","119g_CF"), 
                           c("102c_TF","103c_TF","103d_TM","105d_TM","106c_TF",
                             "106d_TM","107d_TM","112c_TF","112d_TM","116c_TF","116d_TM"),
                           equal.disp = FALSE, smoothing = FALSE)

par(mfrow=c(1,1))
plot(DML_all$diff*100, -log10(DML_all$pval), pch = 20, cex = .5, col = ifelse(DML_all$fdr < 0.05, 'red', 'grey'),
     xlab = "%-diff methylation (Trt-Crt)", ylab = "-log10(p-val)")
