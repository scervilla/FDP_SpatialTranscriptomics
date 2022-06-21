library(dplyr)
library(survival)
library(survminer)
library(RTCGA.clinical)
library(RTCGA.mRNA)

#Load spatial transcriptomics objects
STobject1 <- readRDS("Downloads/OV4A.rds")
STobject2 <- readRDS("Downloads/5merged.rds")

#select top 100 genes for each sample
genes4 <- names(sort(STobject@reductions$NMF@feature.loadings[, 12], decreasing = T))[1:100]
genes5 <- names(sort(STobject@reductions$NMF@feature.loadings[, 14], decreasing = T))[1:100]
#intersect of the genes
genes <- intersect(genes4, genes5)
rm(STobject)

#extract clinical data from ovarian cancer
clin <- survivalTCGA(OV.clinical,
                     extract.cols="admin.disease_code")
#extract expression data from ovarian cancer
expr <- expressionsTCGA( OV.mRNA,
                         extract.cols = genes)

#merge data by patient barcode
expr$bcr_patient_barcode2 <- substr(expr$bcr_patient_barcode, 1, 12)
expr <- expr[!duplicated(expr$bcr_patient_barcode2),]
rownames(clin) <- clin$bcr_patient_barcode
x <- clin[intersect(clin$bcr_patient_barcode, expr$bcr_patient_barcode2),]
rownames(expr) <- expr$bcr_patient_barcode2
y<- expr[intersect(clin$bcr_patient_barcode, expr$bcr_patient_barcode2),]
a <- cbind(x, y)
x <- cbind(x, a[, colnames(y)[4:ncol(y)-1]])

#assign cut points of High and Low groups
cut <- surv_cutpoint(
  x,
  time = "times",
  event = "patient.vital_status",
  variables = colnames(y)[4:ncol(y)-1]
)
x <- surv_categorize(cut)
for (g in colnames(y)[4:ncol(y)-1]) {
  x[,g] <- factor(x[,g], levels = c("low", "high"))
}

#perform survival analysis for each predictor (gene)
#genes found in TCGA
covariates <-  colnames(y)[4:ncol(y)-1]
univ_formulas <- c()
for (g in covariates) {
  univ_formulas <- c(univ_formulas, as.formula(paste('Surv(times, patient.vital_status)~', g)))
}

univ_models <- list()
for (i in 1:18) {
  univ_models[[i]] <- coxph(univ_formulas[[i]], data = x)
}
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         res<-c(p.value, HR, HR.confint.lower, HR.confint.upper)
                         names(res)<-c("p.value", "HR", "CIlower", 
                                       "CIupper")
                         return(res)
                       })
#results
df <- as.data.frame(do.call(rbind, univ_results))
df$group <- factor(covariates, levels = rev(covariates))

#assign significance
df$Label <- NA
df$Label[df$p.value<0.001]<-'***'
df$Label[df$p.value<0.01 & is.na(df$Label)]<-'**'
df$Label[df$p.value<0.05 & is.na(df$Label)]<-'*'
df$Label[is.na(df$Label)] <-""

df$p.value <- format.pval(df$p.value)
df$p <- paste0(df$p.value, df$Label)

#plot hazard ratios
p <- ggplot(df, aes(y=group, x=HR)) +
  geom_point(shape=18, size=4, fill="white") +
  geom_errorbar(width=0.3, aes(xmin=CIlower, xmax=CIupper), colour="black") +
  annotate(geom="text", label=df$p, x = max(df$CIupper)+0.25, y=18:1, size=3) + 
  annotate(geom="text", label="P-value", x = max(df$CIupper)+0.25, y=19, size =3.5) +
  geom_vline(xintercept=1, linetype="dashed", color ="gray") + 
  labs(y="High expression group", x = "Hazard Ratio") + 
  coord_cartesian(clip = 'off') +
  theme_tufte()


pdf("Desktop/HR.pdf", width = 6, height = 7)
plot(p)
dev.off()


