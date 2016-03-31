"/protected/projects/pulmarray/Biollegro/RDS/Percepta/"
cd /protected/projects/pulmarray/Biollegro/RDS/Percepta/
module load R/R-3.1.1
R
.libPaths("/unprotected/projects/cbmhive/R_packages/R-3.0.0")



##Load modules 
library(Biobase)
library("GEOquery")
library("Biobase")
library("prodlim")
library("limma")
library(NMF)
library(Heatplus)
library(heatmap3)
library(GSVA)



allegro = function(outputName, FDR =.05, inputModel = 'simple', ;
smkStatus = 'all', varFilter = 20000){


setwd("/protected/projects/pulmarray/Biollegro/RDS/Percepta/")
data = readRDS("bronch_all_affy_cdf.rds")

###Identify extreme low and high data (do this before stratifying by smk status)
extremeLowPY = quantile(pData(data)$PY, .05, na.rm = TRUE)
extremeHighPY = quantile(pData(data)$PY, .95, na.rm = TRUE)
extremeLowIndices = which(pData(data)$PY <= extremeLowPY)
extremeHighIndices = which(pData(data)$PY >= extremeHighPY)

###Filter out extreme low and high data
data= data[,-which(is.na(pData(data)$PY))]
data = data[,-union(extremeLowIndices, extremeHighIndices)]

###Restrict to currents or former smokers as desired. 
if (smkStatus == 'currents'){data = data[,which(pData(data)$SMK == 1)]
} else if(smkStatus == 'formers'){data = data[,which(pData(data)$SMK == 2)]}

###Restrict gene space to 20k most varable genes
MAD<-apply(exprs(data), 1, mad)
topMad <- order(MAD, decreasing = T)[1:varFilter]
subset = data[topMad,]



###Define model object and contrasts to feed limma
model = model.matrix(~pData(data)$PY )
colnames(model) = c("Intercept", "PY")

###Feed limma model and expressionSet to identify genes that correlate with PY at q < .05
contrast.matrix <- makeContrasts(PY, levels = model)
fit = lmFit(exprs(data), model)
fit = contrasts.fit(fit, contrast.matrix)
fit = eBayes(fit)
limmaRes = topTable(fit, adjust.method = "BH", n = Inf, sort.by = "P")
limmaRes.sig = subset(limmaRes, adj.P.Val < FDR)

###Get gene names
geneNames = fData(data)[which(  rownames(fData(data)) %in% rownames(limmaRes.sig)),]

###Save as RDS
setwd("~")
saveRDS(geneNames, outputName)

