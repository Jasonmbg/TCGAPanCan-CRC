############# MSI Inference Status for ACCC-CRC project concordance #####################

library(PreMSIm)
library(tidyverse)
library(here)

dr_here()

# this
pt <- here("Data/TCGAPanCan.CRC.Exp.Norm.Filt.RNASeq.tsv")
# or equivalent to:
pt <- here("Data", "TCGAPanCan.CRC.Exp.Norm.Filt.RNASeq.tsv")

input_crc_tcga_pancan = data_pre(pt, type = "Symbol")

Out.MSI.PanCanTCGA.CRC <- msi_pre(input_crc_tcga_pancan) # Prediction results  ("1" indicates MSI-high, and "0" MSI-low/microsatellite stability).

# write_tsv(Out.MSI.PanCanTCGA.CRC, file="MSIStatus.PreMSIm.Predict.PanCanTCGA.CRC.NormFilt.tsv") # 450 samples

# evaluate these that have either values, not the ones that were identified as intermediate by MSI MANTIS;

## also load the respective clinical information to fetch also the MSI MANTIS estimates

final.updated.clin.mut.dat <- read_tsv(here("Data",
"PanCanTCGA.CRC.Clinical.MSIReform.MANTIS.MutStatus.450.tsv"))

dim(final.updated.clin.mut.dat)
dim(Out.MSI.PanCanTCGA.CRC)

combo.MSI.eval.dat <- inner_join(final.updated.clin.mut.dat, Out.MSI.PanCanTCGA.CRC, 
by=c("Tumor_Sample_Barcode"="Sample")) %>% dplyr::select(Tumor_Sample_Barcode, MSI_status_MANTIS, 
MSI_status)

# write_tsv(combo.MSI.eval.dat, file="PanCan.TCGA.CRC.MSI.MANTIS.PreMSIm.Eval.450.tsv")

################# continue with extracting the prediction metrics #########################

tcga.pred.sel <- combo.MSI.eval.dat %>% filter(!MSI_status_MANTIS=="Indeterminate")

premsim <- as.character(tcga.pred.sel$MSI_status)
premsim.refined <- recode(premsim, "0"="MSS", "1"="MSI")
tcga.pred.sel$MSI_status <- premsim.refined

cm <- table(tcga.pred.sel$MSI_status, tcga.pred.sel$MSI_status_MANTIS)

accuracy <- sum(cm[1], cm[4]) / sum(cm[1:4]) 
precision <- cm[4] / sum(cm[4], cm[2])
sensitivity <- cm[4] / sum(cm[4], cm[3]) 
fscore <- (2 * (sensitivity * precision))/(sensitivity + precision)
specificity <- cm[1] / sum(cm[1], cm[2]) 

###########################################################################################
