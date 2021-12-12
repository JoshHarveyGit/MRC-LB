#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #  
######### Introduction ################
#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   # 

# Aim: Quality Control EPIC Array Data generated for the MRC Lewy Body Project. 
# Run By: Josh Harvey / Jenny Imm
# Date: 13/12/21


# 1. Load Mset and RGSet Objects
# 2. Median Intensity Check
# 3. BS-Conversion Check
# 4. Sex Check
# 5. Raw PCA Generation (Check Regional Variation between CNG and PFC)
# 6. Relatedness Check
# 7. Dasen Normalisation
# 8. Cell Type Deconvolution and CETYGO quantification
# 9. Final Summary, generate analysis files


#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### Packages ####################
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #

library(minfi)
library(wateRmelon)
library(dplyr)
library(stringr)
library(corrplot)
library(cowplot)
library(rgl)

setwd("/mnt/data1/Josh/DLB_MRC/") #Set working directory


#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### 1. File Generation #########
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #

# MSet, RGSet and preliminary QC and Pheno files had already been generated during array processing by Jenny Imm 
# ("/mnt/data1/Jenny/EPIC_QC/LewyBodyMRC/"). These files will be loaded in and merged below.

# Generate a vector with locations of pheno

pheno_link <- c("/mnt/data1/Jenny/EPIC_QC/LewyBodyMRC/Plate1112/Pheno_1112.csv",
               "/mnt/data1/Jenny/EPIC_QC/LewyBodyMRC/Plate12/Pheno_12.csv",
               "/mnt/data1/Jenny/EPIC_QC/LewyBodyMRC/Plate34/Pheno_34.csv",
               "/mnt/data1/Jenny/EPIC_QC/LewyBodyMRC/Plate56/Pheno_56.csv",
               "/mnt/data1/Jenny/EPIC_QC/LewyBodyMRC/Plate78/Pheno_78.csv",
               "/mnt/data1/Jenny/EPIC_QC/LewyBodyMRC/Pheno_910.csv")

# The below for loop takes msets from the links above and aggregates them into a combined mset object "fullPheno"
for(x in pheno_link){
  y <- which(pheno_link == x)
  if(y == 1){
    fullPheno <- read.csv(file = x, header = T, stringsAsFactors = F)
    print(x)
  }else if(y == 6){
    tempPh <-read.csv(file = x, header = T, stringsAsFactors = F)[,-3]
    colnames(tempPh) <- colnames(fullPheno)
    tempPh$Chip <- str_sub(tempPh$Basename,1,12)
    fullPheno <- rbind(fullPheno, tempPh)
    rm(tempPh)
    print(x)
  }else{
    tempPh <-read.csv(file = x, header = T, stringsAsFactors = F)
    fullPheno <- rbind(fullPheno, tempPh)
    rm(tempPh)
    print(x)
  }
  }

# write.csv(fullPheno,"/mnt/data1/Josh/DLB_MRC/Pheno/FullPheno.csv",row.names = F)

# 
# fullPheno[which(fullPheno$plateNumber %in% c(5,6)),"Basename"]
# 
# msetEPIC <-  readEPIC(idatPath="/mnt/data1/Jenny/EPIC_QC/LewyBodyMRC/idats/", 
#                       barcodes=fullPheno[which(fullPheno$plateNumber %in% c(5,6)),"Basename"], 
#                       parallel = FALSE, force=T)
# 
# save(msetEPIC, file = "/mnt/data1/Josh/DLB_MRC/Meth/LBMRC56_msetEPIC_2.rdat")

"/mnt/data1/Josh/DLB_MRC/Meth/LBMRC56_msetEPIC_2.rdat"

# The below section merges and loads in msets

mset_link <- c("/mnt/data1/Jenny/EPIC_QC/LewyBodyMRC/Plate1112/LBMRC1112_msetEPIC.rdat",
                "/mnt/data1/Jenny/EPIC_QC/LewyBodyMRC/Plate12/LBMRC12_msetEPIC.rdat",
                "/mnt/data1/Jenny/EPIC_QC/LewyBodyMRC/Plate34/LBMRC34_msetEPIC.rdat",
               "/mnt/data1/Josh/DLB_MRC/Meth/LBMRC56_msetEPIC_2.rdat",
                "/mnt/data1/Jenny/EPIC_QC/LewyBodyMRC/Plate78/LewyBodyMRCLBMRC78_msetEPIC.rdat",
                "/mnt/data1/Jenny/EPIC_QC/LBMRC910_msetEPIC.rdat")


for(x in mset_link){
  y <- which(mset_link == x)
  if(y == 1){
    load(x)
    fullMset <- msetEPIC
    rm(msetEPIC)
    print(x)
  }else{
    load(x)
    tempMset <- msetEPIC
    fullMset <- combo(fullMset,tempMset)
    print(x)
    rm(tempMset)
    rm(msetEPIC)
  }
}  

colnames(fullMset) == fullPheno$Basename

save(fullMset, file = "Meth/FullMset.rdat")
