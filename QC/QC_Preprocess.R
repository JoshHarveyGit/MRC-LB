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
# 5. Relatedness Check
# 6. Raw PCA Generation (Check Regional Variation between CNG and PFC)
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

setwd("/mnt/data1/LewyBodyMRC/Methylation/") #Set working directory


#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### 1. File Generation #########
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #

# Load in mSet Object
load()

# Load in RGSet Object
load()

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

write.csv(fullPheno,"/mnt/data1/LewyBodyMRC/Methylation/FullPheno.csv",row.names = F)



unique(pheno$Basename == colnames(FullMset))
unique(pheno$Basename == colnames(FullRGset))

#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 2. Median Intensity Check #####
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

# Median Methylated and Unmethylated intensities will be calculated and plotted across all plates to explore potential
# batch effects across the study. It will also be double checked and validated with previous M / U intensities to 
# validate there were no sample switches during data merging. 

## extract sample intensities 
m_intensities<-methylated(fullMset) ## this gives a matrix where each row is a probe and each column a sample
u_intensities<-unmethylated(ffullMset)

## summarise the intensities of each sample with a single value, the median 
M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)

##Correlate to previous QC


  
#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 3. BS-Conversion Check ########
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 4. Sex Check ##################
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 5. Relatedness Check ##########
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 6. PCA Generation #############
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 7. Dasen Normalisation ########
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 8. Cell Type Decon/Cetygo #####
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 9. Summary / Export ###########
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #
