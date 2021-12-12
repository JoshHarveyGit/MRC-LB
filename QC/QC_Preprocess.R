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
