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

setwd("/mnt/data1/LewyBodyMRC/Methylation/") #Set working directory


#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### 1. File Generation #########
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #

# Load in mSet Object
load("/mnt/data1/LewyBodyMRC/Methylation/QC/mSetFull.rdata")

# Load in pheno file
pheno <- read.csv("/mnt/data1/LewyBodyMRC/Pheno/methPheno.csv", header = T, stringsAsFactors = F)


unique(pheno$Basename == colnames(fullMset))
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
m_intensities <- methylated(fullMset) ## this gives a matrix where each row is a probe and each column a sample
u_intensities <- unmethylated(fullMset)

## summarise the intensities of each sample with a single value, the median 
M.median <- apply(m_intensities, 2, median)
U.median <- apply(u_intensities, 2, median)
M.U.ratio <- M.median / U.median

pheno <- cbind(pheno, M.median, U.median, M.U.ratio)

# Plot histograms of median intensity
pdf("QC/Plots/Intensity_Histograms.pdf")
hist(M.median, xlab = "Median M intensity")
abline(v = 2000, col = "red")
hist(U.median, xlab = "Median U intensity")
abline(v = 1000, col = "red")
dev.off()

# Scatter and boxplots of intensity by plate
pdf("QC/Plots/Intensity_PlateEffect.pdf", width = 12)
ggplot(data = pheno[-which(pheno$SampleID == "Meth_Control"),], aes(x = M.median, y = U.median, color = as.factor(plateNumber)))+
  geom_point()+
  theme_cowplot()+
  ggtitle("Median unMethylated intensity by Plate")+
  labs(color = "Plate")
ggplot(data = pheno[-which(pheno$SampleID == "Meth_Control"),], aes(x = as.factor(plateNumber), y = U.median))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_cowplot()+
  ggtitle("Median unMethylated intensity by Plate")
ggplot(data = pheno[-which(pheno$SampleID == "Meth_Control"),], aes(x = as.factor(plateNumber), y = M.median))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_cowplot()+
  ggtitle("Median Methylated intensity by Plate")
ggplot(data = pheno[-which(pheno$SampleID == "Meth_Control"),], aes(x = as.factor(plateNumber), y = M.U.ratio))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_cowplot()+
  ggtitle("Median Methylated/UnMethylated Ratio by Plate")
dev.off()

# Scatter and boxplots of intensity by plate / chip
pdf("QC/Plots/Intensity_PlateChipEffect.pdf", width = 12)
ggplot(data = pheno[-which(pheno$SampleID == "Meth_Control"),], aes(x = M.median, y = U.median, color = as.factor(chipNumber)))+
  geom_point()+
  theme_cowplot()+
  ggtitle("Median unMethylated intensity by Chip")+
  labs(color = "Chip")+
  facet_wrap(~plateNumber)
ggplot(data = pheno[-which(pheno$SampleID == "Meth_Control"),], aes(x = as.factor(chipNumber), y = U.median))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_cowplot()+
  ggtitle("Median unMethylated intensity by Chip")+
  facet_wrap(~plateNumber)+
  ylab("Median U")+
  xlab("Chip Number")
ggplot(data = pheno[-which(pheno$SampleID == "Meth_Control"),], aes(x = as.factor(chipNumber), y = M.median))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_cowplot()+
  ggtitle("Median Methylated intensity by Chip")+
  facet_wrap(~plateNumber)+
  xlab("Chip Number")+
  ylab("Median M")
ggplot(data = pheno[-which(pheno$SampleID == "Meth_Control"),], aes(x = as.factor(chipNumber), y = M.U.ratio))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_cowplot()+
  ggtitle("Median Methylated/UnMethylated Ratio by Chip")+
  facet_wrap(~plateNumber)+
  xlab("Chip Number")+
  ylab("Median M / U Ratio")
dev.off()

# Store Categorical and Continuous Confounders
catConf <- c("Brain.Bank","Gender","Clinical_Bin")
conConf <- c("Braak.LB.Numeric","Braak.Tangle", "Post.Mortem.Delay","Age")
varList <- c("M.median", "U.median","M.U.ratio")

#The below for loop makes boxplots of intensity metrics per Categorical featuress
for(cat in catConf){
  pdf(file = paste("QC/Plots/IntensityTest",cat,".pdf"), width = 12)
  for(var in varList){
  print(
    ggplot(data = pheno[-which(pheno$SampleID == "Meth_Control"),], aes_string(x = cat, y = var))+
      geom_boxplot()+
      geom_jitter(width = 0.2)+
      theme_cowplot()+
      ggtitle(paste(var,"by",cat))
  )}
  dev.off()
}

#The below for loop makes scatter of intensity metrics per continuous features
for(con in conConf){
  pdf(file = paste("QC/Plots/IntensityTest",con,".pdf"), width = 12)
  for(var in varList){
    print(
      ggplot(data = pheno[-which(pheno$SampleID == "Meth_Control"),], aes_string(x = con, y = var))+
        geom_point()+
        theme_cowplot()+
        ggtitle(paste(var,"by",con))
    )}
  dev.off()
}

#Clean up
rm(cat,catConf,con,conConf,M.median,M.U.ratio,U.median, var,m_intensities,u_intensities,varList)


#!!!!!!SUMMARY  
# Based on visual assessment of intesities, only two samples fail on intensity thresholding these being:
pheno[which(pheno[-which(pheno$SampleID == "Meth_Control"),"U.median"] < 1000),"SampleID"]
# "PD309_PFC" and "C085_CNG" 

# Looking at plate effects: Generally plates are comparable, Plates 1 & 2 have slightly elevated intensities

# The Fully methylated control (203991410002_R01C01) on plate 4 shows low median methylated intensity, although its
# ratio is still high. The Fully methylated control (203991410002_R01C01) on plate 1 has a slightly low ratio, maybe
# resulting from background.


#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 3. BS-Conversion Check ########
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

bs<-bscon(fullMset)

pdf("QC/Plots/BS_Conversion.pdf")
hist(bs, xlab = "Median % BS conversion", breaks = seq(0,100, by = 2))
abline(v = 90, col = "red")
dev.off()

pheno <- cbind(pheno, bs)


#!!!!! SUMMARY
# BS conversion rates are high, with the majority > 90. Only 5 samples show bs rate under
pheno[which(pheno$bs < 90),c("SampleID")]
          # Basename      SampleID  
# 203991460069_R01C01     PD904_CNG 
# 203991410033_R04C01  20163903_PFC 
# 203991410125_R08C01      C085_PFC 
# 203991410064_R07C01  20130400_PFC 
# 203991410064_R08C01 RI95/1030_PFC 


#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 4. Sex Check ##################
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

# Source QC gender clustering functions
source("/mnt/data1/ExeterEWASPipeline/R/findGenderPC.r")
source("/mnt/data1/ExeterEWASPipeline/R/clusterGender.r")

#Make Gender a binary sex variable
pheno[which(pheno$Gender == ""),"Gender"] <- NA
pheno$Sex <- as.numeric(as.factor(pheno$Gender))
pheno <- pheno[-which(pheno$SampleID == "Meth_Control"),]


#Save beta matrix
betas <- betas(fullMset)
betas <- betas[,which(colnames(betas) %in% pheno$Basename)]

pdf("QC/Plots/ClusterGenders.pdf")
predSex1 <-findGenderPC(betas, pheno$Sex)
predSex2 <-clusterGender(betas, pheno$Sex,thres = 0.8)
dev.off()

View(pheno[which(pheno$Sex != predSex1),])


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
