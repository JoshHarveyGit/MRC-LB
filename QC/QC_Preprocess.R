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
# 7. p_filter
# 8. Dasen Normalisation
# 9. Cell Type Deconvolution and CETYGO quantification
# 10. Final Summary, generate analysis files


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
library(FlowSorted.DLPFC.450k)
library(gdsfmt)

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
# unique(pheno$Basename == colnames(FullRGset))


pheno <- pheno[match(pheno$Basename, colnames(fullMset)),]

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
abline(v = 1500, col = "red")
hist(U.median, xlab = "Median U intensity")
abline(v = 1000, col = "red")
dev.off()

# Scatter and boxplots of intensity by plate
pdf("QC/Plots/Intensity_PlateEffect.pdf", width = 12)
ggplot(data = pheno[-which(pheno$SampleID == "Meth_Control"),], aes(x = M.median, y = U.median, color = as.factor(plateNumber)))+
  geom_point()+
  theme_cowplot()+
  ggtitle("Median intensity by Plate")+
  labs(color = "Plate")+
  geom_hline(aes(yintercept =1000), color = "red",linetype='dotted')+
  geom_vline(aes(xintercept =1500), color = "red",linetype='dotted')
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
unique(pheno[which(pheno$M.median < 1500 | pheno$U.median < 1000),"SampleID"])
# "20163903_PFC" "Meth_Control" "C085_PFC"

# 1500 for Methylated
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
######### 4. Outlier Check ##############
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

# Remove the fully methylated controls
controls <- pheno[which(pheno$SampleID == "Meth_Control"),]
write.csv(controls, "QC/MethControlSummary.csv", row.names = F)

pheno <- pheno[-which(pheno$SampleID == "Meth_Control"),]

cMset <- fullMset[,which(colnames(fullMset) %in% controls$Basename)]
save(cMset,file = "QC/ControlMset.R")

fullMset <- fullMset[,which(colnames(fullMset) %in% pheno$Basename)]

pdf("QC/Plots/Outliers.pdf")
outliers <- outlyx(fullMset,iqr=TRUE, iqrP=2, pc=1,
                   mv=TRUE, mvP=0.15, plot=TRUE)
dev.off()

# Outliers flags up 6 cases, 3 of which already failed either due to intensity or bs-conversion

pheno$outliers <- outliers$outliers

write.csv(pheno[which(pheno$outliers == TRUE),], file = "QC/SummaryTables/OutlierSamples.csv",row.names = F)



#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 5. Sex Check ##################
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

# Source QC gender clustering functions
source("/mnt/data1/ExeterEWASPipeline/R/findGenderPC.r")
source("/mnt/data1/ExeterEWASPipeline/R/clusterGender.r")

#Make Gender a binary sex variable
pheno[which(pheno$Gender == ""),"Gender"] <- NA
pheno$Sex <- as.numeric(as.factor(pheno$Gender))
pca <- prcomp(betas)

#Save beta matrix
betas <- betas(fullMset)
betas <- betas[,which(colnames(betas) %in% pheno$Basename)]

pdf("QC/Plots/ClusterGenders.pdf")
predSex1 <-findGenderPC(betas, pheno$Sex)
predSex2 <-clusterGender(betas, pheno$Sex,thres = 0.8)
dev.off()

table(pheno$Sex,predSex1)

pheno$predSex1 <- predSex1

write.csv(pheno[which(pheno$Sex != predSex1),],file = "QC/SummaryTables/SexCheckFails.csv", row.names = F)


#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 6. Relatedness Check ##########
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

# Extract SNP probes
betas.rs<-betas[grep("rs", rownames(betas)),]


#Correlate probes 
snpCor<-cor(betas.rs, use="pairwise.complete.obs")

# Annotate number of expected relations
pheno$n <- NULL
counts <- pheno %>% group_by(ID) %>% tally()
pheno <- left_join(pheno, counts,by = "ID")


##Build an expected actual relationship matrix, annotating actually related samples as 1 and unrelated as 0
actual_relations <- data.frame()
for(i in 1:nrow(pheno)){
  x <- as.integer(pheno$ID[i] == pheno$ID)
  actual_relations <- rbind(actual_relations,x)
}
actual_relations <- as.matrix(actual_relations)
colnames(actual_relations) <- pheno$Basename
rownames(actual_relations) <- pheno$Basename

# Plot actual relationship matrix
# ActMatrix <- corrplot(actual_relations, method = "color", title = "Actual Relationship Matrix",tl.pos = "n")

# Build an observed inferred relationship matrix from the SNP correlations, coded as a boolean 0 for < 0.95 and 1 for > 0.95

inferred_relations <- snpCor
inferred_relations[] <- 0
inferred_relations[which(snpCor > 0.8)] <- 1

#By subtracting the inferred from the actual relationship matrix we can make a similarity matrix where 
#  0 = correct annotation of relationship
#  1 = pheno Related samples / SNP unrelated (cryptUnrelated)
# -1 = pheno unrelated samples / SNP related (cryptRelated)
sim_matrix <- actual_relations - inferred_relations

# Plot discordant samples
# SimMatrix <- corrplot(sim_matrix,method = "color", title = "Difference Relationship Matrix",tl.pos = "n")

# The below for loop takes the similarity matrix and annotates the number of cryptic related and cryptic unrelated samples per individual sample
# If samples show cryptic relatedness it also annotates the inferred related samples
n_cryptRel <- c()
n_cryptUnrel <- c()
for(x in 1:nrow(sim_matrix)){
  n_cryptRel <- append(n_cryptRel,length(which(sim_matrix[x,] < 0)))
  n_cryptUnrel <- append(n_cryptUnrel,length(which(sim_matrix[x,] > 0)))
  if(any(sim_matrix[x,] < 0)){
    pheno[x,"cryptID"] <- unique(pheno[which(sim_matrix[x,] < 0),"ID"])
  }
}

pheno <- cbind(pheno, n_cryptRel,n_cryptUnrel)


cryptics <- pheno[which(pheno$n_cryptRel > 0 | pheno$n_cryptUnrel > 0),]
write.csv(cryptics,"QC/SummaryTables/CrypticRelatedness.csv", row.names = F)


# 17 Samples showed either cryptic relatedness or a lack of genetic correlation between biological replicates. Assessment of relatedness highlights the following
# problem cases to exclude/
# "20130400_PFC"
# "20163903_PFC"
# "C026_CNG"
# "C085_PFC"
# "C090_PFC"
# "PD596_CNG"
# "RI95/1030_PFC"


CrypticBasenames <- c("203991410064_R07C01",
"203991410033_R04C01",
"203968030080_R03C01",
"203991410125_R08C01",
"203991460068_R03C01",
"203991410139_R04C01",
"203991410064_R08C01")

pheno$crypticQC <- FALSE
pheno[which(pheno$Basename %in% CrypticBasenames),"crypticQC"] <- TRUE

#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 7. PCA Generation #############
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

# Run and extract raw principle components of non-outlier samples

phenoB <- pheno[-which(pheno$outliers == TRUE),]

mSetB <- fullMset[,which(colnames(fullMset) %in% phenoB$Basename)]
betas <- betas(mSetB)

betas.com<-betas[complete.cases(betas),]
pca<-prcomp(betas.com)
# save(pca,file = "/mnt/data1/LewyBodyMRC/Methylation/QC/RawPCs.rdata")

plot(pca)

factorPheno <- cbind(phenoB, pca$rotation[,c(1:10)])

factorPheno$Brain_Region <- as.factor(factorPheno$Brain_Region)
factorPheno$Brain.Bank <- as.factor(factorPheno$Brain.Bank)
factorPheno$Braak.LB.Numeric <- as.factor(factorPheno$Braak.LB.Numeric)
factorPheno$Braak.Tangle <- as.factor(factorPheno$Braak.Tangle)

factorPheno$Brain_RegionN <- as.numeric(factorPheno$Brain_Region)
factorPheno$Brain.BankN <- as.numeric(factorPheno$Brain.Bank)
factorPheno$Braak.LB.NumericN <- as.numeric(factorPheno$Braak.LB.Numeric)
factorPheno$Braak.TangleN <- as.numeric(factorPheno$Braak.Tangle)


corPCA <- cor(factorPheno[,c("Sex","Brain.BankN","plateNumber","Brain_RegionN","Braak.LB.NumericN",
                             "Braak.TangleN","PC1","PC2","PC3","PC4","PC5","PC6")],use = "complete.obs")

pdf("/mnt/data1/LewyBodyMRC/Methylation/QC/Plots/Raw_PCA_Corrplot.pdf")
corrplot(corPCA)
dev.off()

pdf("/mnt/data1/LewyBodyMRC/Methylation/QC/Plots/RegionPC.pdf")
ggplot(factorPheno, aes(x = PC1, y = PC2, color = as.factor(Brain_Region)))+
  geom_point()+
  theme_cowplot()
dev.off()


#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 8. P-filter ###################
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #
mSet.pf <- pfilter(fullMset)

# 4 samples having 1 % of sites with a detection p-value greater than 0.05 were removed 
# Samples removed: 203991410033_R04C01 203991410125_R08C01 203991410064_R07C01 203991410064_R08C01 
# 1310 sites were removed as beadcount <3 in 5 % of samples 
# 3353 sites having 1 % of samples with a detection p-value greater than 0.05 were removed 

PF_Fail <- pheno[-which(pheno$Basename %in% colnames(mSet.pf)),"Basename"]

write.csv(pheno[which(pheno$Basename %in% PF_Fail),],"/mnt/data1/LewyBodyMRC/Methylation/QC/SummaryTables/PfilterFail.csv", row.names = F)

pheno$P.filter <- pheno$Basename %in% PF_Fail

# write.csv(pheno, "/mnt/data1/LewyBodyMRC/Methylation/QC/SummaryTables/FullQCMetrics.csv", row.names = F)

mSet <- fullMset[rownames(betas(fullMset)) %in% rownames(betas(mSet.pf)),]

#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 9. Dasen Normalisation ########
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #


FailedSamples <- pheno[which(pheno$M.median < 1500 | pheno$U.median < 1000 | 
                         pheno$bs < 80 | pheno$outliers == TRUE | pheno$crypticQC == TRUE | pheno$P.filter == TRUE),"Basename"]

#Subset passed samples
PassPheno <- pheno[-which(pheno$Basename %in% FailedSamples),]


mSet <- mSet[, which(colnames(mSet) %in% PassPheno$Basename)]

# Subset each region
CNGPheno <- PassPheno[which(PassPheno$Brain_Region == "CNG"),]
PFCPheno <- PassPheno[which(PassPheno$Brain_Region == "PFC"),]


mSetCNG <-mSet[,CNGPheno$Basename]
mSetPFC <-mSet[,PFCPheno$Basename]


#then use dasen to normalise
mSetCNG <-dasen(mSetCNG)
mSetPFC <-dasen(mSetPFC)
mSet <- dasen(mSet)

save(mSetCNG, file = "/mnt/data1/LewyBodyMRC/Methylation/QC/FinalData/Dasen_mSetCNG.rdata")
save(mSetPFC, file = "/mnt/data1/LewyBodyMRC/Methylation/QC/FinalData/Dasen_mSetPFC.rdata")
save(mSet, file = "/mnt/data1/LewyBodyMRC/Methylation/QC/FinalData/Dasen_mSetFull.rdata")


#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 10. Normalisation Check #######
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #
# Load in raw mSet Object
load("/mnt/data1/LewyBodyMRC/Methylation/QC/mSetFull.rdata")

# Load in normalised mSet Objects
load(file = "/mnt/data1/LewyBodyMRC/Methylation/QC/FinalData/Dasen_mSetCNG.rdata")
load(file = "/mnt/data1/LewyBodyMRC/Methylation/QC/FinalData/Dasen_mSetPFC.rdata")
load(file = "/mnt/data1/LewyBodyMRC/Methylation/QC/FinalData/Dasen_mSetFull.rdata")


pheno <- read.csv("QC/SummaryTables/FullQCMetrics.csv", header = T, stringsAsFactors = F)

fullMset <- fullMset[rownames(betas(fullMset)) %in% rownames(betas(mSet)),colnames(betas(fullMset)) %in% colnames(betas(mSet))]
pheno <- pheno[which(pheno$Basename %in% colnames(mSet)),]

CNGPheno <- pheno[which(pheno$Brain_Region == "CNG"),]
PFCPheno <- pheno[which(pheno$Brain_Region == "PFC"),]


plotmset_density<-function(mset, study=""){
  onetwo<-fData(mset)$DESIGN
  mat<-betas(mset)
  
  plot(density(mat[onetwo=="I",1], na.rm=T, bw=0.03), cex.main=0.8, main=paste(study, "Betas"), ylim=c(0, 5.2), xlab="")
  lines(density(mat[onetwo=="II",1], na.rm=T, bw=0.03), col="red")
  
  for(j in 2:ncol(mat)){
    lines(density(mat[onetwo=="I",j], na.rm=T, bw=0.03))
    lines(density(mat[onetwo=="II",j], na.rm=T, bw=0.03), col="red")
  }
  
  legend("topright", legend=c("Type I", "Type II"), lty=1, col=c("black", "red")) 
}

#Plot Normalised sample densities
pdf("/mnt/data1/LewyBodyMRC/Methylation/QC/Plots/Normalisation_BetaDensities.pdf", width = 10)
plotmset_density(fullMset, study="Full Raw data")
plotmset_density(mSet, study="Full Normalised Data")
plotmset_density(fullMset[,PFCPheno$Basename], study="PFC Raw data")
plotmset_density(mSetPFC, study="PFC Normalised Data")
plotmset_density(fullMset[,CNGPheno$Basename], study="CNG Raw data")
plotmset_density(mSetCNG, study="CNG Normalised Data")
dev.off()



#Check normalisation violence
Tvio <- qual(betas(fullMset),betas(mSet))
PFCvio <- qual(betas(fullMset[,PFCPheno$Basename]),betas(mSetPFC))
CNGvio <- qual(betas(fullMset[,CNGPheno$Basename]),betas(mSetCNG))

write.csv(Tvio,"/mnt/data1/LewyBodyMRC/Methylation/QC/SummaryTables/Nviolence_All.csv", row.names = F)
write.csv(PFCvio,"/mnt/data1/LewyBodyMRC/Methylation/QC/SummaryTables/Nviolence_PFC.csv", row.names = F)
write.csv(CNGvio,"/mnt/data1/LewyBodyMRC/Methylation/QC/SummaryTables/Nviolence_CNG.csv", row.names = F)

pheno$rmsd <- Tvio$rmsd
pheno[match(rownames(PFCvio),pheno$Basename),"Region_rmsd"] <- PFCvio$rmsd
pheno[match(rownames(CNGvio),pheno$Basename),"Region_rmsd"] <- CNGvio$rmsd

pdf("/mnt/data1/LewyBodyMRC/Methylation/QC/Plots/Normalisation_violence.pdf")
ggplot(pheno, aes(x = rmsd,  y = Region_rmsd))+
  geom_point()+
  xlab("rmsd (Normalised Together)")+
  ylab("rmsd (Normalised Seperately)")+
  theme_cowplot()+
  ggtitle("Normalisation violence by normalisation approach")
ggplot(pheno, aes(x = as.factor(chipNumber),y = rmsd))+
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~plateNumber)+
  theme_cowplot()+
  xlab("Chip Number")+
  ggtitle("Normalisation violence by Plate and Chip")
dev.off()


hist(Tvio$srms)
plot(CNGvio[,c(1:2)])

#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 10. Cell Type Decon/Cetygo #####
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

# Load finalised RGset
pheno <- read.csv("/mnt/data1/LewyBodyMRC/Methylation/QC/SummaryTables/FullQCMetrics.csv", header = T, stringsAsFactors = F)

FailedSamples <- pheno[which(pheno$M.median < 1500 | pheno$U.median < 1000 | 
                               pheno$bs < 80 | pheno$outliers == TRUE | pheno$crypticQC == TRUE | pheno$P.filter == TRUE),"Basename"]

pheno <- pheno[-which(pheno$Basename %in% FailedSamples),]

# idatPath<-c("/mnt/data1/Jenny/EPIC_QC/LewyBodyMRC/idats/")
# 
# RGset <- read.metharray.exp(base = idatPath, targets = pheno, force = TRUE)
# 
# minfi_props <- estimateCellCounts(rgSet = RGset, 
#                                   compositeCellType = "DLPFC",
#                                   cellTypes = c("NeuN_neg","NeuN_pos"), processMethod = "auto",  
#                                   probeSelect = "any",
#                                   referencePlatform = c("IlluminaHumanMethylation450k"), verbose = TRUE)  

#write.csv(minfi_props, "QC/SummaryTables/Minfi_CellDecon.csv")
# to evaluate the validity of minfi predicted proportions, the watermelon equivilent will be applied and tested

# Try here with the CETYGO error metric will be applied (note, predictions in this circumstance are made with watermelon)

source("/mnt/data1/LewyBodyMRC/Methylation/Scripts/ReferenceFunctions/.normalizeQuantiles2.R")
source("/mnt/data1/LewyBodyMRC/Methylation/Scripts/ReferenceFunctions/estimateCellCounts.wlmn.R")

wlmn_est <- estimateCellCounts.wmln(fullMset,compositeCellType = "DLPFC",platform = "EPIC",
                                    probeSelect = "auto",
                                    cellTypes = c("NeuN_neg","NeuN_pos"),
                                    referencePlatform = "IlluminaHumanMethylation450k")

wlmn_est <- as.data.frame(wlmn_est)

pdf("QC/Plots/CellTypePredictionError.pdf", width = 12)
ggplot(wlmn_est, aes(x = CellTypePredictionError))+
  geom_histogram(color = "white")+
  theme_cowplot()+
  ggtitle("Cell Type Prediction Error: Watermelon, CETs 2 cell, normalised")
dev.off()  

#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 11. Brain Clock ###############
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #
betas <- betas(mSet)

source("/mnt/data1/reference_files/CorticalClock/CorticalClock.r")
braincoef <- read.table("/mnt/data1/reference_files/CorticalClock/CorticalClockCoefs.txt", header = T, stringsAsFactors = F)

overlap<-braincoef[which(braincoef$probe %in% rownames(betas)),]
nrow(overlap) < nrow(braincoef)

## impute function 
imputeNA<-function(betas){
  betas[is.na(betas)]<-mean(betas,na.rm=T)
  return(betas)
}

## apply function 
betasNona<-apply(betas,2,function(x) imputeNA(x))  


## tranform betas - CpG in row


betas <-betas[braincoef$probe,]
braincoef<-braincoef[match(rownames(betas), braincoef$probe),]
brainpred<-braincoef$coef%*%betas+0.577682570446177


anti.trafo<-function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
brainpred<-anti.trafo(brainpred)


pheno<-pheno[match(colnames(betas), pheno$Basename),]
pheno$brainpred<-as.numeric(brainpred)

pdf("QC/Plots/AgePrediction.pdf")
ggplot(pheno, aes(x = Age, y = brainpred))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1, linetype = "dashed")+
  theme_cowplot()+
  xlab("Chronological age (years)")+
  ylab("Predicted age (years) - Shireby2020")
dev.off()

#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### 11. Summary / Export ###########
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #






