####################################################################################

### Elastic Net for the 2 cardiac troponins in GS to train DNAm episcores 

####################################################################################

# Developed by Ola, based on Danni's Episcores script
# Generate episcores for cTnT and cTnI, compare to measured levels

####################################################################################

library(readxl)
library(tidyverse)
library(imputeTS)
library(ggplot2)

m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}

transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

#annotated methylation sites

xtrain = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3_mvals.rds")
xtest = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/norm_mvals_5087.rds")

dim(xtrain)
dim(xtest)

xtrain <- xtrain[which(rownames(xtrain) %in% rownames(xtest)),]
xtest <- xtest[which(rownames(xtest) %in% rownames(xtrain)),]

dim(xtrain)
dim(xtest)

anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]

xtrain <- xtrain[rownames(xtrain) %in% rownames(common_anno),]
xtest <- xtest[rownames(xtest) %in% rownames(common_anno),]

dim(xtrain)
dim(xtest)

#phenotypic (biomarker + patient data)

d1 <- read.delim("/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/GS20K_Troponin_all.PHE")
names(d1)[1] <- "Sample_Name"

target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS10k_Targets.rds")

table(is.na(d1$Troponin_T))
# FALSE  TRUE 
# 10201  9831 

table(is.na(d1$Troponin_I))
# FALSE  TRUE 
# 19130   902 

list <- c(4,5)

i=4
troponin = na.omit(d1[,c(1,2,3,i)])
protein = colnames(d1)[i]
dim(troponin) # this condiders only the rows where both cTnI and cTnT are not null
joined_troponin <- left_join(troponin, target, by = "Sample_Name")
test2 = na.omit(joined_troponin)
dim(test2)

location = "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/yfile"
write.csv(test2, file = paste0(location, "_", protein, "_", as.numeric(count(test2)), "_DNAm_pheno.csv"), row.names = F)

W1 <- test2[which(test2$Set == "W1"),]
W3 <- test2[which(test2$Set == "W3"),]

dim(W1) #2779   17
dim(W3) #2559   17

write.csv(W1, file = paste0(location, "_W1_", protein, "_", as.numeric(count(W1)), "_DNAm_pheno.csv"), row.names = F)
write.csv(W3, file = paste0(location, "_W3_", protein, "_", as.numeric(count(W3)), "_DNAm_pheno.csv"), row.names = F)
  
# Transpose first so cpgs are columns and people are rows
xtr <- t(xtrain)
xte <- t(xtest)

# Now subset the respective waves of GS by the phenotype data present in this file (Sample_Sentrix_IDs for matching to DNAm data)
xtrain_sub_raw <- xtr[which(rownames(xtr) %in% W3$Sample_Sentrix_ID),]
xtest_sub_raw <- xte[which(rownames(xte) %in% W1$Sample_Sentrix_ID),]

# match order of rows in xtrain_sub to the one in phenotypic data which have been subset to w3 or w1
xtrain_sub <- xtrain_sub_raw[match(W3$Sample_Sentrix_ID, rownames(xtrain_sub_raw)),]
xtest_sub <- xtest_sub_raw[match(W1$Sample_Sentrix_ID, rownames(xtest_sub_raw)),]

# convert to beta values - why? They are scaled soon!
m <- xtrain_sub
m <- m2beta(m)

#this shouldnt be neccessary here, as we don't have NAs but I'll leave it 
m <- na_mean(m)

# change numbers to standard deviations (by column!)
scaled <- apply(m, 2, scale)
rownames(scaled) <- rownames(m)
xtrain_scaled <- scaled

m <- xtest_sub
m <- na_mean(m)

scaled <- apply(m, 2, scale)
rownames(scaled) <- rownames(m)
xtest_scaled <- scaled


W3[,4] <- scale(resid(lm(W3[,4] ~ W3$age.x + factor(W3$sex.x), data=W3, na.action="na.exclude")))
W3[,4] <- transform(W3[,4])
W1[,4] <- transform(W1[,4])

demo <- W1[c(2, 3, 5, 10:16)]

W1 <- W1[c(5, 4)]
W3 <- W3[c(5, 4)]

identical(W1$Sample_Sentrix_ID, rownames(xtest_scaled)) # TRUE
identical(W3$Sample_Sentrix_ID, rownames(xtrain_scaled)) # TRUE

inputs_path = "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/"
  
write.csv(W3, paste0(inputs_path, "ytrain_W3_", protein, ".csv"), row.names = F) #troponins for W3
write.csv(W1, paste0(inputs_path, "ytest_W1_", protein, ".csv"), row.names = F) #troponins for W1
saveRDS(xtrain_scaled, paste0(inputs_path, "xtrain_W3_", protein, ".rds")) # cpgs W3
saveRDS(xtest_scaled, paste0(inputs_path, "xtest_W1_", protein, ".rds")) # cpgs
write.csv(demo, paste0(inputs_path, "demo_", protein, ".csv"), row.names = F)

i=5
troponin = na.omit(d1[,c(1,2,3,i)])
protein = colnames(d1)[i]
dim(troponin) # this condiders only the rows where both cTnI and cTnT are not null
joined_troponin <- left_join(troponin, target, by = "Sample_Name")
test2 = na.omit(joined_troponin)
dim(test2)

location = "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/yfile"
write.csv(test2, file = paste0(location, "_", protein, "_", as.numeric(count(test2)), "_DNAm_pheno.csv"), row.names = F)

W1 <- test2[which(test2$Set == "W1"),]
W3 <- test2[which(test2$Set == "W3"),]

dim(W1) #2779   17
dim(W3) #2559   17

write.csv(W1, file = paste0(location, "_W1_", protein, "_", as.numeric(count(W1)), "_DNAm_pheno.csv"), row.names = F)
write.csv(W3, file = paste0(location, "_W3_", protein, "_", as.numeric(count(W3)), "_DNAm_pheno.csv"), row.names = F)

# Transpose first so cpgs are columns and people are rows
xtr <- t(xtrain)
xte <- t(xtest)

# Now subset the respective waves of GS by the phenotype data present in this file (Sample_Sentrix_IDs for matching to DNAm data)
xtrain_sub_raw <- xtr[which(rownames(xtr) %in% W3$Sample_Sentrix_ID),]
xtest_sub_raw <- xte[which(rownames(xte) %in% W1$Sample_Sentrix_ID),]

# match order of rows in xtrain_sub to the one in phenotypic data which have been subset to w3 or w1
xtrain_sub <- xtrain_sub_raw[match(W3$Sample_Sentrix_ID, rownames(xtrain_sub_raw)),]
xtest_sub <- xtest_sub_raw[match(W1$Sample_Sentrix_ID, rownames(xtest_sub_raw)),]

# convert to beta values - why? They are scaled soon!
m <- xtrain_sub
m <- m2beta(m)

#this shouldnt be neccessary here, as we don't have NAs but I'll leave it 
m <- na_mean(m)

# change numbers to standard deviations (by column!)
scaled <- apply(m, 2, scale)
rownames(scaled) <- rownames(m)
xtrain_scaled <- scaled

m <- xtest_sub
m <- na_mean(m)

scaled <- apply(m, 2, scale)
rownames(scaled) <- rownames(m)
xtest_scaled <- scaled



W3[,4] <- scale(resid(lm(W3[,4] ~ W3$age.x + factor(W3$sex.x), data=W3, na.action="na.exclude")))
W3[,4] <- transform(W3[,4])
W1[,4] <- transform(W1[,4])

demo <- W1[c(2, 3, 5, 10:16)]

W1 <- W1[c(5, 4)]
W3 <- W3[c(5, 4)]

identical(W1$Sample_Sentrix_ID, rownames(xtest_scaled)) # TRUE
identical(W3$Sample_Sentrix_ID, rownames(xtrain_scaled)) # TRUE

inputs_path = "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/"

write.csv(W3, paste0(inputs_path, "ytrain_W3_", protein, ".csv"), row.names = F) #troponins for W3
write.csv(W1, paste0(inputs_path, "ytest_W1_", protein, ".csv"), row.names = F) #troponins for W1
saveRDS(xtrain_scaled, paste0(inputs_path, "xtrain_W3_", protein, ".rds")) # cpgs W3
saveRDS(xtest_scaled, paste0(inputs_path, "xtest_W1_", protein, ".rds")) # cpgs
write.csv(demo, paste0(inputs_path, "demo_", protein, ".csv"), row.names = F)


# Data prepared
################################################################################

library(glmnet)
library(tidyverse)
library(foreign)

inputs_path = "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/"
location <- "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/outputs/age_sex_regressed/newest"

# run for troponin T
analysis <- function(xtrain, ytrain, xtest, ytest, demo) {
  # Train weights 
  x <- xtrain # cpgs W3
  q <- ytrain[1] #SentrixID
  p <- ytrain[2] #TroponinT
  name_p <- colnames(p) # here I remember I am ananlysing troponin T
  y <- cbind(q,p)
  names(y)[2] = "pheno" # here I forget about it!
  y <- as.numeric(y$pheno)
  x <- as.matrix(x)
  
  lasso.cv <- cv.glmnet(x, y, family="gaussian", alpha = 0.5, nfolds=10)
  #lasso.cv$lambda.min
  #[1] 0.3801759
  fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, standardize = F, lambda = lasso.cv$lambda.min) # model fit 
  
  coefs <- coef(fit) # Extract coeficients 
  coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
  coefs$Predictor <- name_p # Assign protein identifier
  names(coefs)[1] <- "Coefficient" # Tidy naming 
  coefs$EpiScore <- rownames(coefs) # Create episcores column
  coefs <- coefs[c(3,1,2)] # order 
  coefs3 <- coefs[-1,] # Remove intercept (if there was only intercept, this will now have 0 rows)
  write.csv(coefs3, file = paste0(location, "_", name_p, "_weights.csv"), row.names = F)
  
  # GENERATE SCORES
  q2 <- ytest[1] # Get just basenames for people in the y variable 
  p2 <- ytest[2] # Get the protein data for the iteration of interest from the y variable 
  name_p <- colnames(p2) # Get the name of the protein for this iteration
  y2 <- cbind(q2,p2) # Bind Basename and protein data together into one set 
  names(y2)[2] <- "pheno" # Assign a generic name to the protein variable
  xtest <- t(xtest) # transpose so rows are the episcores and people are columns 
  overlap <- which(rownames(xtest) %in% rownames(coefs3)) # Find the overlap between CpGs in the predictor weights column and the CpGs in the test methylation data 
  xtest <- xtest[overlap,] # Subset methylation CpG sites based on this overlap 
  match <- xtest[match(rownames(coefs3), rownames(xtest)),] # Match up the order of CpGs in methylation file to those in the CpG predictor weights column
  calc <- match * coefs3[,2] # Multiply beta predictor weights to the CpG methylation values for each person (column) in the methylation dataset 
  sum <- colSums(calc) # Sum the score for each person (column) in the methylation dataset 
  export_sum <- as.data.frame(sum) # Get the scores ready to be written to file
  names(export_sum)[1] <- "Scores"
  export_sum$Predictor <- name_p
  export_sum$Sample_Name <- rownames(export_sum)
  write.csv(export_sum, file = paste0(location, "_", name_p, "_scores.csv"), row.names = F)
  
  # GENERATE INCREMENTAL R2
  names(export_sum)[3] <- "Sample_Sentrix_ID"
  join <- left_join(export_sum, ytest, by = "Sample_Sentrix_ID")
  join <- left_join(join, demo, by = "Sample_Sentrix_ID")
  print(names(join)[4])
  names(join)[4] <- "outcome"
  null <- summary(lm(outcome ~ age.x + sex.x, data=join))$r.squared
  full <- summary(lm(outcome ~ age.x + sex.x + Scores, data=join))$r.squared
  print(round(100*(full - null), 3))
}

studied_proteins = c("Troponin_T", "Troponin_I")
set.seed(1234) # set seed to ensure fold variation minimised 

i = 1
ytrain_t <- read.csv(paste0(inputs_path, "ytrain_W3_", studied_proteins[i], ".csv"))
ytest_t <- read.csv(paste0(inputs_path, "ytest_W1_", studied_proteins[i], ".csv"))
xtrain_t <- readRDS(paste0(inputs_path, "xtrain_W3_", studied_proteins[i], ".rds"))
xtest_t <- readRDS(paste0(inputs_path, "xtest_W1_", studied_proteins[i], ".rds"))
demo_t <- read.csv(paste0(inputs_path, "demo_", studied_proteins[i], ".csv"))
analysis(xtrain = xtrain_t, ytrain = ytrain_t, xtest = xtest_t, ytest = ytest_t, demo = demo_t)


i = 2 
ytrain_t <- read.csv(paste0(inputs_path, "ytrain_W3_", studied_proteins[i], ".csv"))
ytest_t <- read.csv(paste0(inputs_path, "ytest_W1_", studied_proteins[i], ".csv"))
xtrain_t <- readRDS(paste0(inputs_path, "xtrain_W3_", studied_proteins[i], ".rds"))
xtest_t <- readRDS(paste0(inputs_path, "xtest_W1_", studied_proteins[i], ".rds"))
demo_t <- read.csv(paste0(inputs_path, "demo_", studied_proteins[i], ".csv"))
analysis(xtrain = xtrain_t, ytrain = ytrain_t, xtest = xtest_t, ytest = ytest_t, demo = demo_t)

