library(imputeTS)
library(readxl)
library(tidyverse)

## Prep x files
m2beta <- function(m) {
  beta <- 2^m/(2^m + 1)
  return(beta)
}

transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

# Read in the DNAm data
# wave 3
xtrain = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3_mvals.rds")
# wave 1
xtest = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/norm_mvals_5087.rds")

# Dimensions
dim(xtest)
# [1] 860926   5087
dim(xtrain)
# [1] 773860   4450

xtrain <- xtrain[which(rownames(xtrain) %in% rownames(xtest)),]
xtest <- xtest[which(rownames(xtest) %in% rownames(xtrain)),]

# Check dimensions again to ensure common sites are in both files (772,667)
dim(xtest)
# [1] 772667   5087
dim(xtrain)
# [1] 772667   4450

# Load EPIC array file and subset to probes common to 450k and EPIC array
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]

# Subset methylation files to those in the common annotation file
xtrain <- xtrain[rownames(xtrain) %in% rownames(common_anno),]
xtest <- xtest[rownames(xtest) %in% rownames(common_anno),]

# Check dimensions for final meth sites
dim(xtrain)
# [1] 398422   4450
dim(xtest)
# [1] 398422   5087

# Read in new variable data from Riccardo (taken from GS folder on datastore)
d1 <- read.delim("/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/GS20K_Troponin_all.PHE")
names(d1)[1] <- "Sample_Name"

target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS10k_Targets.rds")
table(is.na(d1$Troponin_T)) # here!
#FALSE  TRUE
#10201  9831

table(is.na(d1$Troponin_I))
#FALSE  TRUE
#19130   902

test <- na.omit(d1)
dim(test) # here we loose half troponin records, this needs rewritten
#[1] 10201     5

d1 <- left_join(test, target, by = "Sample_Name")
test2 <- na.omit(d1)
# 5338   17 ! Only 5330 rows left with phenotypic data. Bit less than what Danni had.
write.csv(test2, "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/yfile_all_troponin_DNAm_5338.csv", row.names = F)

W1 <- test2[which(test2$Set == "W1"),]
W3 <- test2[which(test2$Set == "W3"),]
dim(W1) #2779   17
# [1] 2779   17
dim(W3) #2559   17
# [1] 2559   17

write.csv(W1, "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/yfile_W1_DNAm_2779.csv", row.names = F)
write.csv(W3, "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/yfile_W3_DNAm_2559.csv", row.names = F)

# W1 will form ytest and W3 will form ytrain

# Transpose first so cpgs are columns and people are rows
xtr <- t(xtrain)
xte <- t(xtest)

# Now subset the respective waves of GS by the phenotype data present in this file (Sample_Sentrix_IDs for matching to DNAm data)
xtrain_sub <- xtr[which(rownames(xtr) %in% W3$Sample_Sentrix_ID),]
xtest_sub <- xte[which(rownames(xte) %in% W1$Sample_Sentrix_ID),]

# Number of people is now the same across test and train files
# > dim(xtrain_sub)
# [1]   2559 398422
# > dim(xtest_sub)
# [1]   2779 398422
# > dim(W3)
# [1] 2559   17
# > dim(W1)
# [1] 2779   17

# Match order of x files to y files
xtrain_sub <- xtrain_sub[match(W3$Sample_Sentrix_ID, rownames(xtrain_sub)),]
xtest_sub <- xtest_sub[match(W1$Sample_Sentrix_ID, rownames(xtest_sub)),]

m <- xtrain_sub
m <- m2beta(m)

#this shouldnt be neccessary here, as we don't have NAs but I'll leave it
m <- na_mean(m)

# change numbers to standard deviations (by column!)
scaled <- apply(m, 2, scale)
rownames(scaled) <- rownames(m)
xtrain <- scaled


m <- xtest_sub
m <- na_mean(m)

scaled <- apply(m, 2, scale)
rownames(scaled) <- rownames(m)
xtest <- scaled


# Regress ytrain file by age/sex/PCs and select y protein data in each case. W3 - test
list <- c(4,5)
for(i in list) {
  W3[,i] <- scale(resid(lm(W3[,i] ~ W3$age.x + factor(W3$sex.x), data=W3, na.action="na.exclude")))
}

for(i in list) {
  W3[,i] <- transform(W3[,i])
}


# here we may want to add transforming W3
# I don't have much demo data in this dataset, I'll put there only the rest of the target file
demo <- W1[c(2,3,6,11:17)]

W1 <- W1[c(6, 4, 5)]
W3 <- W3[c(6, 4, 5)]

identical(W1$Sample_Sentrix_ID, rownames(xtest)) # TRUE
identical(W3$Sample_Sentrix_ID, rownames(xtrain)) # TRUE


# write.csv(W3, "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/ytrain_W3_file.csv", row.names = F) #troponins for W3
# write.csv(W1, "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/ytest_W1_file.csv", row.names = F) #troponins for W1
# saveRDS(xtrain, "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/xtrain_W3_file.rds") # cpgs W3
# saveRDS(xtest, "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/xtest_W1_file.rds") # cpgs
# write.csv(demo, "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/demo_file.csv", row.names = F) # stuff

############################################# Prepared Dataset ####################################################

library(glmnet)
library(tidyverse)
library(foreign)

# ytrain <- read.csv("/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/ytrain_W3_file.csv")
# ytest <- read.csv("/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/ytest_W1_file.csv")
# xtrain <- readRDS("/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/xtrain_W3_file.rds")
# xtest <- readRDS("/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/xtest_W1_file.rds")
# demo <- read.csv("/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/inputs/demo_file.csv")

# alternatively
ytrain <- W3
ytest <- W1

location <- "/Cluster_Filespace/Marioni_Group/Ola/GS_added_biomarkers/elnets/outputs/newest"

xtest_copy <- xtest
x <- xtrain # cpgs W3
q <- ytrain[1] #SentrixID
p <- ytrain[2] #TroponinT
name_p <- colnames(p) # here I remember I am ananlysing troponin T
y <- cbind(q,p)
names(y)[2] = "pheno" # here I forget about it!
y <- as.numeric(y$pheno)
x <- as.matrix(x)

set.seed(1234)

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


xtest = xtest_copy

x <- xtrain # cpgs W3
q <- ytrain[1] #SentrixID
p <- ytrain[3] #TroponinT
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
p2 <- ytest[3] # Get the protein data for the iteration of interest from the y variable
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


