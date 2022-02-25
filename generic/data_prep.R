library("optparse")
library("glue")
library("imputeTS")
library("dplyr")
library("clusterSim")
library("jsonlite")

# Get arguments
######################################################

option_list = list(
  make_option(c("-s", "--settings"), type="character", default=NULL, 
              help="Settings file path (settings.json)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

url = '/Users/shirin/Projects/R/troponin_episcores/generic/settings/test_settings/settings_local.json'

if (!is.null(opt$settings)) {
  url = opt$settings
}

settings <- fromJSON(txt=url, flatten = FALSE)

cat(settings$transform_test)
quit()

cat("\n#### Dimensions - test:\n")
test = list()
test_datasets = names(settings$test)
for(dataset in test_datasets) {
  test[[dataset]] = readRDS(settings$test[[dataset]])
  cat(paste(dataset, "\n", sep = ""))
  cat(paste("\t", dim(test[[dataset]])))
  cat("\n")
}

cat("\n#### Dimensions - train:\n")
train=list()
train_datasets = names(settings$train)
for(dataset in train_datasets) {
  train[[dataset]] = readRDS(settings$train[[dataset]])
  cat(paste(dataset, "\n", sep = ""))
  cat(paste("\t", dim(train[[dataset]])))
  cat("\n")
}

target = read.csv(settings$target, sep='\t')
unique = read.csv(settings$cpg_subset)
pheno = read.csv(settings$pheno)
feature = settings$feature

# Some functions
######################################################

m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}

transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

# Remove related (train/test)
######################################################

cat("\n#### Remove related - train:\n")
unrelated = unique %>% left_join(target, by="Sample_Name")

for(dataset in train_datasets) {
    train[[dataset]] = train[[dataset]][,colnames(train[[dataset]]) %in% unrelated$ID]
    cat(paste(dataset, "\n", sep = ""))
    cat(paste("\t", dim(train[[dataset]])))
    cat("\n")
}

# Read in pheno + fold data
######################################################

cat("\n#### Import pheno data:\n")
cat("Before removing NA and reducing col no:\n")
cat(paste("\t", dim(pheno)))

# pheno = pheno %>% left_join(target, by=c("id" = "Sample_Name"))
# pheno = pheno[,c("id", "ID", "sex.y", "age.y", "Batch", feature)]
# pheno = na.omit(pheno)
# names(pheno) = c("id", "sentrix_ID", "sex", "age", "batch", feature)

pheno = na.omit(pheno)
#Remove related from pheno too

cat("After filtering:\n")
cat(paste("\t", dim(pheno)))

pheno_train = list()
for(train_ids in settings$train_identifiers) {
  pheno_train[[train_ids]] = subset(pheno, pheno$set == train_ids)
}
pheno_train = do.call("rbind", pheno_train)
pheno_train = pheno_train[pheno_train$Sample_Sentrix_ID %in% unrelated$ID,]

cat("Train after filtering:\n")
cat(paste("\t", dim(pheno_train)))

pheno_test = list()
for(test_ids in settings$test_identifiers) {
  pheno_test[[test_ids]] = subset(pheno, pheno$set == test_ids)
}

#SPRAWDZIC!
pheno_test = do.call("rbind", pheno_test)

if (settings$transform_test == "log") {
  cat("Log transforming test!")
  pheno_test[feature] = log(pheno_test[feature])
} 
if (settings$transform_test == "log+1") {
  cat("Log+1 transforming test!")
  pheno_test[feature] = log(pheno_test[feature]+1)
}
if (settings$transform_test == "rank") {
  cat("Rank transforming test!")
  pheno_test[feature] = transform(pheno_test[feature]+1)
}

cat("Test after filtering:\n")
cat(paste("\t", dim(pheno_test)))

cat("\nDivided into train and test. Imported pheno data.\n")
rm(pheno)
gc()

# Clean up methylation data a lil'
######################################################
# Get common probes and keep common
cat("\n#### DNAm filtering: common probes in train + test; EPIC filtering \n")

probes_wave = list()
for(dataset in train_datasets) {
  probes_wave[[dataset]] = rownames(train[[dataset]])
}
for(dataset in test_datasets) {
  probes_wave[[dataset]] = rownames(test[[dataset]])
}
common <- Reduce(intersect, probes_wave)

# for (dataset in train_datasets) {
#   train[[dataset]] <- train[[dataset]][which(rownames(train[[dataset]]) %in% common),]
# }
# w3_w4 <- do.call("cbind", train)
# saveRDS(w3_w4, "/Volumes/marioni-lab/Ola/Lab/Test_sets/w3_w4_local.rds")

cat("Train+Test: Common CpGs:\n")
cat(paste("\t", length(common)))

# Get EPIC array file and subset to probes common to 450k and EPIC array
epic = readRDS(settings[["epic"]])
common_array <- rownames(epic[which(epic$Methyl450_Loci == "TRUE"),])
rm(epic)
cat("\nEPIC RAM clean up...\n\n")
gc()

# Filter (also making sure sample in pheno file)
for (dataset in train_datasets) {
  train[[dataset]] <- train[[dataset]][which(rownames(train[[dataset]]) %in% common),]
  train[[dataset]] <- train[[dataset]][which(rownames(train[[dataset]]) %in% common_array),]
  train[[dataset]] <- train[[dataset]][,which(colnames(train[[dataset]]) %in% pheno_train$Sample_Sentrix_ID)]
}

# Filter (also making sure sample in pheno file)
for (dataset in test_datasets) {
  test[[dataset]] <- test[[dataset]][which(rownames(test[[dataset]]) %in% common),]
  test[[dataset]] <- test[[dataset]][which(rownames(test[[dataset]]) %in% common_array),]
  test[[dataset]] <- test[[dataset]][,which(colnames(test[[dataset]]) %in% pheno_test$Sample_Sentrix_ID)]
}

cat("\nKept only probes in EPIC array and those common to all waves.\n")
# Print dimensions to make sure everything is alright
cat("\n#### Train Dimensions:\n")
for (dataset in train_datasets) {
  cat(paste(dataset, "\n", sep = ""))
  cat(paste("\t", dim(train[[dataset]])))
  cat("\n")
}

cat("\n#### Test Dimensions:\n")
for (dataset in test_datasets) {
  cat(paste(dataset, "\n", sep = ""))
  cat(paste("\t", dim(test[[dataset]])))
  cat("\n")
}

cat("\nTranspose datasets.\n")
for (dataset in train_datasets) {
  train[[dataset]] <- t(train[[dataset]])
  gc()
}

for (dataset in test_datasets) {
  test[[dataset]] <- t(test[[dataset]])
  gc()
}

# Convert to beta values, impute NAs, match order
######################################################
dataset = train_datasets[1]
cols <- colnames(train[[dataset]])

for (dataset in train_datasets) {
  train[[dataset]] = m2beta(train[[dataset]])
  train[[dataset]] = na_mean(train[[dataset]])
  rownames = rownames(train[[dataset]])
  train[[dataset]] = apply(train[[dataset]], 2, scale)
  rownames(train[[dataset]]) = rownames
  train[[dataset]] = train[[dataset]][,cols]
  cat("\nRAM clean up...\n\n")
  gc()
}

for (dataset in test_datasets) {
  test[[dataset]] = m2beta(test[[dataset]])
  test[[dataset]] = na_mean(test[[dataset]])
  rownames = rownames(test[[dataset]])
  test[[dataset]] = apply(test[[dataset]], 2, scale)
  rownames(test[[dataset]]) = rownames
  test[[dataset]] = test[[dataset]][,cols]
  cat("\nRAM clean up...\n\n")
  gc()
}

# > identical(colnames(train[["wave3"]]), colnames(train[["wave4"]]))
# [1] TRUE
# > identical(colnames(train[["wave3"]]), colnames(test[["wave1"]]))
# [1] TRUE

# Fuse
######################################################

train_df <- do.call("rbind", train)
cat("\n#### Dimensions after fusing waves - train:\n")
cat(paste("\t", dim(train_df), "\n"))

rm(train)
cat("\nRAM clean up...\n\n")
gc()

train_df <- train_df[match(pheno_train$Sample_Sentrix_ID, rownames(train_df)),] # Match order of methylation table (x) and phenotype table (y) 
cat("\nRAM clean up...\n\n")
gc()
identical(pheno_train$Sample_Sentrix_ID, rownames(train_df))

test_df = do.call("rbind", test)
cat("\n#### Dimensions after fusing waves - test:\n")
cat(paste("\t", dim(test_df), "\n"))
rm(test)

test_df <- test_df[match(pheno_test$Sample_Sentrix_ID, rownames(test_df)),] # Match order of methylation table (x) and phenotype table (y) 
identical(pheno_test$Sample_Sentrix_ID, rownames(test_df))
cat("\nRAM clean up...\n\n")
gc()
# > dim(train_df)
# [1] 9846  713
# > dim(test_df)
# [1] 5080  713
# > dim(pheno_train)
# [1] 9846    6
# > dim(pheno_test)
# [1] 5080    6


# Export data
######################################################

cat("\nExporting prepped data for models...\n")
o_name_rds <- paste0(opt$out, "methylation_training_", opt$name, ".rds")
saveRDS(train_df, o_name_rds, compress = FALSE)

# 
write.csv(pheno_train, settings$o_pheno_train, row.names = F) #troponins for W3
write.csv(pheno_test, settings$o_pheno_test, row.names = F) #troponins for W1
saveRDS(train_df, settings$o_train_df) # cpgs W3
saveRDS(test_df, settings$o_test_df) # cpgs


# If external, export just external data without GS
######################################################
# 
# (!is.null(opt$lbc)) {
#   pheno_noGS <- pheno[!(pheno$cohort %in% c("W1", "W3", "W4")),]
#   train_df <- train_df[rownames(train_df) %in% rownames(pheno_noGS),]
#   
#   cat("\nExporting prepped data for models (just external)...\n")
#   o_name_rds <- paste0(opt$out, "methylation_training_", opt$name, "_noGS.rds")
#   saveRDS(train_df, o_name_rds, compress = FALSE)
# }