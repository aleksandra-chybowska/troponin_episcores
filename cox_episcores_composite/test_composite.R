#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Ola
library("survival")
library("optparse")
library("stringr")
library("imputeTS")
library("jsonlite")
library("dplyr")
box::use(../../modules/transformations[...])

# Get arguments
######################################################

option_list = list(
  make_option(c("-s", "--settings"), type="character", default=NULL, 
              help="Settings file path (settings.json)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

url = '/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/test/test_composite_settings.json'

if (!is.null(opt$settings)) {
  url = opt$settings
}

settings <- fromJSON(txt=url, flatten = FALSE)

# Import + prep data
######################################################

set.seed(1234) # Set seed to ensure fold variation minimised 
seed <- 1234

# folds <- read.delim(opt$folds, header = TRUE, row.names = 2)
df <- read.csv(settings$input, check.names = FALSE)
# Remove people that have missing or strange time-to-event values (negative)
df = df[!is.na(df$tte) & df$tte>0, ]
df = na.omit(df)

if (settings$transform == "rank") {
  df$assign = transform(df$assign)
} else if (settings$transform == "log") {
  df$assign = log(df$assign + 1)
}

rownames(df) = df$Sentrix_ID

# Variables to keep: age, sex, grimage components, and 109 episcores (everything but dead status and tte, first two columns)
x = subset(df, select = -c(Sentrix_ID, id, assign, set, dead, event, tte))
x = as.matrix(x)


weights <-read.csv(settings$scores, check.names = FALSE)
cat("\nHave imported weights for CpGs.\n")

episcores = t(x)
episcores = episcores[which(rownames(episcores) %in% weights$Variable),]
episcores = episcores[match(weights$Variable, rownames(episcores)),] 
calc = episcores * weights$Coefficient
sum = as.data.frame(colSums(calc)) # Sum the score for each person (column) in the methylation dataset 
export_sum = data.frame("sample" = rownames(sum) , "summed_scores" = sum, "trait" = "EpiScore")
colnames(export_sum) = c("sample", "summed_scores", "trait")


# Output
######################################################

write.csv(export_sum, paste0(settings$output, "export_sum_", settings$run, ".csv"), row.names = F)

# Performance Check
######################################################

x = df[c("age", "sex", "assign", "tte", "event")]
x$summed_scores = export_sum$summed_scores

riskFactorsOnlyCoxPH = coxph(Surv(tte, event) ~ age + sex, x)
proteinCoxPH = coxph(Surv(tte, event) ~ age + sex + assign, x)

models = list(r = riskFactorsOnlyCoxPH, p = proteinCoxPH)

# GENERATE INCREMENTAL R2

predictCoxPHOnset <- function(dataDF, coxPHModel, threshold = 12) {
  uniqueTimes <- sort(unique(dataDF$tte))
  thresholdIndex <- match(threshold, uniqueTimes)
  cumulativeBaseHaz <- gbm::basehaz.gbm(dataDF$tte, dataDF$event, predict(coxPHModel), uniqueTimes)
  survivalPredictions <- exp(-cumulativeBaseHaz[[thresholdIndex]]) ^ exp(predict(coxPHModel))# check what this is doing
  onsetPredictions <- 1 - survivalPredictions
  
  # event should be 0 if tte is > 10
  dataDF$event <- sapply(1:nrow(dataDF), function(i) {
    if (dataDF$tte[[i]] > threshold) {
      0
    } else {
      dataDF$event[[i]]
    }
  })
  
  auc <- MLmetrics::AUC(y_pred = onsetPredictions, y_true = dataDF$event)
  prauc <- MLmetrics::PRAUC(y_pred = onsetPredictions, y_true = dataDF$event)
  roc <- pROC::roc(response = dataDF$event, predictor = onsetPredictions)
  list(cumulativeBaseHaz = cumulativeBaseHaz, onsetPredictions = onsetPredictions, auc = auc, prauc = prauc, roc = roc)
}

testResults <- lapply(models, function(m) {predictCoxPHOnset(x, m)})

aucs <- sapply(testResults, function(r) {r$auc})
praucs <- sapply(testResults, function(r) {r$prauc})
metricsTable <- data.frame(AUC = aucs, PRAUC = praucs)
metricsTable
