library(dplyr)
library(parallel)

box::use(../../modules/transformations[...])

cox = read.csv('/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/cox_covars.csv')
min_set = subset(cox, !is.na(tte) & tte>0)
min_set = na.omit(min_set)

episcores_w1_w3 = read.csv('/Volumes/marioni-lab/Ola/Lab/EpiScores/Protein_projections/EpiScore_projections_GS_9537.csv', check.names = FALSE)
episcores_w4 = read.csv('/Volumes/marioni-lab/Ola/Lab/EpiScores/Protein_projections/EpiScore_projections_W4_8877_220221.csv', check.names = FALSE)
names(episcores_w4)[1] = 'ID'
identical(colnames(episcores_w1_w3), colnames(episcores_w4))
episcores = union(episcores_w1_w3, episcores_w4)
merged = merge(min_set, episcores, by.x="Sentrix_ID", by.y="ID")

#join with cox
W1 = subset(merged, set == "W1")
W3_W4 = subset(merged, set == "W4" | set == "W3")

#keep only unrelated
unrelated = read.csv('/Volumes/marioni-lab/Ola/Lab/Test_sets/Unique_W3_W4_and_W1.csv')
W3_W4 = W3_W4[which(W3_W4$id %in% unrelated$Sample_Name),]

start = length(min_set) + 1
end = ncol(merged)

head(W3_W4[,1:10])
head(W3_W4[,100:ncol(W3_W4)])
W3_W4[start:end] = mclapply(W3_W4[start:end], transform)
head(W3_W4[,100:ncol(W3_W4)])
head(W3_W4[,1:10])

head(W1[,1:10])
head(W1[,100:ncol(W1)])
W1[start:end] = mclapply(W1[start:end], transform)
head(W1[,100:ncol(W1)])
head(W1[,1:10])

write.csv(W1, '/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/cox_covars_episcores_W1.csv', row.names = F)
write.csv(W3_W4, '/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/cox_covars_episcores_W3_W4.csv', row.names = F)
