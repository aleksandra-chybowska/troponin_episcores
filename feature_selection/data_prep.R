library(dplyr)
library(parallel)

# Rank Inverse Based Normalisation of the data
transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

cox = read.csv('/Volumes/marioni-lab/Ola/Lab/Cox/Full_Dataset.csv')
min_set = cox[,c("ID", "assign", "agemonths", "sex", "tte", "event")]
min_set = subset(min_set, !is.na(tte) & tte>0)
min_set = na.omit(min_set)

episcores_w1_w3 = read.csv('/Volumes/marioni-lab/Ola/Lab/EpiScores/Protein_episcores/Projections/EpiScore_projections_GS_9537.csv', check.names = FALSE)
episcores_w4 = read.csv('/Volumes/marioni-lab/Ola/Lab/EpiScores/Protein_episcores/Projections/EpiScore_projections_W4_8877_220221.csv', check.names = FALSE)
names(episcores_w4)[1] = 'ID'
dim(episcores_w1_w3)
dim(episcores_w4)
identical(colnames(episcores_w1_w3), colnames(episcores_w4))
episcores = union(episcores_w1_w3, episcores_w4)
# here it is ready for the first round of the analysis - bar transformations

# divide into three sets
target = read.csv('/Volumes/marioni-lab/Ola/Lab/Test_sets/gs20ktargets.tsv', sep='\t')
episcores = episcores %>% left_join(target, by=c("ID" = "ID"))
episcores = episcores[-c(111:114)]

W1_epi = subset(episcores, Set == "W1")
W3_W4_epi = subset(episcores, Set == "W4" | Set == "W3")
count(W3_W4_epi) + count(W1_epi) == count(episcores)
#join with cox
episcores = min_set %>% merge(episcores, by.x="ID", by.y="ID")
W1 = min_set %>% merge(W1_epi, by.x="ID", by.y="ID")
W3_W4 = min_set %>% merge(W3_W4_epi, by.x="ID", by.y="ID")

head(episcores[,100:ncol(episcores)])
colnames(episcores)
episcores[7:(ncol(episcores)-1)] = mclapply(episcores[7:(ncol(episcores)-1)], transform)
head(episcores[,100:ncol(episcores)])

head(W3_W4[,100:ncol(W3_W4)])
colnames(W3_W4)
W3_W4[7:(ncol(W3_W4)-1)] = mclapply(W3_W4[7:(ncol(W3_W4)-1)], transform)
head(W3_W4[,100:ncol(W3_W4)])
head(W3_W4[,1:10])

head(W1[,100:ncol(W1)])
colnames(W1)
W1[7:(ncol(W1)-1)] = mclapply(W1[7:(ncol(W1)-1)], transform)
head(W1[,100:ncol(W1)])
head(W1[,1:10])

write.csv(episcores, '/Volumes/marioni-lab/Ola/Lab/EpiScores/Improving_assign/runs/w1_w3_w4_normalised_sep/age_sex_only/ASSIGN_and_Episcores_W1_W3_W4_joined.csv', row.names = F)
write.csv(W1, '/Volumes/marioni-lab/Ola/Lab/EpiScores/Improving_assign/ASSIGN_and_Episcores_W1.csv', row.names = F)
write.csv(W3_W4, '/Volumes/marioni-lab/Ola/Lab/EpiScores/Improving_assign/ASSIGN_and_Episcores_W3_W4_joined.csv', row.names = F)
