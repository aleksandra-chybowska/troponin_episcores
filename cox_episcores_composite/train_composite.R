library("optparse")
library("glmnet")
library("survival")
library("jsonlite")

transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

# Get arguments
######################################################

option_list = list(
  make_option(c("-s", "--settings"), type="character", default=NULL, 
              help="Settings file path (settings.json)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

url = '/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/composite_cox_settings.json'

if (!is.null(opt$settings)) {
  url = opt$settings
}

settings <- fromJSON(txt=url, flatten = FALSE)


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

if (settings$include_set == 1) {
  dummy_set = model.matrix( ~ set - 1, df)
  df = data.frame(df[ ,!colnames(df) %in% "set"], dummy_set)
} else {
  df = subset(df, select = -c(set))
}

# Variables to keep: age, sex, grimage components, and 109 episcores (everything but dead status and tte, first two columns)
x = subset(df, select = -c(Sentrix_ID, id, assign, dead, event, tte))
x = as.matrix(x)
y = Surv(df$tte,df$event) # Time to event, and wether the event has happened or not
cat("Imported and prepped data!\n")


# Get CV'd lambda and obtain effects of each CpG
######################################################

cat("Fitting elastic net model (Cox PH) for time to CVD...\n")
# Obtain lambda
cv <- cv.glmnet(x, y, family = "cox", type.measure = "C", seed = seed, nfolds = 10) # Harrell's C index to obtain best parameters (kind of like residuals to minimize in OLS)
lambda <- cv$lambda.min
# Obtain coefs
fit <- glmnet(x, y, family = "cox", lambda = lambda)

# Get the good stuff!
######################################################

cat("Now extracting info of interest.\n")
coefs <- coef(fit) # Get betas
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
coefs["Variable"] <- rownames(coefs)
names(coefs)[1] <- "Coefficient"
coefs <- coefs[,c("Variable", "Coefficient")]

# Export 
######################################################

cat("Exporting!\n")
filename <- paste0(settings$output, "elnet_coefficients_", settings$run, ".csv")
#filename <- paste0("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/elasticnet_models/elnet_train/w1w3/random/", "elnet_coefficients_random", ".tsv")
write.csv(coefs, file = filename, row.names = FALSE, quote = FALSE)

