library(dplyr)

#box::use(../../modules/cox[...])
make_cox_dataset = function(ds) {
  
  ds[which(is.na(ds$mob)), "mob"] = 1
  # there is another file - targets file - where this has been corrected
  
  ds$dead = 0
  # add natural deaths
  ds[which(!is.na(ds$dod_ym)), "dead"] = 1
  
  #check for multiple events!!!
  ds$event = 0
  ds[which(!is.na(ds$Incident.Event)), "event"] = 1
  
  # year and month of event
  ds$yoe = as.numeric(substring(ds$Event_Date.Event, 1, 4))
  ds$moe = as.numeric(substring(ds$Event_Date.Event, 5, 6))
  
  # year and month of death
  ds$yod = as.numeric(substring(ds$Event_Date.Death, 1, 4))
  ds$mod = as.numeric(substring(ds$Event_Date.Death, 5, 6))
  
  #censoring
  ds$yr_diff = ifelse(ds$event==0, 2019 - ds$yob, ds$yoe - ds$yob)
  ds$m_diff = ifelse(ds$event==0, ((9 - ds$mob)/12), ((ds$moe - ds$mob)/12))
  ds$age_event = ds$yr_diff + ds$m_diff
  ds$tte = ds$age_event - ds$age
  ds$tte = ifelse(ds$tte < -1, NA, ds$tte)
  ds$tte = ifelse(ds$tte < 0, 0, ds$tte)
  
  return(ds)
}

####################################### Input data #######################################

agesex = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/2022-03-10_age_alive.csv');
agemonths = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/agemonths.csv');
cvd_deaths = read.csv('/Volumes/marioni-lab/Ola/Lab/Cox/2022-01-18/cvd_deaths_withcode.csv'); # 696
non_cvd_deaths = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/2021-09-27_death.csv');
hosp_events = read.csv('/Volumes/marioni-lab/Ola/Lab/Cox/2022-01-18/cvd_hospitalisations_and_operations_withcode.csv') # 3666
assign =  read.csv('/Volumes/marioni-lab/Ola/Lab/ASSIGN/GS_assign_data.csv')
target = readRDS("/Volumes/marioni-lab/Ola/Lab/Test_sets/GS20k_Targets.rds")

####################################### CVD subsets #######################################

events = hosp_events
events = subset(events, Incident==1 & Event_Type != "GP" & Event_Type != "Death") # 1297

main_ds = assign %>%
  left_join(agesex, c("id" = "id")) %>%
  left_join(target, c("id" = "Sample_Name")) %>%
  left_join(non_cvd_deaths, c("id" = "id")) %>%
  left_join(cvd_deaths, c("id" = "id"))

ds = main_ds %>% left_join(events, c("id" = "id"))

ds = ds[c("id", "Sample_Sentrix_ID", "sex.x", "age.x", "assign", "Troponin_T", "Troponin_I", "cTnI_corrected", "Set", "yob", "mob", "dod_ym", 
          "Incident.x", "Event_Type.x", "Event_Date.x", "Incident.y", "Event_Type.y", "Event_Date.y")] 
colnames(ds) = c("id", "Sentrix_ID", "sex", "age", "assign", "Troponin_T", "Troponin_I", "cTnI_corrected", "set", "yob", "mob", "dod_ym", 
                 "Incident.Death", "Event_Type.Death", "Event_Date.Death", "Incident.Event", "Event_Type.Event", "Event_Date.Event")


####################################### Cox dataset #######################################

cox = make_cox_dataset(ds)
cox = cox[c("id", "Sentrix_ID", "sex", "age", "assign", "Troponin_T", "Troponin_I", "cTnI_corrected", "set", "dead", "event", "tte")] 
cox = subset(cox, !is.na(tte) & tte>0 & tte<20)
cox = subset(cox, !is.na(assign))

# write.csv(cox, '/Volumes/marioni-lab/Ola/Lab/EpiScores/For_Riccardo/cox_covars.csv', row.names = F)

assign_cox_var = transform(cox$assign)
mod1 = coxph(Surv(tte, event) ~ age + sex + transform(cox$assign), data=cox)

ggcoxfunctional(mod1, ylim = c(-0.5,1))
martingale = residuals(mod1, type="martingale")

ggcoxdiagnostics(mod1, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())

ggcoxdiagnostics(mod1, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())

residual_test = cox.zph(mod1)
residual_test
ggcoxzph(residual_test)

#write.csv(cox, '/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/runs/episcores_only_13k_na_corrected/train/cox_covars.csv', row.names = F)

threshold = 12
cox$event <- sapply(1:nrow(cox), function(i) {
  if (cox$tte[[i]] > threshold) {
    0
  } else {
    cox$event[[i]]
  }
})


mod1 <- survfit(Surv(tte, event) ~ factor(sex), data = cox)
# -loglog plot
plot(mod1, fun = "cloglog", xlab = "Time (in days) using log",
     ylab = "log-log survival", main = "log-log curves by clinic")

ggsurvplot(mod1,
           conf.int=TRUE, # add confidence intervals
           pval=F, # show the p-value for the log-rank test
           risk.table=F, # show a risk table below the plot
           legend.labs=c("Female", "Male"), # change group labels
           legend.title="Sex",  # add legend title
           palette=c("dodgerblue4", "orchid2"), # change colors of the groups
           title="Kaplan-Meier Curve for incident CVD", # add title to plot
           risk.table.height=.2)


cox$at_risk = 0
cox[which(cox$assign>20), "at_risk"] = 1

mod1 <- survfit(Surv(tte, event) ~ factor(at_risk), data = cox)
ggsurvplot(mod1,
           conf.int=TRUE, # add confidence intervals
           pval=FALSE, # show the p-value for the log-rank test
           risk.table=FALSE, # show a risk table below the plot
           legend.labs=c("N", "Y"), # change group labels
           legend.title="ASSIGN > 20",  # add legend title
           palette=c("dodgerblue4", "orchid2"), # change colors of the groups
           title="Kaplan-Meier Curve for incident CVD", # add title to plot
           risk.table.height=.2)


cox$age_risk = 0
top_25 = max(cox$age) * 0.75
bottom_25 = max(cox$age) * 0.75
cox[which(cox$age>top_25), "age_risk"] = 1


mod1 <- survfit(Surv(tte, event) ~ factor(age_risk), data = cox)
ggsurvplot(mod1,
           conf.int=TRUE, # add confidence intervals
           pval=FALSE, # show the p-value for the log-rank test
           risk.table=FALSE, # show a risk table below the plot
           legend.labs=c("N", "Y"), # change group labels
           legend.title="ASSIGN > 20",  # add legend title
           palette=c("dodgerblue4", "orchid2"), # change colors of the groups
           title="Kaplan-Meier Curve for incident CVD", # add title to plot
           risk.table.height=.2)



mod1 <- survfit(Surv(tte, event) ~ factor(at_risk), data = cox)
ggsurvplot(mod1,
           conf.int=TRUE, # add confidence intervals
           pval=FALSE, # show the p-value for the log-rank test
           risk.table=FALSE, # show a risk table below the plot
           legend.labs=c("N", "Y"), # change group labels
           legend.title="ASSIGN > 20",  # add legend title
           palette=c("dodgerblue4", "orchid2"), # change colors of the groups
           title="Kaplan-Meier Curve for incident CVD", # add title to plot
           risk.table.height=.2)

mod1 = coxph(Surv(cox$tte, cox$event) ~ cox$age + factor(cox$sex))
summary(mod1)
residual_test = cox.zph(mod1)
residual_test

plot(mod1, fun = "cloglog", xlab = "Time (in days) using log",
     ylab = "log-log survival", main = "log-log curves by clinic")
