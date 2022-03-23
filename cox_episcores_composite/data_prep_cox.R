library(dplyr)

box::use(../../modules/cox[...])

####################################### Input data #######################################

agesex = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/2022-03-10_age_alive.csv');
agemonths = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/agemonths.csv');
cvd_deaths = read.csv('/Volumes/marioni-lab/Ola/Lab/Cox/2022-01-18/cvd_deaths_withcode.csv'); # 696
non_cvd_deaths = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/2021-09-27_death.csv');
hosp_events = read.csv('/Volumes/marioni-lab/Ola/Lab/Cox/2022-01-18/cvd_hospitalisations_and_operations_withcode.csv') # 3666
assign =  read.csv('/Volumes/marioni-lab/Ola/Lab/ASSIGN/GS_assign_data.csv')
target = read.csv('/Volumes/marioni-lab/Ola/Lab/Test_sets/gs20ktargets.tsv', sep='\t')

####################################### CVD subsets #######################################

events = hosp_events
events = subset(events, Incident==1 & Event_Type != "GP" & Event_Type != "Death") # 1297

main_ds = assign %>%
  left_join(agesex, c("id" = "id")) %>%
  left_join(target, c("id" = "Sample_Name")) %>%
  left_join(non_cvd_deaths, c("id" = "id")) %>%
  left_join(cvd_deaths, c("id" = "id"))

main_ds = main_ds[which(!is.na(main_ds$age.y)),]

ds = main_ds %>% left_join(events, c("id" = "id"))

ds = ds[c("id", "ID", "sex.x", "age.y", "assign", "Set", "yob", "mob", "dod_ym", 
            "Incident.x", "Event_Type.x", "Event_Date.x", "Incident.y", "Event_Type.y", "Event_Date.y")] 
colnames(ds) = c("id", "Sentrix_ID", "sex", "age", "assign", "set", "yob", "mob", "dod_ym", 
                 "Incident.Death", "Event_Type.Death", "Event_Date.Death", "Incident.Event", "Event_Type.Event", "Event_Date.Event")


####################################### Cox dataset #######################################

cox = make_cox_dataset(ds)
cox = cox[c("id", "Sentrix_ID", "sex", "age", "assign", "set", "dead", "event", "tte")] 
cox = subset(cox, !is.na(tte) & tte>0)
write.csv(cox, '/Volumes/marioni-lab/Ola/Lab/EpiScores/Cox_episcores_composite/cox_covars.csv', row.names = F)
