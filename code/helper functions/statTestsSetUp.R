#all data wrangling and set up needed for all statistical testing in`statstical_testing.Rmd`

library(tidyverse)
library(readxl)
library(infer)

#all paths assume that this file is being run from `statstical_testing.Rmd`. 
#If you need to run it from a file with a different path length, this will need to be changed.
source(file = "../helper functions/onsetTimesHelperFunctions.R")

sensitivity <- read_excel("../../data/other/Sensitivity.xlsx")
all_genes <- read_excel("../../data/other/zjb999094348sd2(1).xlsx") 

#read sicegar results from file
osmo_sicegar_list <-read.csv("../../outputs/sicegar/osmo_sicegar_data.csv") %>%
  mutate(Experiment = "Osmotic Shock")
all_cold_sicegar_list <-read.csv("../../outputs/sicegar/recount_cold_sicegar_data.csv") %>%
  dplyr::select(geneName, value)%>%
  mutate(Experiment = "Cold Shock")
all_SP_sicegar_list <-read.csv("../../outputs/sicegar/recount_SP_sicegar_data.csv") %>%
  dplyr::select(geneName, value)%>%
  mutate(Experiment = "Stationary Phase")

all_profile <- all_genes %>%
  filter(sensitivity != "nsRegulation")%>%
  filter(sensitivity != "linear")
sens <- all_profile[,c(2,8)] 

#use helper function to get only genes with sensitivity profiles & add column with that data
SPonTimes<-sensitivitySep(sensitivity, all_SP_sicegar_list)
OsmoonTimes<-sensitivitySep(sensitivity, osmo_sicegar_list)
AColdonTimes<-sensitivitySep(sensitivity, all_cold_sicegar_list)

#get list of all genes with sensitivity profiles for each stressor
Osmogenes <- semi_join(osmo_sicegar_list, all_profile, by = "geneName" )
Osmogenes <- left_join(Osmogenes, sens, by = "geneName") %>%
  filter(!is.na(value))

SPgenes <- semi_join(all_SP_sicegar_list, all_profile, by = "geneName" )
SPgenes <- left_join(SPgenes, sens, by = "geneName") %>%
  filter(!is.na(value))

AColdgenes <- semi_join(all_cold_sicegar_list, all_profile, by = "geneName" )
AColdgenes <- left_join(AColdgenes, sens, by = "geneName")%>%
  filter(!is.na(value))

#make tables that combines 2 stressors into 1 data set 
#change value column so join will work
Osmogenes_val <- Osmogenes %>% dplyr::rename(valOsmo = value)
SPgenes_val <- SPgenes %>% dplyr::rename(valSP = value)
Coldgenes_val <- AColdgenes %>% dplyr::rename(valCold = value)

#join
Osmo_Cold_genes <- left_join(Osmogenes_val, Coldgenes_val, by = "geneName") %>%
  dplyr::rename(sensitivity = sensitivity.x) %>%
  dplyr::select(-sensitivity.y)%>%
  filter(!is.na(valCold))
Osmo_SP_genes <- left_join(Osmogenes_val, SPgenes_val, by = "geneName")%>%
  dplyr::rename(sensitivity = sensitivity.x) %>%
  dplyr::select(-sensitivity.y)%>%
  filter(!is.na(valSP))
Cold_SP_genes <- left_join(Coldgenes_val, SPgenes_val, by = "geneName")%>%
  dplyr::rename(sensitivity = sensitivity.x) %>%
  dplyr::select(-sensitivity.y)%>%
  filter(!is.na(valSP))

#find onset times
osmo_onset <- OsmoonTimes %>%
  filter(!is.na(value))

#cold - there are only 3 insensitive genes(that aren't ambiguous) so we won't really be able to say much about that 
cold_onset <- AColdonTimes %>%
  filter(!is.na(value))

#SP
sp_onset <- SPonTimes %>%
  filter(!is.na(value))


