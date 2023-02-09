#this file does all data wrangling and set up needed for western figure making and analysis

library(tidyverse)
library(readxl)
theme_set(theme_classic())


#The paths of all the files assume you are calling this file from the code/figures/fig1 directory
#It also works for the analysis file in the code/results/western analysis directory 
source("../helper functions/westernHelperFunctions.R") #This loads the helper functions for normalization, etc.

# pulling western data from excel files
#in the names of the read_excel data frames, the number refers to a specific day the blot was done, and the letter refers to the replicate. If there is no letter that means there was only one blot done that day. 

##osmo data 
#The osmotic shock western data is coming from Johnson Hoangs repo called westernAnalysis. These western from four cold shocks which were blotted on July 1, 2019. There is a biological replicate of trials A and B from June 19, 2019. However, this western blot is dirty and speckled and not as clean as the July 1 replicate and will not be used here. 

osmo_revert2a <- read_excel("../../data/westerns/osmo/revert OSa OSb 2019 07 01.xls")
osmo_rpoS2a <- read_excel("../../data/westerns/osmo/stain OSa OSb 2019 07 01.xls")

osmo_revert2b <- read_excel("../../data/westerns/osmo/revert OSc OSd 2019 07 01.xls")
osmo_rpoS2b <- read_excel("../../data/westerns/osmo/stain OSc OSd 2019 07 01.xls")


##cold data 
# The stationary phase western data is from 2 different repos. The short course data (to 150 mins) is from Johnson Hoang's repo called westernAnalysis. The long course data (300 mins) is from Josephine Adams' branch of Dan Stoebels repo called westernAnalysis. 
# There were three biological replicates of the third short western blot from June 2019. We chose the blot from June 18 because it was the cleanest and had the most clear lines out of any of the three blots. 

#short 
cold_revertS1a <- read_excel("../../data/westerns/cold/revert CSa CSb 2019 07 09.xls")
cold_rpoSS1a <- read_excel("../../data/westerns/cold/stain CSa CSb 2019 07 09.xls")

cold_revertS1b <- read_excel("../../data/westerns/cold/revert CSc CSd 2019 07 09.xls")
cold_rpoSS1b <- read_excel("../../data/westerns/cold/stain CSc CSd 2019 07 09.xls")

cold_revertS2 <- read_excel("../../data/westerns/stationary phase/revert sp3 cs1 2019 06 18.xls")%>%
  filter(Experiment == "cold shock")
cold_rpoSS2 <- read_excel("../../data/westerns/stationary phase/stain_SP3_CS1_2019-06-18.xls")%>%
  filter(Experiment == "cold shock")

#long 
cold_revertL1 <- read_excel("../../data/westerns/cold/ShapesTableExportRevert_2021-06-30.xls")
cold_rpoSL1 <- read_excel("../../data/westerns/cold/ShapesTableExportBlot_21_2021-06-30.xls")

cold_revertL2 <- read_excel("../../data/westerns/cold/ShapesTableExportrEVERT_7_22_2021-07-27.xls")
cold_rpoSL2 <- read_excel("../../data/westerns/cold/ShapesTableExportBlot_7_22_2021-07-27.xls")

cold_revertL3 <- read_excel("../../data/westerns/cold/ShapesTableExportBlot_0_3_REVERT.xls")
cold_rpoSL3 <- read_excel("../../data/westerns/cold/ShapesTableExportBlot_9_3_BLOT.xls")

cold_revertL4 <- read_excel("../../data/westerns/cold/ShapesTableExportBlot_9_17_revert.xls")
cold_rpoSL4 <- read_excel("../../data/westerns/cold/ShapesTableExportBlot_9_17_BLOT.xls")

cold_revertL5 <- read_excel("../../data/westerns/cold/ShapesTableExportRevert_2021-06-29.xls")
cold_rpoSL5 <- read_excel("../../data/westerns/cold/ShapesTableExportBlot_2021-06-29 copy.xls")

cold_revertL6 <- read_excel("../../data/westerns/cold/ShapesTableExport_2021-06-17.xls")
cold_rpoSL6 <- read_excel("../../data/westerns/cold/ShapesTableExportBlot_2021-06-17.xls")


##SP data 
# The stationary phase western data is from 2 different repos. The short course data (to 150 mins) is from Johnson Hoang's repo called westernAnalysis. The long course data (300 mins) is from his repo called johnson_longer_stationary. 
# There were three replicates of the short-course western blots from June 2019. In order to have just one set of data from each stationary phase trial, we picked the SP1 and SP2 western blot from June 12 and the SP3 Blot from June 18. This is because it was these were the cleanest and most clear of the three replicates for each sample. Part of the last lane from the June 12, which was for Δrpos is cut off. However, since Δrpos is control and its REVERT levels are not used in this analysis, we decided that this was the best blot to use.  

#short 
SP_revertS1a <- read_excel("../../data/westerns/stationary phase/revert sp1 sp2 2019 06 12.xls")
SP_rpoSS1a <- read_excel("../../data/westerns/stationary phase/stain_sp1_sp2_2019-06-12.xls")

SP_revertS1b <- read_excel("../../data/westerns/stationary phase/revert sp3 cs1 2019 06 18.xls")%>%
  filter(Experiment == "stationary phase") %>%
  mutate(Trial = "F")
SP_rpoSS1b <- read_excel("../../data/westerns/stationary phase/stain_SP3_CS1_2019-06-18.xls")%>%
  filter(Experiment == "stationary phase")%>%
  mutate(Trial = "F")

#long 
SP_revertL1a <- read_excel("../../data/westerns/stationary phase/revert LSPa2019 11 11.xlsx")
SP_rpoSL1a <- read_excel("../../data/westerns/stationary phase/western LSPa 2019 11 11.xlsx")

SP_revertL1b <- read_excel("../../data/westerns/stationary phase/revert LSPc 2019 11 11.xlsx")
SP_rpoSL1b <- read_excel("../../data/westerns/stationary phase/western LSPc 2019 11 11.xlsx")

SP_revertL1c <- read_excel("../../data/westerns/stationary phase/revert LSPd 2019 11 11.xls")
SP_rpoSL1c <- read_excel("../../data/westerns/stationary phase/western LSPd 2019 11 11.xls")


#calling the helper  function in order to get the normalized signal for every set of Western blots 

##osmo
osmo_blotData2a <- normalizeSingleBlot(osmo_revert2a, osmo_rpoS2a)
osmo_blotData2b <- normalizeSingleBlot(osmo_revert2b, osmo_rpoS2b)

#make into one table
osmo_allData <- bind_rows(osmo_blotData2a, osmo_blotData2b) %>%
  mutate(Experiment = "High Osmolarity")


##cold
cold_blotDataS1a <- normalizeSingleBlot(cold_revertS1a, cold_rpoSS1a) %>%
  mutate(exp_length = "short") # to keep track of which length of experiment data points are from
cold_blotDataS1b <- normalizeSingleBlot(cold_revertS1b, cold_rpoSS1b)%>%
  mutate(exp_length = "short")
cold_blotDataS2 <- normalizeSingleBlot(cold_revertS2, cold_rpoSS2) %>%
  mutate(exp_length = "short") 

cold_blotDataL1 <- normalizeSingleBlot(cold_revertL1, cold_rpoSL1)%>%
  mutate(exp_length = "long")
cold_blotDataL2 <- normalizeSingleBlot(cold_revertL2, cold_rpoSL2)%>%
  mutate(exp_length = "long")
cold_blotDataL3 <- normalizeSingleBlot(cold_revertL3, cold_rpoSL3)%>%
  mutate(exp_length = "long")
cold_blotDataL4 <- normalizeSingleBlot(cold_revertL4, cold_rpoSL4)%>%
  mutate(exp_length = "long")
cold_blotDataL5 <- normalizeSingleBlot(cold_revertL5, cold_rpoSL5)%>%
  mutate(exp_length = "long") %>%
  mutate(Trial = "I")
cold_blotDataL6 <- normalizeSingleBlot(cold_revertL6, cold_rpoSL6)%>%
  mutate(exp_length = "long")

#make into one table
cold_allData <- bind_rows(cold_blotDataS1a, cold_blotDataS1b, cold_blotDataS2, cold_blotDataL1, cold_blotDataL2, cold_blotDataL3, cold_blotDataL4, cold_blotDataL5, cold_blotDataL6) %>%
  mutate(Experiment = "Low Temperature")


##SP
SP_blotDataS1a <- normalizeSingleBlot(SP_revertS1a, SP_rpoSS1a)%>%
  mutate(exp_length = "short") # to keep track of which length of experiment data points are from
SP_blotDataS1b <- normalizeSingleBlot(SP_revertS1b, SP_rpoSS1b)%>%
  mutate(exp_length = "short") 

SP_blotDataL1a <- normalizeSingleBlot(SP_revertL1a, SP_rpoSL1a)%>%
  mutate(exp_length = "long") %>%
  mutate(Trial = "C")
SP_blotDataL1b <- normalizeSingleBlot(SP_revertL1b, SP_rpoSL1b)%>%
  mutate(exp_length = "long") %>%
  mutate(Trial = "D")
SP_blotDataL1c <- normalizeSingleBlot(SP_revertL1c, SP_rpoSL1c)%>%
  mutate(exp_length = "long") %>%
  mutate(Trial = "E")

#make into one table
SP_allData <- bind_rows(SP_blotDataS1a, SP_blotDataS1b, SP_blotDataL1a, SP_blotDataL1b, SP_blotDataL1c)%>%
  mutate(Experiment = "Stationary Phase")

#make one table with all data 
allData <- merge(SP_allData,cold_allData, all = TRUE)
allData <- merge(allData, osmo_allData, all = TRUE) 
