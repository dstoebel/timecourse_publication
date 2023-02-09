#Helper function and some data wrangling for looking at onset times from sicegar

library(readxl)
library(tidyverse)
library(sicegar)
##set up for sensitivity analysis 

#getting sicegar results from file
osmo_sicegar_list <-read.csv("../../outputs/sicegar/osmo_sicegar_data.csv") %>%
  select(geneName, value) %>% 
  mutate(Experiment = "High Osmolarity")
all_cold_sicegar_list <-read.csv("../../outputs/sicegar/recount_cold_sicegar_data.csv") %>%
  select(geneName, value)%>%
  mutate(Experiment = "Low Temperature")
all_SP_sicegar_list <-read.csv("../../outputs/sicegar/recount_SP_sicegar_data.csv") %>%
  select(geneName, value)%>%
  mutate(Experiment = "Stationary Phase")

#combine into one data fram
all_sicegar_list <- merge(osmo_sicegar_list, all_cold_sicegar_list, all = TRUE)
all_sicegar_list <- merge(all_sicegar_list, all_SP_sicegar_list, all = TRUE)

#all genes & profiles (sens, insens, linear, etc) from wong et al paper
all_genes <- read_excel("../../data/other/zjb999094348sd2(1).xlsx") 

#Getting sensitivity profiles from wong et al 2017 paper
sensitivity <- read_excel("../../data/other/Sensitivity.xlsx")



#helper function for onset time figure - taken from Dan Stoebel's repo emily_summer_2020
#takes in list of genes and data frame containing genes identified as sensitive/insensitive
#outputs the same dataframe with new column for sensitivity profiles
sensitivitySep<- function(sensitivitydf, tdf){
  sensitivitydf<-split(sensitivitydf, sensitivitydf$sensitivity)
  sensitive<-sensitivitydf$sensitive
  linear<-sensitivitydf$linear
  insensitive<-sensitivitydf$insensitive
  
  
  sensList<- tdf %>% filter(geneName %in% sensitive$geneName)
  
  insensList<- tdf %>% filter(geneName %in% insensitive$geneName)
  
  linearList <-tdf %>% filter(geneName %in% linear$geneName)
  
  #creates type column, combines sens and insens into one df
  sensList$type<-rep("Sensitive", dim(sensList)[1])
  insensList$type<-rep("Insensitive", dim(insensList)[1])
  linearList$type<-rep("Linear", dim(linearList)[1])
  on_times<-rbind(linearList, sensList, insensList)
  on_times
}


##From: emily_summer_2020/helper functions/sicegar_map for graphing sicegar results
##to be used with the "map" tidyr function for running Sicegar and generating graphs for every unit of observation of a dataset
##input: dat=data corresponding to a single gene (tidy data frame), AIC=threshold AIC score (default is -10), thresholdRatio=threshold intensity ratio (default is 0.75)
##output: a graph of the fitted curve with onset time labeled
sicegarGraphMap <- function(dat, type="", AIC=-10, thresholdRatio=0.75){
  gene=dat$geneName[1]
  dat2<-dat
  dat<-dat %>% select(intensity,time)
  dat <- dat %>% arrange(time) #arranges by increasing time for reversal if counts are decreasing
  #model fitting and categorization
  Model <- fitAndCategorize(dataInput=dat, threshold_t0_max_int=1E10, threshold_dsm_tmax_IntensityRatio = thresholdRatio, threshold_AIC=AIC)
  #if ambiguous, check to see if decreasing (avg intensity at end is smaller than beginning), if so, then flip the counts so they are increasing and fit the model again
  if(Model$summaryVector$decision=="ambiguous"){
    if(mean(dat$intensity[dat$time==0])>mean(dat$intensity[dat$time==max(dat$time)])){
      dat$intensity=rev(dat$intensity) #reverses counts
      Model <- fitAndCategorize(dataInput=dat, threshold_t0_max_int=1E10, threshold_dsm_tmax_IntensityRatio = thresholdRatio, threshold_AIC=AIC)
      #if sigmoidal, graph the sigmoidal curve
      if(Model$summaryVector$decision=="sigmoidal"){
        a<-figureModelCurves(dataInput=Model$normalizedInput, sigmoidalFitVector = Model$sigmoidalModel, showParameterRelatedLines = TRUE)+ggtitle(paste(gene, "Sigmoidal:", type, sep=" "))
      }
      #if sigmoidal, graph the double sigmoidal curve
      if(Model$summaryVector$decision=="double_sigmoidal"){
        a<-figureModelCurves(dataInput=Model$normalizedInput, doubleSigmoidalFitVector = Model$doubleSigmoidalModel, showParameterRelatedLines = TRUE)+ggtitle(paste(gene, "Impulse:", type, sep=" "))
      }
      #if ambiguous, just graph points with no curve (color corresponds to short or long time course)
      if(Model$summaryVector$decision=="ambiguous"){
        a<-ggplot(dat2, aes(x=time, y=intensity, color=replicate)) + geom_point() + ggtitle(paste(gene, "(ambiguous)", sep=" "))
      }}
    else {
      a<-ggplot(dat2, aes(x=time, y=intensity, color=replicate)) + geom_point() + ggtitle(paste(gene, "(ambiguous)", sep=" "))
    }
  }
  if(Model$summaryVector$decision=="sigmoidal"){
    a<-figureModelCurves(dataInput=Model$normalizedInput, sigmoidalFitVector = Model$sigmoidalModel, showParameterRelatedLines = TRUE)+ggtitle(paste(gene, "Sigmoidal:", type, sep=" "))
  }
  if(Model$summaryVector$decision=="double_sigmoidal"){
    a<-figureModelCurves(dataInput=Model$normalizedInput, doubleSigmoidalFitVector = Model$doubleSigmoidalModel, showParameterRelatedLines = TRUE)+ggtitle(paste(gene, "Impulse:", type, sep=" "))
  }
  a
}