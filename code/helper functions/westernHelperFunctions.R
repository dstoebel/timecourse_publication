#Helper function for western blot analysis for western RpoS levels figure - taken from Dan Stoebel's repo westernAnalysis

library(tidyverse)

#This function reads the raw data for both the antibody of interest (usually RpoS), and values to normalize by (usually RpoD or revert stain)
#normalizes, and also records other critical information to set up a tidy data table:

normalizeSingleBlot <- function(normalizeValues, antibodyValues) {
  #order tables by lane
  
  #select normalize signal, keep two values, and rename one of them.
  
  normalizeData <- dplyr::select(normalizeValues, Signal, Lane) %>%
    dplyr::rename(normalize = Signal)

  
  #select western signal
  
  blotData <- dplyr::select(antibodyValues, -`Image Name`, -Channel, -Name, -Total, -Area, -Bkgnd., -Type,
                     -`Conc. Std.`, -Concentration) %>%
    dplyr::rename(RpoS = Signal) %>%
    left_join(normalizeData, by="Lane") %>%
    mutate(RpoS/normalize) %>%
    dplyr::rename(RpoSratio = `RpoS/normalize`)
  
  return(blotData)
  
}

