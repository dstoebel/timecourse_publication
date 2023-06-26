# Code, data, and outputs for *The timing of transcription of RpoS-dependent genes varies across multiple stresses* by Adams et al.





[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8083605.svg)](https://doi.org/10.5281/zenodo.8083605)




## Quick Start: 
1. code - contains the code for basic data wrangling, figures, other results & analysis such as statistical testing, and helper functions.
2. data - contains all the data files for data analysis, including raw data in htseq-count tables, western blot qunatifcation data, and files containing information of sensitivity 3. profiles.
outputs - contains the graphs made for figures, sicegar results, GO analysis results (although that is from PANTHER), and other data files we produced.
4. renv - contains necessary files for the R package `renv` to work so that the project is reproducible and all the needed libraries are downloaded.

## Specifics 

###Code directory: 
 - This is split into four folders that create all of the figures and the analysis for the paper.
 - `helper functions` contains data wrangling scripts used by almost all of the `.Rmd` files that do the analysis.
- Supplementary figures contains the code to create a supplemental figure plotting the onset timing of all DE genes. 

###Data Directory: 
This directory is split into three folders each of which contain a different kind of data that the paper uses. 
- RNAseq contains the counts files from each data set, as output by htseq-count. moc1533 holds all stationary phase counts (short and long - LSP), moc1691 contains high osmolarity (osmotic shock - OS), and mocp0022 has all low temperature counts (short and long). These files are called in the data_wrangling.R file in the code directory.
- westerns contains .xls files from the imaging and quantification of western blots. There are two for each western - the total REVERT stain quantification that is used to normalize RpoS and the RpoS quantification. Apart from the amount of fluorecence detected these files also usually contain lane, genotype, time point, stress, OD, and experiment number information. These files are called in westernAnalysisSetUp.R These files are broken into folders for each of the three stresses.
- other 
  - Contains all other data needed in different data anlysis or figures. Data in this file includes sensitivity info from wong et al, the gff with basic info and names of genes, and more. 
  - If a bunch of files of one type of data need to be uploaded they should get there own folder. What goes here are more miscellaneous         data that isn't easily group with others. 
  
###Outputs directory: 
This is split into three folders and some outlying .csv files. The .csv files are either individual DEseq result files, a combined DEseq results file with all information like ajusted p-vals and log2 fold change, or a file with the DESeq data and sicgear onset times combined.
- GO contains the .txt files containing just the GeneIDs of genes in a group of DE genes (such as DE in stationary phase). These files are all tab delimited to work with the PANTHER website. 
- graphs contains the figures as output by R for the paper, as well as annotated versions of the venn diagram.
- sicegar contains the onset times for genes in each environment. These files are output by R code. The same numbers are present ing `combined_deseq_sicegar.csv`
- Supplemental contains a few results that are supplemental figures for the paper.
	- Supp western figure.pdf shows example western blots.
	- PANTHER results holds .txt files of the outputs of the statistically significant functions for each geneset. 
	- sicegar contains the plots of each DE gene with time in each stress. The sicegar fit is shown if it exists. The `*_all_graphs.pdf` files show all of the plots for a given stress.

###renv folder and renv.lock
These are used by the R package `renv` for keeping track of versions of packages used in this analysis to support a reproducible analysis.
  
