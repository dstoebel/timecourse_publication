#This file does all the necessary set up for DEseq and sicegar as well as running those. Sicegar is commented out because it can take a while to run and the outputs have been saved to file so you don't have to run it.
#The paths in this file are setup to work when its run from `combined_DESeq_analysis.rmd` and thus currently won't work if you just run the file directly. 
#setup
knitr::opts_chunk$set(echo = TRUE)

library(DESeq2)
library(tidyverse)
library(sicegar)
library(gdata)
library(utils)
library(limma)
theme_set(theme_classic())


#Read in the gff file to get basic information about the genes, then format it using `dplyr` functions
gff <- read_delim("../../data/other/GCF_000005845.2_ASM584v2_genomic.gff", 
                  "\t", escape_double = FALSE, col_names = FALSE, 
                  trim_ws = TRUE, skip = 7)

gene_info <- gff %>%
  filter(X3 == "gene") %>%
  dplyr::rename(identifiers = X9) %>%
  extract(identifiers, into = c("bNum", "GeneID" ,"name"), 
          regex="ID=gene-([[:alnum:]]+);Dbxref=ASAP.ABE-[[:digit:]]+,ECOCYC:[[:alnum:]]+,(GeneID:[[:digit:]]+);Name=([[:alnum:]]+)",
          remove = FALSE) %>%
  dplyr::select(bNum,  name, GeneID) %>%
  filter(!is.na(bNum))

###stationary phase data wrangling

##Set up for DESeq
#importing data

#moved long SP data from moc1816 to moc1533 so that all SP data is in one place

##needs to add long data
SPFiles <- c( "JH01_A01_counts.txt", "JH01_A02_counts.txt", "JH01_A03_counts.txt", "JH01_A04_counts.txt", "JH01_A05_counts.txt","JH01_A06_counts.txt","JH01_B01_counts.txt", "JH01_B02_counts.txt", "JH01_B03_counts.txt", "JH01_B04_counts.txt", "JH01_B05_counts.txt","JH01_B06_counts.txt", "JH01_C01_counts.txt", "JH01_C02_counts.txt",  "JH01_C03_counts.txt", "JH01_C04_counts.txt", "JH01_C05_counts.txt","JH01_C06_counts.txt",
              "JH01_LSP_A_01_counts.txt","JH01_LSP_A_02_counts.txt","JH01_LSP_A_03_counts.txt","JH01_LSP_A_04_counts.txt","JH01_LSP_A_05_counts.txt","JH01_LSP_A_06_counts.txt","JH01_LSP_C_01_counts.txt","JH01_LSP_C_02_counts.txt","JH01_LSP_C_03_counts.txt","JH01_LSP_C_04_counts.txt","JH01_LSP_C_05_counts.txt","JH01_LSP_C_06_counts.txt","JH01_LSP_E_01_counts.txt","JH01_LSP_E_02_counts.txt","JH01_LSP_E_03_counts.txt","JH01_LSP_E_04_counts.txt","JH01_LSP_E_05_counts.txt","JH01_LSP_E_06_counts.txt",
              "JH02_A01_counts.txt", "JH02_A02_counts.txt", "JH02_A03_counts.txt", "JH02_A04_counts.txt", "JH02_A05_counts.txt","JH02_A06_counts.txt","JH02_B01_counts.txt", "JH02_B02_counts.txt", "JH02_B03_counts.txt", "JH02_B04_counts.txt", "JH02_B05_counts.txt","JH02_B06_counts.txt", "JH02_C01_counts.txt", "JH02_C02_counts.txt",  "JH02_C03_counts.txt", "JH02_C04_counts.txt", "JH02_C05_counts.txt","JH02_C06_counts.txt",
              "JH02_LSP_A_01_counts.txt","JH02_LSP_A_02_counts.txt","JH02_LSP_A_03_counts.txt","JH02_LSP_A_04_counts.txt","JH02_LSP_A_05_counts.txt","JH02_LSP_A_06_counts.txt","JH02_LSP_C_01_counts.txt","JH02_LSP_C_02_counts.txt","JH02_LSP_C_03_counts.txt","JH02_LSP_C_04_counts.txt","JH02_LSP_C_05_counts.txt","JH02_LSP_C_06_counts.txt","JH02_LSP_E_01_counts.txt","JH02_LSP_E_02_counts.txt","JH02_LSP_E_03_counts.txt","JH02_LSP_E_04_counts.txt","JH02_LSP_E_05_counts.txt","JH02_LSP_E_06_counts.txt")

SPNames <-c("wt_SP_A_0_short", "wt_SP_A_30_short", "wt_SP_A_60_short", "wt_SP_A_90_short", "wt_SP_A_120_short", "wt_SP_A_150_short","wt_SP_B_0_short", "wt_SP_B_30_short", "wt_SP_B_60_short", "wt_SP_B_90_short", "wt_SP_B_120_short", "wt_SP_B_150_short","wt_SP_C_0_short", "wt_SP_C_30_short", "wt_SP_C_60_short", "wt_SP_C_90_short", "wt_SP_C_120_short", "wt_SP_C_150_short",
            "wt_SP_D_0_long", "wt_SP_D_60_long", "wt_SP_D_120_long", "wt_SP_D_180_long", "wt_SP_D_240_long", "wt_SP_D_300_long", "wt_SP_E_0_long", "wt_SP_E_60_long", "wt_SP_E_120_long", "wt_SP_E_180_long", "wt_SP_E_240_long", "wt_SP_E_300_long","wt_SP_F_0_long", "wt_SP_F_60_long", "wt_SP_F_120_long", "wt_SP_F_180_long", "wt_SP_F_240_long", "wt_SP_F_300_long",
            "del_SP_A_0_short", "del_SP_A_30_short", "del_SP_A_60_short", "del_SP_A_90_short", "del_SP_A_120_short", "del_SP_A_150_short", "del_SP_B_0_short", "del_SP_B_30_short", "del_SP_B_60_short", "del_SP_B_90_short", "del_SP_B_120_short", "del_SP_B_150_short","del_SP_C_0_short", "del_SP_C_30_short", "del_SP_C_60_short", "del_SP_C_90_short", "del_SP_C_120_short", "del_SP_C_150_short",
            "del_SP_D_0_long", "del_SP_D_60_long", "del_SP_D_120_long", "del_SP_D_180_long", "del_SP_D_240_long", "del_SP_D_300_long", "del_SP_E_0_long", "del_SP_E_60_long", "del_SP_E_120_long", "del_SP_E_180_long", "del_SP_E_240_long", "del_SP_E_300_long","del_SP_F_0_long", "del_SP_F_60_long", "del_SP_F_120_long", "del_SP_F_180_long", "del_SP_F_240_long", "del_SP_F_300_long")

SPTime <- c(rep(c(rep(c("0","30","60", "90", "120", "150"),3), rep(c("0", "90", "120", "180", "240", "300"),3)),2))

SPGenotype <- c(rep("wt", 36), rep("del", 36))

SPBatch <- c(rep(c(rep("short",18), rep("long", 18)),2))

SPSampleTable <- data.frame(sampleNames = SPNames,
                            fileNames = SPFiles,
                            genotype = SPGenotype,
                            time  = SPTime,
                            batch = SPBatch)

SPDESeqTable <- DESeqDataSetFromHTSeqCount(sampleTable  = SPSampleTable,
                                           directory = "../../data/RNAseq/moc1533",
                                           design =  ~  batch + time + genotype + time:genotype)

##running DESeq & getting results
ecoliddsAllSP<-DESeq(SPDESeqTable, test="LRT", reduced= ~batch + time + genotype)

#get significant genes
SP_results <-  results(ecoliddsAllSP, alpha = 0.01)

SP_sig_genes_no_gene <- SP_results %>%
  as.data.frame() %>%
  as_tibble(rownames = "bNum")

#fixing coverging in beta issue when running DESeq
#using the command `results(ecoliddsAllSP)[!mcols(ecoliddsAllSP)$fullBetaConv,]` we see that there are 8 genes that never converge using during the LRT test. Since there are only 8 genes with very few counts and, thus, very little power, we will remove these genes from the data set
bad_bNums = c("b0790", "b2302", "b3513", "b4149", "b4433", "b4435", "b4704", "b4755")

SP_sig_genes_no_gene <- SP_sig_genes_no_gene %>%
  filter(!bNum %in% bad_bNums)

SP_sig_genes <- full_join(gene_info,SP_sig_genes_no_gene, by = "bNum")%>%
  filter(!is.na(padj))

normecounts_SP <- counts(ecoliddsAllSP, normalized=TRUE) %>%
  as.data.frame() %>%
  as_tibble(rownames = "bNum") %>%
  left_join(SP_sig_genes, by = "bNum")

tidy_normecounts_SP <- normecounts_SP %>%
  filter(padj < 0.01) %>%
  select(-baseMean, -log2FoldChange, -lfcSE, -stat, -padj, -pvalue) %>%
  pivot_longer(cols = contains("wt_"), names_to = "sample_name", values_to = "intensity") %>%
  tidyr::separate(sample_name, into = c("replicate","stress", "genotype", "time"), sep = "_") %>%
  mutate(time = as.numeric(time))

#changing column names to match previous
tidy_normecounts_names_SP<-tidy_normecounts_SP %>%
  dplyr::rename(geneName = name)

##using `vst` & `removeBatchEffects`
#vst transformation
vsd_deseq_SP <- vst(ecoliddsAllSP , blind=FALSE)

mat_deseq_SP <- assay(vsd_deseq_SP)
mm_deseq_SP <- model.matrix(~time+genotype + time:genotype, colData(vsd_deseq_SP))
no_batch_matrix_deseq_SP <- limma::removeBatchEffect(mat_deseq_SP, batch=vsd_deseq_SP$batch, design=mm_deseq_SP)
assay(vsd_deseq_SP) <- no_batch_matrix_deseq_SP
no_batch_tibble_SP <- no_batch_matrix_deseq_SP %>%
  as.data.frame() %>%
  as_tibble(rownames = "gene")%>%
  dplyr::rename(bNum = gene)

no_batch_tibble_SP <- full_join(gene_info,no_batch_tibble_SP, by = "bNum")

vst_sig_genes_SP <- semi_join(no_batch_tibble_SP, SP_sig_genes, by = "name" )%>%
  pivot_longer(cols = contains("wt_"), names_to = "sample_name", values_to = "intensity") %>%
  tidyr::separate(sample_name, into = c("genotype","stress", "replicate","time", "batch"), sep = "_") %>%
  mutate(time = as.numeric(time)) %>%
  dplyr::rename(geneName = name)

vst_sig_genes_SP <- semi_join(vst_sig_genes_SP, tidy_normecounts_names_SP, by  = "geneName")

#reverse transform (treat vst as a log2 transform)
vst_sig_genes_SP <- vst_sig_genes_SP %>% mutate(reverseIntensity = 2 ^ intensity)

#select columns of interest
vst_sig_genes_SP <- vst_sig_genes_SP %>% select(-intensity) %>% dplyr::rename(intensity = reverseIntensity) %>%
  select(bNum, geneName, replicate, genotype,time, intensity, batch) %>%
  filter(!is.na(bNum))

#save this to file so it can be used - commented because it only needed to be run once
write_csv( vst_sig_genes_SP, file="../../outputs/DESeq_results_SP.csv")


###osmo data wrangling 

##importing data
#file names
osmoFiles <- c("JH01_OS_A_01_counts.txt","JH01_OS_A_02_counts.txt","JH01_OS_A_03_counts.txt","JH01_OS_A_04_counts.txt","JH01_OS_A_05_counts.txt","JH01_OS_A_06_counts.txt","JH01_OS_C_01_counts.txt","JH01_OS_C_02_counts.txt","JH01_OS_C_03_counts.txt","JH01_OS_C_04_counts.txt","JH01_OS_C_05_counts.txt","JH01_OS_C_06_counts.txt","JH01_OS_D_01_counts.txt","JH01_OS_D_02_counts.txt","JH01_OS_D_03_counts.txt","JH01_OS_D_04_counts.txt","JH01_OS_D_05_counts.txt","JH01_OS_D_06_counts.txt","JH02_OS_A_01_counts.txt","JH02_OS_A_02_counts.txt","JH02_OS_A_03_counts.txt","JH02_OS_A_04_counts.txt","JH02_OS_A_05_counts.txt","JH02_OS_A_06_counts.txt","JH02_OS_C_01_counts.txt","JH02_OS_C_02_counts.txt","JH02_OS_C_03_counts.txt","JH02_OS_C_04_counts.txt","JH02_OS_C_05_counts.txt","JH02_OS_C_06_counts.txt","JH02_OS_D_01_counts.txt","JH02_OS_D_02_counts.txt","JH02_OS_D_03_counts.txt","JH02_OS_D_04_counts.txt","JH02_OS_D_05_counts.txt","JH02_OS_D_06_counts.txt")

#replicate names
osmoNames <-c("wt_OS_A_0", "wt_OS_A_30","wt_OS_A_60", "wt_OS_A_90","wt_OS_A_120","wt_OS_A_150","wt_OS_C_0", "wt_OS_C_30","wt_OS_C_60", "wt_OS_C_90","wt_OS_C_120","wt_OS_C_150","wt_OS_D_0", "wt_OS_D_30","wt_OS_D_60", "wt_OS_D_90","wt_OS_D_120","wt_OS_D_150","del_OS_A_0", "del_OS_A_30","del_OS_A_60", "del_OS_A_90","del_OS_A_120","del_OS_A_150","del_OS_C_0", "del_OS_C_30","del_OS_C_60", "del_OS_C_90","del_OS_C_120","del_OS_C_150","del_OS_D_0", "del_OS_D_30","del_OS_D_60", "del_OS_D_90","del_OS_D_120","del_OS_D_150")

#time points for sample collection
osmoTime <- c(rep(c("0","30","60", "90", "120", "150"),6))

osmoGenotype <- c(rep(c(rep("wt", 18), rep("del", 18))))

#make all osmo info into one table
osmoSampleTable <- data.frame(sampleNames = osmoNames,
                              fileNames = osmoFiles, 
                              genotype = osmoGenotype, 
                              time  = osmoTime)

#set up DEseq table
osmoDESeqTable <- DESeqDataSetFromHTSeqCount(sampleTable  = osmoSampleTable, 
                                             directory = "../../data/RNAseq/moc1691",
                                             design =  ~  time + genotype + time:genotype)

##running DESeq & getting results 

ecoliddsOsmo<-DESeq(osmoDESeqTable, test="LRT", reduced= ~ time+genotype) 

#get significant genes 
osmo_results <-  results(ecoliddsOsmo, alpha = 0.01)


osmo_sig_genes_no_gene <- osmo_results %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "bNum")

osmo_sig_genes <- full_join(gene_info,osmo_sig_genes_no_gene, by = "bNum")


##Get the normalized counts of expression at each time point, and add the p-values.
normecounts_osmo <- counts(ecoliddsOsmo, normalized=TRUE) %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "bNum") %>% 
  left_join(osmo_sig_genes, by = "bNum")

##Convert to tidy format, create columns for time, genotype, and replicate, and save only wt data
tidy_normecounts_osmo <- normecounts_osmo %>% 
  filter(padj < 0.01) %>%  #picking the DE ones using our threshold 
  select(-baseMean, -log2FoldChange, -lfcSE, -stat, -padj, -pvalue) %>% 
  pivot_longer(cols = contains("wt_"), names_to = "sample_name", values_to = "intensity") %>% 
  tidyr::separate(sample_name, into = c("replicate","stress", "genotype", "time"), sep = "_") %>% 
  mutate(time = as.numeric(time)) %>%
  select(bNum, name, replicate, stress, genotype,time, intensity)

#changing column names to match previous 
tidy_normecounts_names_osmo<-tidy_normecounts_osmo %>%
  dplyr::rename(geneName = name)

#save this to file so it can be used 
write_csv( tidy_normecounts_names_osmo, file="../../outputs/DESeq_results_osmo.csv")

###cold data wrangling
##Set up for DESeq
#importing data

#moved short CS data from moc1691 to mocp0022 so that all cold shock data is in one place (so only stuff left in moc1691 is osmo data)

coldFiles <- c("A_del_0_counts.txt", "A_del_90_counts.txt","A_del_120_counts.txt","A_del_180_counts.txt","A_del_240_counts.txt","A_del_300_counts.txt","A_wt_0_counts.txt", "A_wt_90_counts.txt", "A_wt_120_counts.txt" ,"A_wt_180_counts.txt", "A_wt_240_counts.txt","A_wt_300_counts.txt","B_del_0_counts.txt", "B_del_90_counts.txt", "B_del_120_counts.txt", "B_del_180_counts.txt","B_del_240_counts.txt", "B_del_300_counts.txt", "B_wt_0_counts.txt", "B_wt_90_counts.txt", "B_wt_120_counts.txt" ,"B_wt_180_counts.txt", "B_wt_240_counts.txt","B_wt_300_counts.txt", "C_del_0_counts.txt", "C_del_90_counts.txt", "C_del_120_counts.txt","C_del_180_counts.txt","C_del_240_counts.txt","C_del_300_counts.txt","C_wt_0_counts.txt", "C_wt_90_counts.txt", "C_wt_120_counts.txt" , "C_wt_180_counts.txt", "C_wt_240_counts.txt", "C_wt_300_counts.txt", "D_del_0_counts.txt", "D_del_90_counts.txt", "D_del_120_counts.txt", "D_del_180_counts.txt", "D_del_240_counts.txt","D_del_300_counts.txt", "D_wt_0_counts.txt", "D_wt_90_counts.txt", "D_wt_120_counts.txt" ,"D_wt_180_counts.txt", "D_wt_240_counts.txt", "D_wt_300_counts.txt", "JH01_CS_1_01_counts.txt", "JH01_CS_1_02_counts.txt", "JH01_CS_1_03_counts.txt", "JH01_CS_1_04_counts.txt","JH01_CS_1_05_counts.txt","JH01_CS_1_06_counts.txt","JH01_CS_A_01_counts.txt","JH01_CS_A_02_counts.txt","JH01_CS_A_03_counts.txt","JH01_CS_A_04_counts.txt","JH01_CS_A_05_counts.txt","JH01_CS_A_06_counts.txt","JH01_CS_B_01_counts.txt","JH01_CS_B_02_counts.txt","JH01_CS_B_03_counts.txt","JH01_CS_B_04_counts.txt","JH01_CS_B_05_counts.txt","JH01_CS_B_06_counts.txt","JH02_CS_1_01_counts.txt","JH02_CS_1_02_counts.txt","JH02_CS_1_03_counts.txt","JH02_CS_1_04_counts.txt","JH02_CS_1_05_counts.txt","JH02_CS_1_06_counts.txt","JH02_CS_A_01_counts.txt","JH02_CS_A_02_counts.txt","JH02_CS_A_03_counts.txt","JH02_CS_A_04_counts.txt","JH02_CS_A_05_counts.txt","JH02_CS_A_06_counts.txt","JH02_CS_B_01_counts.txt","JH02_CS_B_02_counts.txt","JH02_CS_B_03_counts.txt","JH02_CS_B_04_counts.txt","JH02_CS_B_05_counts.txt","JH02_CS_B_06_counts.txt")

coldNames <-c("del_CS_A_0_long", "del_CS_A_90_long", "del_CS_A_120_long", "del_CS_A_180_long", "del_CS_A_240_long", "del_CS_A_300_long","wt_CS_A_0_long", "wt_CS_A_90_long","wt_CS_A_120_long", "wt_CS_A_180_long","wt_CS_A_240_long","wt_CS_A_300_long","del_CS_B_0_long", "del_CS_B_90_long", "del_CS_B_120_long", "del_CS_B_180_long", "del_CS_B_240_long", "del_CS_B_300_long","wt_CS_B_0_long", "wt_CS_B_90_long","wt_CS_B_120_long", "wt_CS_B_180_long", "wt_CS_B_240_long", "wt_CS_B_300_long", "del_CS_C_0_long", "del_CS_C_90_long", "del_CS_C_120_long", "del_CS_C_180_long", "del_CS_C_240_long", "del_CS_C_300_long","wt_CS_C_0_long", "wt_CS_C_90_long","wt_CS_C_120_long", "wt_CS_C_180_long", "wt_CS_C_240_long","wt_CS_C_300_long","del_CS_D_0_long", "del_CS_D_90_long", "del_CS_D_120_long", "del_CS_D_180_long", "del_CS_D_240_long", "del_CS_D_300_long","wt_CS_D_0_long", "wt_CS_D_90_long", "wt_CS_D_120_long", "wt_CS_D_180_long", "wt_CS_D_240_long","wt_CS_D_300_long","wt_CS_E_0_short", "wt_CS_E_30_short","wt_CS_E_60_short", "wt_CS_E_90_short","wt_CS_E_120_short","wt_CS_E_150_short", "wt_CS_F_0_short", "wt_CS_F_30_short","wt_CS_F_60_short", "wt_CS_F_90_short", "wt_CS_F_120_short", "wt_CS_F_150_short","wt_CS_G_0_short", "wt_CS_G_30_short","wt_CS_G_60_short", "wt_CS_G_90_short", "wt_CS_G_120_short","wt_CS_G_150_short","del_CS_E_0_short", "del_CS_E_30_short","del_CS_E_60_short", "del_CS_E_90_short","del_CS_E_120_short", "del_CS_E_150_short","del_CS_F_0_short", "del_CS_F_30_short", "del_CS_F_60_short", "del_CS_F_90_short","del_CS_F_120_short","del_CS_F_150_short","del_CS_G_0_short", "del_CS_G_30_short","del_CS_G_60_short", "del_CS_G_90_short","del_CS_G_120_short","del_CS_G_150_short")

coldTime <- as.factor(c(rep(c("0", "90", "120", "180", "240", "300"),8),rep(c("0","30","60", "90", "120", "150"),6)))

coldGenotype <- c(rep(c(rep("del",6), rep("wt",6)),4),rep(c(rep("wt", 18), rep("del", 18))))

coldBatch <- c(rep("long", 48), rep("short", 36))

coldSampleTable <- data.frame(sampleNames = coldNames,
                              fileNames = coldFiles,
                              genotype = coldGenotype,
                              time  = coldTime,
                              batch = coldBatch)

coldDESeqTable <- DESeqDataSetFromHTSeqCount(sampleTable  = coldSampleTable,
                                             directory = "../../data/RNAseq/mocp0022",
                                             design =  ~  batch + time + genotype + time:genotype)

##running DESeq & getting results

ecoliddsAllCold<-DESeq(coldDESeqTable, test="LRT", reduced= ~ batch+time+genotype)

#get significant genes
cold_results <-  results(ecoliddsAllCold, alpha = 0.01)

cold_sig_genes_no_gene <- cold_results %>%
  as.data.frame() %>%
  as_tibble(rownames = "bNum")

cold_sig_genes <- full_join(gene_info,cold_sig_genes_no_gene, by = "bNum")

normecounts_cold <- counts(ecoliddsAllCold, normalized=TRUE) %>%
  as.data.frame() %>%
  as_tibble(rownames = "bNum") %>%
  left_join(cold_sig_genes, by = "bNum")

#filter for significant ajusted p values & make the data set tidy
tidy_normecounts_cold <- normecounts_cold %>%
  filter(padj < 0.01) %>%
  select(-baseMean, -log2FoldChange, -lfcSE, -stat, -padj, -pvalue) %>%
  pivot_longer(cols = contains("wt_"), names_to = "sample_name", values_to = "intensity") %>%
  tidyr::separate(sample_name, into = c("replicate","stress", "genotype", "time"), sep = "_") %>%
  mutate(time = as.numeric(time))


#changing column names to match previous
tidy_normecounts_names_cold<-tidy_normecounts_cold %>%
  dplyr::rename(geneName = name)

##using `vst` & `removeBatchEffects`
#vst transformation
vsd_deseq_cold <- vst(ecoliddsAllCold , blind=FALSE)

mat_deseq_cold <- assay(vsd_deseq_cold)
mm_deseq_cold <- model.matrix(~time+genotype + time:genotype, colData(vsd_deseq_cold))
no_batch_matrix_deseq_cold <- limma::removeBatchEffect(mat_deseq_cold, batch=vsd_deseq_cold$batch, design=mm_deseq_cold)
assay(vsd_deseq_cold) <- no_batch_matrix_deseq_cold
no_batch_tibble_cold <- no_batch_matrix_deseq_cold %>%
  as.data.frame() %>%
  as_tibble(rownames = "gene")%>%
  dplyr::rename(bNum = gene)

no_batch_tibble_cold <- full_join(gene_info,no_batch_tibble_cold, by = "bNum")

vst_sig_genes_cold <- semi_join(no_batch_tibble_cold, cold_sig_genes, by = "bNum" )%>%
  pivot_longer(cols = contains("wt_"), names_to = "sample_name", values_to = "intensity") %>%
  tidyr::separate(sample_name, into = c( "genotype","stress", "replicate","time", "batch"), sep = "_")  %>%
  dplyr::rename(geneName = name)

vst_sig_genes_cold <- semi_join(vst_sig_genes_cold, tidy_normecounts_names_cold, by  = "geneName") %>%
  mutate(time = as.numeric(time))

#reverse transform (treat vst as a log2 transform)
vst_sig_genes_cold <- vst_sig_genes_cold %>% mutate(reverseIntensity = 2 ^ intensity)

#select columns of interest
vst_sig_genes_cold <- vst_sig_genes_cold %>% select(-intensity)
vst_sig_genes_cold <- vst_sig_genes_cold%>% dplyr::rename(intensity =reverseIntensity) %>%
  filter(!is.na(bNum))%>%
  select(bNum, geneName, replicate, genotype,time, intensity, batch)

#save this to file so it can be used
write_csv( vst_sig_genes_cold, file="../../outputs/DESeq_results_cold.csv")

###Uncomment to run sicegar.
#Because it is time consuming to run sicegar, the results have been saved to file, which can be found in outputs/sicegar
#
#
# ##stationary
# all_SP_sicegar_list<- vst_sig_genes_SP %>%
#   split(.$geneName) %>%
#   map(safely( ~ sicegarDataMapT(
#     dat = .,
#     AIC = 5,
#     thresholdRatio = 0.65
#   ))) %>%
#   purrr::transpose() %>%
#   .$result %>%
#   unlist()
#
# #making sicegar product tidy
# all_SP_sicegar_list <- as_tibble(all_SP_sicegar_list, rownames = "geneName")
#
# #write to file to avoid rerunning sicegar
# write.csv( all_SP_sicegar_list, file="outputs/sicegar/recount_SP_sicegar_data.csv")
#
# ##osmo
# osmo_sicegar_list <- tidy_normecounts_names_osmo %>%
#   split(.$geneName) %>%
#   map(safely( ~ sicegarDataMapT(
#     dat = .,
#     AIC = 5,
#     thresholdRatio = 0.65
#   ))) %>%
#   purrr::transpose() %>%
#   .$result %>%
#   unlist()
#
# #making sicegar product tidy
# osmo_sicegar_list <- as_tibble(osmo_sicegar_list, rownames = "geneName")
#
# #write to file to avoid rerunning sicegar
# write.csv( osmo_sicegar_list, file="outputs/sicegar/osmo_sicegar_data.csv")
#
# ##cold
# all_cold_sicegar_list<- vst_sig_genes_cold %>%
#   split(.$geneName) %>%
#   map(safely( ~ sicegarDataMapT(
#     dat = .,
#     AIC = 5,
#     thresholdRatio = 0.65
#   ))) %>%
#   purrr::transpose() %>%
#   .$result %>%
#   unlist()
#
# #making sicegar product tidy
# all_cold_sicegar_list <- as_tibble(all_cold_sicegar_list, rownames = "geneName")
#
# #write to file to avoid rerunning sicegar
# write.csv( all_cold_sicegar_list, file="outputs/sicegar/recount_cold_sicegar_data.csv")


#get sicegar results from file
all_SP_sicegar_list <-read.csv("../../outputs/sicegar/recount_SP_sicegar_data.csv") %>%
  select(geneName, value)%>%
  mutate(Experiment = "Stationary Phase")
osmo_sicegar_list <-read.csv("../../outputs/sicegar/osmo_sicegar_data.csv") %>%
  mutate(Experiment = "Osmotic Shock")
all_cold_sicegar_list <-read.csv("../../outputs/sicegar/recount_cold_sicegar_data.csv") %>%
  select(geneName, value)%>%
  mutate(Experiment = "Cold Shock")
