#making all different groups of genesets for GO analysis and comparison to other papers

#get genomic info
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

all_SP_sicegar_list <-read.csv("../../outputs/sicegar/recount_SP_sicegar_data.csv") %>%
  select(geneName, value)%>%
  mutate(Experiment = "Stationary Phase")
osmo_sicegar_list <-read.csv("../../outputs/sicegar/osmo_sicegar_data.csv") %>%
  mutate(Experiment = "Osmotic Shock")
all_cold_sicegar_list <-read.csv("../../outputs/sicegar/recount_cold_sicegar_data.csv") %>%
  select(geneName, value)%>%
  mutate(Experiment = "Cold Shock")

#genes shared by a just 2 of the stressors 
shared_osmo_sp <- semi_join(osmo_sicegar_list, all_SP_sicegar_list, by = "geneName")
shared_just_osmo_sp <- anti_join(shared_osmo_sp, all_cold_sicegar_list, by = "geneName")

shared_osmo_cold <- semi_join(osmo_sicegar_list, all_cold_sicegar_list, by = "geneName")
shared_just_osmo_cold <- anti_join(shared_osmo_cold, all_SP_sicegar_list, by = "geneName")

shared_SP_cold <- semi_join(all_SP_sicegar_list,all_cold_sicegar_list , by = "geneName")
shared_just_SP_cold <- anti_join(shared_SP_cold, osmo_sicegar_list, by = "geneName")

#genes shared in all three
shared_all_three <- semi_join(all_cold_sicegar_list, shared_osmo_sp, by = "geneName") %>%
  dplyr::rename(name = geneName)
shared_all_three <- semi_join(gene_info,shared_all_three, by = "name" )

#genes in just cold
cold_not_sp <- anti_join(all_cold_sicegar_list, all_SP_sicegar_list, by = "geneName")
just_cold <- anti_join(cold_not_sp, osmo_sicegar_list, by = "geneName") %>% 
  rename(name = geneName)
just_cold <- semi_join(gene_info,just_cold, by = "name" ) %>% 
  select(GeneID)



