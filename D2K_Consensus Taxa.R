#D2K
#Taxa summary table
#Identify taxa that have been pulled out of multiple analysis tools. I
#Then identify overlap in sig taxa across multiple tools.
#Aldex, RF, β-diversity
#MDD, MDD_GAD


#Updated for PhyseqV3-paper 6jan22
#Updated for Paper_v2 9feb22
#Updated for phyloseqV4 12may22

#SA updated 23aug22




# Install packages --------------------------------------------------------
devtools::install_github("SarahAsbury/BioDataTools")
detach("package:BioDataTools", unload = TRUE)
# Load Packages -----------------------------------------------------------
#Load packages
library(plyr)
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(scales)
library(gtools)
library(BioDataTools)


# Clear Environment -------------------------------------------------------
#Clear environment
rm(list = ls(all.names = TRUE))


# Load glom to genus ------------------------------------------------------
load("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/Physeq_glommed_objects/physeq_ResolveGenus.RData")
dat_glom #Filtered phyloseq object - Counts
dat_rel_glom #Filtered phyloseq object - Rel Abundance
dat_glom_aldex #Contains OTHER taxa so that seq depth is preserved for CLR calculation


# Functions - Microbiome data toolbox ---------------------------------------------------------------
sp_genus <- function(physeq, #object with reference taxa table
                     taxa #To convert. Must contain species code (sp) assigned by phyloseq (character list)
)
  #Use species code (e.g Sp1) to extract clean genus-level names from phyloseq object
  #Extracts family_NA when genus is not available
{
  taxdf <- physeq %>% tax_table %>% data.frame() %>% rownames_to_column(var = "OTU") %>%
    mutate(sp = gsub("_.*", "", OTU)) %>%
    mutate(genus_level = ifelse(is.na(Genus),
                                paste(sp, Family, Genus, sep = "_"),
                                paste(sp, Genus, sep = "_")
    )
    ) %>% select(sp, genus_level)

  taxa.in <- taxa %>% data.frame %>% mutate(sp = gsub("_.*", "", taxa))
  taxa.out <- left_join(taxa.in, taxdf, by = "sp")
  return(taxa.out$genus_level)
}

aldex.contglm.extract <- function(aldex_result, #aldex glm output
                                  voi #variable of interest
)
  #Extracts BH p-value, glm coefficient (slope), and glm x-intercept from aldex results for each taxon
  #Renames columns to human read-able columns
  #Taxa names must be rownames
{
  cont.glm <- aldex_result %>%
    rename(p = paste0("model.", voi, ".Pr...t.."),
           p.BH = paste0("model.", voi, ".Pr...t...BH"),
           X.Intercept = X.Intercept..Estimate,
           coef = paste0("model.", voi, ".Estimate"))%>%
    select(p, p.BH, X.Intercept, coef) %>%
    rownames_to_column(var = "OTU")
  return(cont.glm)
}





# Import WCNA ------------------------------------------------------------
#Import WCNA module membership
#Brown module
module.edit <- function(df, color)
  #filter to module colour
  #format edits
  {
  df <- df %>% rename(taxa = X) %>% filter(moduleColors == color) %>% rename(moduleColor = moduleColors) %>%
    mutate(across(starts_with("MM"), ~round(.x,digits = 2))) %>%
    mutate(across(starts_with("p.MM"), ~scientific(.x,digits = 3))) %>%
    mutate(taxa = sp_genus(taxa = taxa, physeq = dat_glom))
  return(df)
}

brown.member <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/brown/browntaxa_membership.csv") %>% module.edit("brown")
head(brown.member)


#Import WCNA correlations - all brown module taxa
wcna.alledit <- function(df){
  p.name <- df %>% select(starts_with("p.GS")) %>% colnames
  p.sub <- p.name %>% sub("GS.", "", .)
  sig.name <- p.name %>% sub("p.GS.", "sig.", .)

  df.return <- df %>% rename(taxa = X) %>% select(-moduleColors) %>%
    mutate(across(starts_with("GS"), ~round(.x,digits = 2))) %>%
    rename_at(vars(starts_with("GS")), funs(sub("GS.", "", .))) %>%
    mutate(p = scientific(eval(parse(text = p.name)), digits =3)) %>%
    relocate(p, .before = paste0(p.name)) %>%
    rename_with(.fn = ~paste0(p.sub), .cols = p) %>%
    mutate(across(starts_with("p.GS"), ~p.adjust(.x, method = "BH"), digits = 3)) %>%
    mutate(sig = stars.pval(eval(parse(text = p.name)))) %>%
    rename_with(.fn = ~paste0(sig.name), .cols = sig) %>%
    mutate(across(starts_with("p.GS"), ~ifelse(.x < 0.001, "<0.001", round(.x, digits = 3)))) %>%
    rename_at(vars(starts_with("p.GS")), funs(sub("p.GS.", "p.BH.", .))) %>%
    mutate(taxa = sp_genus(taxa = taxa, physeq = dat_glom)) %>%

    return(df.return)
}

wcna.all.gad <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/brown/gad7_total/brown_gad7_total_ALLclinicaltaxa.csv") %>% wcna.alledit
wcna.all.phq <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/brown/phq9_total/brown_phq9_total_ALLclinicaltaxa.csv") %>% wcna.alledit
wcna.all.dars <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/brown/dars17_total/brown_dars17_total_ALLclinicaltaxa.csv") %>% wcna.alledit


head(wcna.all.gad)
head(wcna.all.phq)
head(wcna.all.dars)

# WCNA - Brown Module Membership and Taxa Significance --------------------
#Supplementary table
#Module Membership and Taxa Significance table
rm.wcna.mmsigtab <- function(df){
  out <- df %>% select(-c(MMturquoise, MMblue, moduleColor,
                           starts_with("p.MM")))
}

wcna.mmsigtab <- join_all(list(brown.member, wcna.all.gad, wcna.all.phq, wcna.all.dars), by = "taxa", type = "left") %>%
  rm.wcna.mmsigtab
head(wcna.mmsigtab)




# Export Module Membership Taxa Significance Table ------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25")
write.csv(wcna.mmsigtab, "brownMM_taxsig_table.csv", row.names = FALSE)



# Import WCNA hub taxa ----------------------------------------------
#Brown module
wcna.edit <- function(df){
  return(df %>% data.frame() %>% rename(taxa = X) %>% mutate(taxa = sp_genus(taxa = taxa, physeq = dat_glom)))
}
wcna.hub <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/brown/brown_hubtaxa.csv") %>% wcna.edit %>% select(taxa)


# Import RF results --------------------------------------------------------------
rf.taxa20 <- function(df){
  return(df %>% slice_max(n = 20, order_by = IncNodePurity_name) %>% select(-c(X.IncMSE_name)) %>% mutate(taxa = sp_genus(taxa = predictors, physeq = dat_glom))
         )
}

rf.rank <- function(df){
  return(df %>% data.frame %>% arrange(desc(IncNodePurity_name)) %>% select(-c(X.IncMSE_name)) %>% mutate(taxa = sp_genus(taxa = predictors, physeq = dat_glom)) %>%
    rownames_to_column(var = "rank")
  )
}

rf.gad.dir <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV5/rf/D2K_2022-08-16_gad7_total_gad7_total/imp_aggregate.csv"
rf.phq.dir <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV5/rf/D2K_2022-08-17_phq9_total_phq9_total/imp_aggregate.csv"
rf.dars.dir <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV5/rf/D2K_2022-08-17_dars17_total_dars17_total/imp_aggregate.csv"

rf.gad <- read.csv(rf.gad.dir) %>% rf.taxa20
rf.phq <- read.csv(rf.phq.dir) %>% rf.taxa20
rf.dars <- read.csv(rf.dars.dir) %>% rf.taxa20


rf.gad.rank <- read.csv(rf.gad.dir) %>% rf.rank
rf.phq.rank <- read.csv(rf.phq.dir) %>% rf.rank
rf.dars.rank <- read.csv(rf.dars.dir) %>% rf.rank


head(rf.gad)
head(rf.dars.rank)

# Import Aldex  ------------------------------------------------------------
aldex.p5 <- function(df){
  return(df %>% rename(taxa = OTU) %>% filter(p <= 0.05) %>% mutate(taxa = sp_genus(taxa = taxa, physeq = dat_glom)))
}

aldex.gad <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/aldex/gad_aldex.csv") %>% data.frame %>% column_to_rownames(var = "X") %>%
  aldex.contglm.extract(voi = "gad7_total") %>% aldex.p5
aldex.phq <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/aldex/phq_aldex.csv") %>% data.frame %>% column_to_rownames(var = "X") %>%
  aldex.contglm.extract(voi = "phq9_total")%>% aldex.p5
aldex.dars <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/aldex/anh_aldex.csv") %>% data.frame %>% column_to_rownames(var = "X") %>%
  aldex.contglm.extract(voi = "dars17_total") %>% aldex.p5


# Retain-resolve taxa taxonomy, mra, and prev - Combine and Export ---------------------------------------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5")
glom.taxa.levels <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/retain-resolve/glom_spID_levels.csv")
mra.prev.genus <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/retain-resolve/mra_prevalence/mraprev_glomgenus.csv") %>% rename(taxa = OTU)
mra.prev.asv <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/retain-resolve/mra_prevalence/mraprev_original.csv") %>% rename(taxa = OTU)

glom.taxonomy <- glom.taxa.levels %>%
  left_join(mra.prev.genus, by = "taxa") %>%
  left_join(mra.prev.asv, by = "taxa", suffix = c("", ".asv")) %>%
  left_join(dat_glom %>% tax_table %>% data.frame %>% rownames_to_column(var = "taxa"), by = "taxa") %>%
  mutate(mean_relabund = ifelse(is.na(mean_relabund), mean_relabund.asv, mean_relabund),
         prev = ifelse(is.na(prev), prev.asv, prev)
  ) %>% relocate(prev, .before = mean_relabund) %>%
  select(-c(prev.asv, mean_relabund.asv)) %>%
  mutate(spID = gsub(x = taxa, pattern = "_.*", replacement = "")) %>% relocate(spID, .before = Kingdom)

write.csv(glom.taxonomy, "glom_taxonomy.csv", row.names = FALSE)

# Taxa Consensus Functions ----------------------------------------------------------
convert_numeric <- function(x){
  ifelse(x == 1, "YES",
         ifelse(x == 0, "", NA))
}

aldex.format <- function(df){
  df <- df %>% mutate(p = ifelse(p < 0.001, "<0.001", as.character(round(p, digits = 3)))) %>%
    mutate(p.BH = round(p.BH, digits = 3)) %>% mutate(coef = round(coef, digits = 2)) %>%
    relocate(coef, .before = p) %>% relocate(X.Intercept, .before = p)

  return(df)
}

wcna.filter <- function(wcna = wcna.mmsigtab,
                        module.members,
                        module,
                        voi.in){
  print(voi.in)
  #Naming and filtering variables
  MMmodule <- paste0("MM", module)
  quant.50 <- quantile(module.members %>% pull(eval(parse(text = MMmodule)))) %>% data.frame() %>%
    rename(val = ".") %>% rownames_to_column(var = "quartile") %>% filter(quartile == "50%") %>% pull(val)
  sig.name <- paste0("sig.", voi.in)
  sig.symbols <- c("***", "**", "*")

  #Edit df
  wcna <- wcna %>%
    rename_at((vars(starts_with(paste(voi.in)))), funs(paste0("r"))) %>%
    mutate(r = abs(r)) %>%
    rename_at(MMmodule, funs(paste0("MM"))) %>%
    rename_at(vars(contains(sig.name)), funs(paste0("sig"))) %>%
    filter(MM >= quant.50 & r >= 0.2 & sig %in% sig.symbols)

  return(wcna)

}


consensus.taxa <- function(voi #voi must the way it appears in dataframe names
                           )
  #Extracts consensus taxa from WCNA, RF, and Aldex results
  {

  aldex.df <- get(paste0("aldex.", voi))
  aldex.taxa <- aldex.df %>% select(taxa) %>% pull(taxa)
  print("Extracted ALDEx2 results")
  print(head(aldex.df))

  rf.taxa <- get(paste0("rf.", voi)) %>% select(taxa) %>% pull(taxa)
  rf.rank <- get(paste0("rf.", voi, ".rank")) %>% select(-IncNodePurity_name)
  print("Extracted RF results")
  print(head(rf.rank))

  wcna.module.colour <- "brown"
  wcna.taxa <- wcna.mmsigtab %>% wcna.filter(wcna = ., module.members = brown.member, module = wcna.module.colour, voi.in = paste(voi)) %>% select(taxa) %>% pull(taxa)
  print(paste("Extracted WCNA", wcna.module.colour, "module taxa"))
  print(head(wcna.taxa))

  summary <- union(wcna.taxa, rf.taxa) %>% union(., aldex.taxa) %>%
    pipe_message("1. Joined WCNA, RF, and ALDEx2 taxa df") %>%
    data.frame %>% rename(taxa = ".") %>%
    pipe_message("2. Test for WCNA, RF, and ALDEx2 criteria") %>%
    mutate(WCNA = ifelse(taxa %in% wcna.taxa, 1, 0)) %>%
    mutate(RandomForest = ifelse(taxa %in% rf.taxa, 1, 0)) %>%
    mutate(Aldex = ifelse(taxa %in% aldex.taxa, 1, 0)) %>%
    pipe_message("3. Filter taxa achieving 2/3 analyses' criteria") %>%
    rowwise %>% mutate(criteria = sum(WCNA, RandomForest, Aldex)) %>% filter(criteria >= 2) %>% select(-criteria) %>%
    mutate(across(c("WCNA", "RandomForest", "Aldex"),convert_numeric )) %>%
    pipe_message("4. Add ALDEx2 results (β and p.BH)") %>%
    left_join(aldex.df, by = "taxa") %>%
    aldex.format %>% select(-X.Intercept) %>%
    pipe_message("5. Add Random Forest feature importance rank to taxa outside of Top20") %>%
    left_join(rf.rank, by = "taxa") %>% mutate(RandomForest = ifelse(RandomForest == "YES", RandomForest, rank)) %>% select(-c(rank, predictors)) %>%
    pipe_message("6. Add taxonomic resolution - i.e. ASV or genus") %>%
    left_join(glom.taxonomy %>% select(taxa, level) %>% rename(res = level) %>% mutate(taxa = sp_genus(physeq = dat_glom, taxa = taxa)),
              by = "taxa") %>%
    relocate(res, .after = taxa)
  return(summary)
}


# Export consensus results ----------------------------------------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/consensus")
write.csv(consensus.taxa(voi = "gad"), "consensus_gad.csv", row.names = FALSE)
write.csv(consensus.taxa(voi = "phq"), "consensus_phq.csv", row.names = FALSE)
write.csv(consensus.taxa(voi = "dars"), "consensus_anh.csv", row.names = FALSE)





# Import Aldex - Full Table -----------------------------------------------
aldex.full <- function(df, voi){
  return(df %>% data.frame %>% rename(taxa = X) %>% filter(taxa != "other") %>%
           mutate(taxa = sp_genus(taxa = taxa, physeq = dat_glom)) %>%
           column_to_rownames(var = "taxa") %>%
           aldex.contglm.extract(voi = voi) %>%
           filter(p < 0.05)
  )
}

aldex.gad <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/aldex/gad_aldex.csv") %>% aldex.full(voi = "gad7_total") %>% aldex.format() %>% rename(β = coef)
head(aldex.gad)
aldex.phq <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/aldex/phq_aldex.csv") %>% aldex.full(voi = "phq9_total") %>% aldex.format() %>%
  rename(β = coef)
head(aldex.phq)
aldex.dars <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/aldex/anh_aldex.csv") %>% aldex.full(voi = "dars17_total") %>% aldex.format() %>%
  rename(β = coef)
head(aldex.dars)


# Export aldex tables -----------------------------------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/aldextables")
write.csv(aldex.gad, "aldexgadtab.csv", row.names = FALSE)
write.csv(aldex.phq, "aldexphqtab.csv", row.names = FALSE)
write.csv(aldex.dars, "aldexdarstab.csv", row.names = FALSE)




