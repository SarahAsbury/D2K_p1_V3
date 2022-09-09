#PhyloseqV5
#16aug22

# Install packages --------------------------------------------------------
devtools::install_github("SarahAsbury/retainresolve")
detach("package:retainresolve", unload = TRUE)

# Import libraries --------------------------------------------------------
library(retainresolve)
library(tidyverse)
library(phyloseq)

# Clear Environment -------------------------------------------------------
rm(list = ls())

# Set wd (physeq input) ------------------------------------------------------------------
#Set wd
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5")


#Data paths
asvfile <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV5/PhyloseqInput/seqtab_nochim_transposed_JFUTSouthwestern_v34_SA21feb21.csv"
taxfile <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV5/PhyloseqInput/taxa_JFUTSouthwestern_v34_silva132.csv"
mapfile <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV5/mapdf/mapdf.csv"


# Import data -------------------------------------------------------------
# === map ===
mapdf <- read.csv(mapfile,
                 strip.white = TRUE,
                 na.strings = 'NA',
                 stringsAsFactors = FALSE) %>%
  mutate(depression = ifelse(phq9_total >= 10, "depressed", ifelse(phq9_total <=9, "control", NA)), #add  categorical clinical variables
         anxiety = ifelse(gad7_total >=10, "anxious", ifelse(gad7_total <=9, "control", NA)),
         anhedonia = ifelse(dars17_total <= 44, "anhedonia", ifelse(dars17_total > 44, "control", NA)))


dim(mapdf)
head(mapdf)

#Check that D2K IDs removed for phyloseqV5 are absent in dataset
removed.samples.v5 <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/mapdf/removed_samples.csv")
qc <- mapdf %>% filter(record %in% removed.samples.v5$record) %>% nrow()
ifelse(qc == 0, "QC passed - correct df for physeqV5 imported", "QC failed")



# === asvs ===
asvdf <- read.csv(asvfile, row.names=1) %>% select(mapdf$JF_ID)
str(asvdf)
dim(asvdf)


# === taxanomic classification ===
taxdf <- read.csv(taxfile, row.names = 1)
head(taxdf)



# Create phyloseq object --------------------------------------------------
dat <- physeq.create(asvdf = asvdf, taxdf = taxdf, mapdf = mapdf, ID.col = "JF_ID", export.asvsp = TRUE, export.wd = "/Users/sar/Dropbox/Dallas2K/PhyloseqV5/")




# Add full taxonomic ID to phyloseq taxa labels ---------------------------
dat <- sp_taxa_ID(dat)


# Add clinical phenotypes (MDD, MDD_GAD, NONE) to Phyloseq object ------------------------------------
#Remove anhedonia from analysis, but use anhedonia as exclusion criteria for NONE
clin.phen <- sample_data(dat) %>% data.frame() %>%
  mutate(mood = ifelse(depression == "control" & anxiety == "control" & anhedonia == "control", "NONE",
                       ifelse(depression == "depressed" & anxiety == "control", "MDD",
                              ifelse(depression == "control" & anxiety == "anxious", "GAD",
                                     ifelse(depression == "depressed" & anxiety == "anxious", "MDD_GAD", NA)
                              )
                       )
  )
  )
nrow(clin.phen) #n = 378

#Subset samples in Phyloseq object
  #V5: This is legacy code, no samples are filtered based on clinical scores.
dat.prep2 <- dat %>% subset_samples(StudyID %in% clin.phen$StudyID) #Only include StudyIDs that made it through previous clin phenotype filtering on mapdf_alpha (QC was done on just df; not phyloseq )
dat.prep2
#QC that StudyIDs are aligned between mapdf_alpha (which contains mood column) and Phyloseq object
qc.df <- sample_data(dat.prep2) %>% data.frame() %>%
  mutate(StudyID_clin.phen = clin.phen$StudyID) %>%
  mutate(QC = ifelse(StudyID_clin.phen == StudyID, "Yes", "No")) %>% count('QC')
qc.df
qc.pass <- ifelse("No" %in% qc.df$QC, FALSE, TRUE)

#Add mood column to phyloseq object
if (qc.pass == TRUE){
  print("QC passed. Mood column will be added to Phyloseq object (dat.prep3)")
  dat.prep3 <- dat.prep2
  sample_data(dat.prep3)$mood <- clin.phen$mood
}
if (qc.pass == FALSE){
  print("Warning: QC failed. Mood column not added.")
}
dat.prep3
dat.prep3 %>% sample_data() %>% head()
dat <- dat.prep3 #assign newly generated phyloseq object containig mood in sample data to dat


# Remove QC'ed Samples from Pipeline --------------------------------------
#Samples:
#D2K000240 Reads >10,000 & Low alpha diversity
#D2K000106 Reads ~ 17.5k & low alpha diversity
dat
dat_pr <- subset_samples(dat, !StudyID %in% c("JF1789", "JF1793")) #Remove QC'ed samples (Low counts and/or α diversity)
dat_pr


# Remove Host Sequences ---------------------------------------------------
dat_pr = subset_taxa(dat_pr,
                     !(Kingdom == 'Eukaryota' |
                         is.na(Phylum) |
                         (!is.na(Family) &
                            Family == 'Mitochondria')))
dat_pr





# Subset visit = 1 --------------------------------------------------------
dat_all <- dat_pr
dat_pr <- subset_samples(dat_pr, visit == 1)
write.csv(file = "/Users/sar/Dropbox/Dallas2K/PhyloseqV5/final_sample_list.csv",
  dat_pr %>% sample_data %>% data.frame %>% select(record) %>% arrange(record))


# Retain resolve ----------------------------------------------------------
retain.resolve_genus(physeq = dat_pr, dir = "/Users/sar/Dropbox/Dallas2K/PhyloseqV5", export = TRUE)

