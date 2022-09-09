#Updated 16aug22 PhyseqV5

# Install packages --------------------------------------------------------
devtools::install_github("SarahAsbury/BioDataTools")
detach("package:BioDataTools", unload = TRUE)


devtools::install_github("SarahAsbury/retainresolve")
detach("package:retainresolve", unload = TRUE)

# Clear environment -------------------------------------------------------
rm(list = ls())


# Load packages -----------------------------------------------------------
library(BioDataTools)
library(janitor)
library(phyloseq)
library(retainresolve)

# Load packages - RF -------------------------------------------------------
myPackages <- c("GGally", "e1071", "cowplot","randomForest","tidyverse","nnet","ROCR",
                "Hmisc", "NCmisc", "phyloseq", "forecast", "janitor", "ggpubr", "cowplot", "ViewPipeSteps",
                "lubridate", "logr")
tryCount <- 0

while( !all(myPackages %in% (.packages())) ){

  try(require(GGally))
  try(require(e1071))
  try(require(cowplot))
  try(require(randomForest))
  try(require(tidyverse))
  try(require(nnet))
  try(require(ROCR))
  try(require(Hmisc))
  try(require(NCmisc))
  try(require(phyloseq))
  try(require(forecast))
  try(require(janitor))
  try(require(ggpubr))
  try(require(cowplot))
  try(require(ViewPipeSteps))
  try(require(logr))
  try(require(lubridate))

  tryCount <- tryCount + 1

  if( !all(myPackages %in% (.packages()))  ){
    cat(paste0("Failure: ", tryCount, "\n"))
    cat("Failed to load: ")
    cat(myPackages[ !myPackages %in% (.packages()) ])
    cat("\n")
  } else {
    print(paste0("Packages loaded successfully!"))
  }

  Sys.sleep(5)

}



# Load glom ---------------------------------------------------------------
load("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/Physeq_glommed_objects/physeq_ResolveGenus.RData")
dat_glom_aldex


# Prepare df --------------------------------------------------------------
dat_clr <- dat_glom_aldex %>% microbiome::transform("clr") %>% sp_taxa_ID(taxa.col = c("Family", "Genus"), remove = TRUE, remove.delimiter = "_")

x <- psmelt(dat_clr) %>%
  group_by(StudyID) %>%
  dplyr::select(c(OTU,Sample, Abundance, StudyID, gad7_total, phq9_total, dars17_total, age)) %>%
  arrange(Sample, desc(Abundance))%>%
  group_by(Sample) %>%
  tidyr::pivot_wider(values_from = Abundance, names_from = OTU) %>%
  column_to_rownames("Sample") %>%
  dplyr::select(-StudyID) %>% select(-other_Other_Other) %>% clean_names()

response <- c("gad7_total", "phq9_total", "dars17_total", "age")
response.list <- response




# Run rf ------------------------------------------------------------------

for(i in 1:length(response)){
  response.var <- response[i]
  other.var<- response[-i]
  D2K <- x %>% select(-other.var)

  rf_standard(rf.type = "reg", vpred = response.var, df = D2K, dir = "/Users/sar/Dropbox/Dallas2K/PhyloseqV4/rf",
              experiment.note = paste("Predict", response.var, "(continous) from CLR transformed 16S abundances. PhyloseqV5."),
              rf.param = c(dataframe.name = "D2K", predictors.name = response.var)
              )

}


# Edit variable importance graph formatting ----------------------------------------
# === STOPPED HERE; wait for RF to run === 

vi.gad7 <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/rf/D2K_2022-08-16_gad7_total_gad7_total/imp_aggregate.csv")
vi.phq9 <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/rf/D2K_2022-08-17_phq9_total_phq9_total/imp_aggregate.csv")
vi.age <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/rf/D2K_2022-08-18_age_age/imp_aggregate.csv")
vi.dars <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/rf/D2K_2022-08-17_dars17_total_dars17_total/imp_aggregate.csv")


head(vi.gad7)
head(vi.phq9)
head(vi.age)
head(vi.dars)




# === Export varimp plots === 

#Set wd
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/rf/var_importance")

#Set export parameters
w <- 2500
h <- 2500
res <- 300

#plot parameters
top.taxa <- 20

#plots
tiff("gad7_varimp.tiff", width = w, height = h, res = res)
varimport_plot(varimp = vi.gad7, rf.type = "reg", metric = "gini", selection_type = "random_top", top = top.taxa)
dev.off()

tiff("phq9_varimp.tiff", width = w, height = h, res = res)
varimport_plot(varimp = vi.phq9, rf.type = "reg", metric = "gini", selection_type = "random_top", top = top.taxa)
dev.off()

tiff("age_varimp.tiff", width = w, height = h, res = res)
varimport_plot(varimp = vi.age, rf.type = "reg", metric = "gini", selection_type = "random_top", top = top.taxa)
dev.off()

tiff("dars_varimp.tiff", width = w, height = h, res = res)
varimport_plot(varimp = vi.dars, rf.type = "reg", metric = "gini", selection_type = "random_top", top = top.taxa)
dev.off()


#MSE plots
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/rf/var_importance/mse")
tiff("gad7_varimp.tiff", width = w, height = h, res = res)
varimport_plot(varimp = vi.gad7, rf.type = "reg", metric = "mse", selection_type = "random_top", top = top.taxa)
dev.off()

tiff("phq9_varimp.tiff", width = w, height = h, res = res)
varimport_plot(varimp = vi.phq9, rf.type = "reg", metric = "mse", selection_type = "random_top", top = top.taxa)
dev.off()

tiff("age_varimp.tiff", width = w, height = h, res = res)
varimport_plot(varimp = vi.age, rf.type = "reg", metric = "mse", selection_type = "random_top", top = top.taxa)
dev.off()

tiff("dars_varimp.tiff", width = w, height = h, res = res)
varimport_plot(varimp = vi.dars, rf.type = "reg", metric = "mse", selection_type = "random_top", top = top.taxa)
dev.off()

