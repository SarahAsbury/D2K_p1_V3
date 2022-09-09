
# Import libraries --------------------------------------------------------
library(tidyverse)


# Clear Environment -------------------------------------------------------
rm(list = ls())

# setwd  ------------------------------------------------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25")

# Import data -------------------------------------------------------------
load("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/Physeq_glommed_objects/physeq_ResolveGenus.RData")
dat_glom #Filtered phyloseq object - Counts
dat_rel_glom #Filtered phyloseq object - Rel Abundance
dat_glom_aldex #Contains OTHER taxa so that seq depth is preserved for CLR calculation



meta <- dat_glom %>% sample_data %>% data.frame
meta %>% group_by(psych_tx) %>% count()


# Import WCNA -------------------------------------------------------------
wcna <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/eigentaxa.csv") %>% rename(StudyID = X)
# Subgroup analysis -------------------------------------------------------
#inspect datasets
head(meta)
head(wcna)

#join
df <- wcna %>% left_join(meta, by = "StudyID")


#plot
tiff("brown_TxSubgroup.tiff", width = 1500, height = 1500 , res = 300)
ggboxplot(data = df, x = "psych_tx", y = "MEbrown", color = "psych_tx", add = "jitter") +
          stat_compare_means(method = "wilcox.test")
dev.off()

