#D2K_Glom_Aldex
#Aldex differential abundance testing of Resolve to Genus glomming
#Phyloseq V4 6may21
#Phyloseq V5 16aug22

# Clear Environment -------------------------------------------------------
#Clear environment
rm(list = ls(all.names = TRUE))


# Load Packages -----------------------------------------------------------
library(ALDEx2)
library(tidyverse)
library(microbiome)
library(janitor)
library(ggpubr)
# Load colour palette -----------------------------------------------------
colors_contrast <- c("#114f89ff", "#c11319ff", "#41721eff", "#e8b91eff", "#5c1f7cff", "#140909ff",
                     "#D991BA", "#65743A", "#b1f8f2", "#32E875", "#90A9B7", "#820263")



# Glom to genus  ----------------------------------------------------------
load("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/Physeq_glommed_objects/physeq_ResolveGenus.RData")
dat_glom #Filtered phyloseq object - Counts
dat_rel_glom #Filtered phyloseq object - Rel Abundance
dat_glom_aldex #Contains OTHER taxa so that seq depth is preserved for CLR calculation 

# Set wd -------------------------------------------------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/aldex")


# Setup physeq for aldex ------------------------------------------------------------


factors <- c("Illumina.Plate.Number", "Date.on.Sample", "dsex", "hisp_ethnicity", "race", "depression", "anxiety", "anhedonia", "mood", "dx")
numeric <- c("visit", "age", "dars17_total", "gad7_total", "phq9_total")


#Inspect phyloseq object
dat_glom <- dat_glom_aldex
dat_glom
otu_table(dat_glom)[1:5,1:5]
dat_glom %>% sample_data() %>% head()

# Aldex - Continuous ------------------------------------------------------
factors <- c("Illumina.Plate.Number", "Date.on.Sample", "dsex", "hisp_ethnicity", "race", "depression", "anxiety", "anhedonia", "mood", "dx")
numeric <- c("visit", "age", "dars17_total", "gad7_total", "phq9_total")
#Inspect phyloseq object
dat_glom <- dat_glom_aldex
dat_glom
otu_table(dat_glom)[1:5,1:5]
dat_glom %>% sample_data() %>% head()

#Age
aldex.dat <- dat_glom %>% phyloseq::subset_samples(!is.na(age))
otu <- aldex.dat %>% otu_table()
cond <- aldex.dat %>% sample_data() %>% data.frame() %>% 
  mutate_at(factors, as.factor) %>% mutate_at(numeric, as.numeric) %>%
  dplyr::select(StudyID, age)

age.result <- aldex(otu,
                    model.matrix(~1 + age, data = cond), 
                    test = "glm") 

write.csv(age.result, "age_aldex.csv")

#GAD
aldex.dat <- dat_glom %>% phyloseq::subset_samples(!is.na(gad7_total))
otu <- aldex.dat %>% otu_table()
cond <- aldex.dat %>% sample_data() %>% data.frame() %>% 
  mutate_at(factors, as.factor) %>% mutate_at(numeric, as.numeric) %>%
  dplyr::select(StudyID, gad7_total)

gad.result <- aldex(otu,
                     model.matrix(~1 + gad7_total, data = cond), 
                     test = "glm") 

write.csv(gad.result, "gad_aldex.csv")

#GAD + age
aldex.dat <- dat_glom %>% phyloseq::subset_samples(!is.na(gad7_total) & !is.na(age))
otu <- aldex.dat %>% otu_table()
cond <- aldex.dat %>% sample_data() %>% data.frame() %>% 
  mutate_at(factors, as.factor) %>% mutate_at(numeric, as.numeric) %>%
  dplyr::select(StudyID, gad7_total, age)

gad_age_interact.result <- aldex(otu,
                               model.matrix(~1 + age + gad7_total + age:gad7_total, data = cond), 
                               test = "glm") 
write.csv(gad_age_interact.result, "gadage_aldex.csv")

#PHQ 
aldex.dat <- dat_glom %>% phyloseq::subset_samples(!is.na(phq9_total))
otu <- aldex.dat %>% otu_table()
cond <- aldex.dat %>% sample_data() %>% data.frame() %>% 
  mutate_at(factors, as.factor) %>% mutate_at(numeric, as.numeric) %>%
  dplyr::select(StudyID, phq9_total)

phq.result <- aldex(otu,
                    model.matrix(~1 + phq9_total, data = cond), 
                    test = "glm") 
write.csv(phq.result, "phq_aldex.csv")

#ANH
aldex.dat <- dat_glom %>% phyloseq::subset_samples(!is.na(dars17_total))
otu <- aldex.dat %>% otu_table()
cond <- aldex.dat %>% sample_data() %>% data.frame() %>% 
  mutate_at(factors, as.factor) %>% mutate_at(numeric, as.numeric) %>%
  dplyr::select(StudyID, dars17_total)


anh.result <- aldex(otu,
                    model.matrix(~1 + dars17_total, data = cond), 
                    test = "glm") 

write.csv(anh.result, "anh_aldex.csv")


#Continuous mood - depression/anxiety
aldex.dat <- dat_glom %>% phyloseq::subset_samples(!is.na(gad7_total) & !is.na(phq9_total))
otu <- aldex.dat %>% otu_table()
cond <- aldex.dat %>% sample_data() %>% data.frame() %>% 
  mutate_at(factors, as.factor) %>% mutate_at(numeric, as.numeric) %>%
  dplyr::select(StudyID, gad7_total, phq9_total)

cont_mood_int.result <- aldex(otu,
                          model.matrix(~1 + gad7_total + phq9_total + gad7_total:phq9_total, data = cond), 
                          test = "glm") 

write.csv(cont_mood_int.result, "gadphq_aldex.csv")

#Continuous mood - depression/anxiety AND anhedonia
aldex.dat <- dat_glom %>% phyloseq::subset_samples(!is.na(gad7_total) & !is.na(phq9_total) & !is.na(dars17_total))
otu <- aldex.dat %>% otu_table()
cond <- aldex.dat %>% sample_data() %>% data.frame() %>% 
  mutate_at(factors, as.factor) %>% mutate_at(numeric, as.numeric) %>%
  dplyr::select(StudyID, gad7_total, phq9_total, dars17_total)

cont_mood_anh_effect.result <- aldex(otu,
                                 model.matrix(~1 + gad7_total + phq9_total + dars17_total, data = cond), 
                                 test = "glm") 
write.csv(cont_mood_anh_effect.result, "gadphqanh_aldex.csv")


#Continuous mood - depression/anhedonia interaction 
aldex.dat <- dat_glom %>% phyloseq::subset_samples(!is.na(gad7_total) & !is.na(phq9_total) & !is.na(dars17_total))
otu <- aldex.dat %>% otu_table()
cond <- aldex.dat %>% sample_data() %>% data.frame() %>% 
  mutate_at(factors, as.factor) %>% mutate_at(numeric, as.numeric) %>%
  dplyr::select(StudyID, gad7_total, phq9_total, dars17_total)

cont_mood_anh_int.result <- aldex(otu,
                                     model.matrix(~1 + phq9_total + dars17_total + phq9_total:dars17_total, data = cond), 
                                     test = "glm") 
write.csv(cont_mood_anh_int.result, "phqanh_aldex.csv")


# Microbiome CLR abundances --------------------------------------
#List of all taxa except other
taxa <- dat_glom %>% tax_table() %>% data.frame() %>% 
  rownames_to_column(var = "OTU") %>% dplyr::select(OTU)

#Convert to CLR
dat_clr <- microbiome::transform(dat_glom_aldex, "clr") %>% 
  prune_taxa(taxa$OTU, .)

otu_table(dat_clr)[1:5, 1:5]

#ASV wide df
df_clr <- dat_clr %>% psmelt() %>%
  group_by(StudyID) %>%
  arrange(desc(Abundance)) %>% 
  dplyr::select(c(OTU,Sample, Abundance, StudyID, mood, age, dsex, race, phq9_total, gad7_total, dars17_total)) %>%
  arrange(Sample, desc(Abundance))%>%
  group_by(Sample) %>%
  tidyr::pivot_wider(values_from = Abundance, names_from = OTU) %>%
  dplyr::select(-StudyID)

# Functions - Continuous Taxa graphs -------------------------------------------------
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
{
  cont.glm <- aldex_result %>% 
    rename(p.BH = paste0("model.", voi, ".Pr...t...BH"),
           X.Intercept = X.Intercept..Estimate,
           coef = paste0("model.", voi, ".Estimate"))%>% 
    select(p.BH, X.Intercept, coef) %>% 
    rownames_to_column(var = "OTU") 
  return(cont.glm)
}

aldex.glm.plot <- function(aldex, #aldex results 
                           trait, #variable of interest
                           export_params = c(4000, 4000, 300, 1), #width, height, res, ncol
                           tsize = 10,#textsize 
                           p.cutoff = 0.7, #Minimum adjusted p-value required for plotting
                           eqnsize = 2 #text size of equation
                           )
  #Scatter plots with linear model line, equation and p-value
  #Plots all taxa that meet the minimum adjusted p-value cut-off 
  #Dependant on aldex.contglm.extract function
  #Exports output to TIFF file 
{
  cont_graphing <- aldex.contglm.extract(aldex_result = aldex, 
                                         voi = paste(trait)) %>% 
    filter(p.BH < p.cutoff) %>% arrange(p.BH)
    
  plot.df <- data.frame()
  plot.list <- list()
  plotlist <- list()
  for (i in 1:nrow(cont_graphing)){
    #Introduce taxa
    y.name <- cont_graphing[i,1]
    print(y.name)
    
    #Abbreviate taxa for axis label
    taxa_abbreviated <- sp_genus(dat_glom_aldex, y.name)
    
    #Clean variable names 
    y.name <- y.name %>% make_clean_names()
    plot_df <- rename_with(df_clr, make_clean_names)
    
    #Select columns necessary for plotting df 
    plot_df <- plot_df %>% select(trait, y.name) 
    
    
    
    #Intercept sign
    int.sign <- ifelse(cont_graphing[i, 3] >= 0, 
                       "+", 
                       "-")
    
    #Plot
    plot_name <- paste(trait, i, sep = ".")
    plot.list <- rbind(plot.list, plot_name)
    
    p <- ggscatter(plot_df, 
                   x = paste(trait), 
                   y = paste(y.name),
                   ylab = "Abundance (CLR)",
                   title = paste(taxa_abbreviated),
                   color = "white",
                   fill = "white") + 
      font("title", size = tsize, face = "bold") +
      font("xlab", size = tsize) + 
      font("ylab", size = tsize) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.margin = unit(c(-0.5,2,-0.5,2), "cm"))
    
    p <- p %>% ggadd("point", alpha = 0.4, size = 1.8, color = "gray26")
    p <- p + geom_abline(intercept = cont_graphing[i, 3], slope = cont_graphing[i, 4], color = "grey26") + 
      annotate(geom = "text", 
               label = paste0("y=", round(cont_graphing[i,4], digits = 2), "x", int.sign, abs(round(cont_graphing[i, 3], digits = 1)),
                              "  ", "p=", 
                              ifelse(cont_graphing[i, 2] < 0.001,
                                     "<0.001",
                                     round(cont_graphing[i, 2], 3))
                              ),
               y = max(eval(parse(text = paste0("plot_df$", paste(y.name)))), na.rm = TRUE)* 1.2,
               x = max(eval(parse(text = paste0("plot_df$", paste(trait)))), na.rm = TRUE) * 0.5, 
               size = eqnsize)
  
    #plot list for loop
    plotname <- paste("plot", i, sep  = ".")
    assign(paste(plotname), p)
    plot.df <- rbind(plot.df, plotname)
    plotlist[[i]] <- p
  }
  
  p <- cowplot::plot_grid(plotlist = plotlist, scale = 0.8, ncol = export_params[4])
  print(p)
  
  tiff(paste0(trait, "_glm.tiff"), width = export_params[1], height = export_params[2], res = export_params[3])
  print(p)
  dev.off()
  
}





# ALDEX continous taxa graphs - paper ---------------------------------------------------------------
p.cutoff.fun <- function(df, #object; aldex results df
                         BH.col.name, #characteric; name of BH corrected p-value  column
                         n #numeric; number of graphs to include
){
  cutoff <- df %>% arrange(get(BH.col.name)) %>% slice_min(get(BH.col.name), n = n) %>% pull(get(BH.col.name)) %>% .[n]
  return(cutoff)
}

#Export params:
w <- 1750

#Graphs: 
aldex.glm.plot(aldex = phq.result, 
               trait = "phq9_total", 
               export_params = c(w, 5000, 300, 1), #width, height, res, ncol
               tsize = 15,
               eqnsize = 5, 
               p.cutoff = p.cutoff.fun(phq.result, "model.phq9_total.Pr...t...BH", n = 7)) 


aldex.glm.plot(aldex = gad.result, 
               trait = "gad7_total", 
               export_params = c(w, 5000, 300, 1), #width, height, res, ncol
               p.cutoff = p.cutoff.fun(gad.result, "model.gad7_total.Pr...t...BH", n = 7), 
               tsize = 15, 
               eqnsize = 5) 

aldex.glm.plot(aldex = age.result, 
               trait = "age", 
               export_params = c(w, 5000, 300, 1), #width, height, res, ncol
               p.cutoff = p.cutoff.fun(age.result, "model.age.Pr...t...BH", n = 7),
               tsize = 15, 
               eqnsize = 5) 

aldex.glm.plot(aldex = anh.result, 
               trait = "dars17_total", 
               export_params = c(w, 5000, 300, 1), #width, height, res, ncol
               p.cutoff = p.cutoff.fun(anh.result, "model.dars17_total.Pr...t...BH", n = 7),
               tsize = 15, 
               eqnsize = 5)
