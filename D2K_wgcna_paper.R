#Co-occurence networks


#Updated 3jan22
#Use new paper physeq glom
#High resolution exports

#Updated 6may22
#Phyloseq V4

#Update 16aug22
#Phyloseq V5
  #Removed hub taxa stability (unstandardized)



# Load Packages ------------------------------------------------------------
library(WGCNA)
library(cowplot)
library(tidyverse)
library(ggpubr)
library(phyloseq)
library(stats)
library(janitor)
library(microbiome)
library(RCy3)
library(Hmisc)

# Clear Environment -------------------------------------------------------
rm(list = ls())


# Retain/reso to genus  ----------------------------------------------------------
load("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/Physeq_glommed_objects/physeq_ResolveGenus.RData")
dat_glom #Filtered phyloseq object - Counts
dat_rel_glom #Filtered phyloseq object - Rel Abundance
dat_glom_aldex #Contains OTHER taxa so that seq depth is preserved for CLR calculation

# Set wd -------------------------------------------------------------------
parent.dir <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna"
setwd(parent.dir)




# Abbreviate taxa names  -------------------------------------------
otu_abbrv <- dat_glom_aldex %>% tax_table() %>% data.frame() %>% select(Family, Genus) %>%
  rownames_to_column(var = "taxa") %>% mutate(sp = gsub("_.*", "", taxa)) %>%
  mutate(name = paste(sp, Family, Genus, sep = "_")) %>% select(taxa, name) %>%
  mutate(taxa = make_clean_names(taxa))



# Prep CLR df ----------------------------------------------------------------
dat_clr <- dat_glom_aldex %>% microbiome::transform("clr")
rm_col <- c("Order", "Class", "Phylum", "Kingdom")
metadata <- c("Illumina.Plate.Number", "Date.on.Sample", "record", "visit", "age", "dsex" , "hisp_ethnicity",
              "race", "depression", "anxiety", "anhedonia", "dars17_total", "gad7_total", "phq9_total", "mood", "dx", "psych_tx", "BMI")
numeric_trait <-c("dars17_total", "gad7_total", "phq9_total", "age", "BMI")


clr <- psmelt(dat_clr) %>%
  group_by(StudyID) %>% select(-all_of(rm_col)) %>%
  arrange(Sample, desc(Abundance))%>%
  group_by(Sample) %>%
  filter(OTU != "other") %>%
  mutate(sp = gsub("_.*", "", OTU)) %>% mutate(taxa = paste(sp, Family, Genus, sep = "_") %>% make_clean_names) %>%
  select(-c(sp, OTU, Family, Genus)) %>%
  tidyr::pivot_wider(values_from = Abundance, names_from = taxa) %>%
  select(-StudyID)

# Prep relabund df --------------------------------------------------------
rel <- psmelt(dat_rel_glom) %>%
  group_by(StudyID) %>% select(-all_of(rm_col)) %>%
  arrange(Sample, desc(Abundance))%>%
  group_by(Sample) %>%
  filter(OTU != "other") %>%
  mutate(sp = gsub("_.*", "", OTU)) %>% mutate(taxa = paste(sp, Family, Genus, sep = "_") %>% make_clean_names) %>%
  select(-c(sp, OTU, Family, Genus)) %>%
  tidyr::pivot_wider(values_from = Abundance, names_from = taxa) %>%
  select(-StudyID)


qc.rel_clr_taxalist <- (rel %>% colnames) %in% (clr %>% colnames) %>% all
ifelse(qc.rel_clr_taxalist == TRUE, "QC passed. Taxa names same in rel abundance and CLR datasets.", "QC failed, rel and CLR have different taxa names")


# Setup df for WGCNA ------------------------------------------------------
datExpr <- clr %>% select(-metadata) %>%
  as.matrix() %>% as.data.frame() %>% column_to_rownames(var = "Sample") %>% mutate_all(as.numeric)

datTraits <- clr %>% select(all_of(numeric_trait)) %>%
  as.matrix() %>% as.data.frame() %>% column_to_rownames(var = "Sample") %>% mutate_all(as.numeric)
traitColors <- WGCNA::numbers2colors(datTraits, signed = FALSE)

datMeta <- clr %>% select(metadata) %>%
  as.matrix() %>% as.data.frame() %>% column_to_rownames(var = "Sample")

# Identify sample outliers with hclust -----------------------------------------------------

sampletree <- hclust(dist(datExpr, method = "euclidean"), method = "ward.D2")
plot(sampletree, cex = 0.6, hang = -1)

plotDendroAndColors(sampletree, traitColors, groupLabels = colnames(datTraits))



# Function: Network Setup -----------------------------------------------------------
allowWGCNAThreads()

#Network parameters:

network_selection <- function(net_type){
  #Create directory for exports
  wd <- paste0(parent.dir, "/", net_type)
  dir.create(wd)
  setwd(wd)

  #Soft-thresholding selection
  powers = c(c(1:10), seq(from = 12, to=20, by=2)) # Choose a set of soft-thresholding powers
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = net_type) # Call the network topology analysis function


  #Graph - Setup R Window
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2))
  cex1 = 0.9

  # Graph: Scale-free topology fit
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, R²",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  abline(h=0.90,col="red") # this line corresponds to using an R^2 cut-off of h

  # Graph: Mean connectivity
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



  # === Export plots ===
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2))
  cex1 = 0.9

  #Graph - Setup R Window

  # Graph: Scale-free topology fit
  tiff("softpower_r2.tiff", width = 1200, height = 1200, res = 300)
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, R²",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  abline(h=0.90,col="red") # this line corresponds to using an R^2 cut-off of h
  dev.off()

  # Graph: Mean connectivity
  tiff("softpower_k.tiff", width = 1200, height = 1200, res = 300)
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()


}

#Adjacency, TOM and merged module calculation
network_modules <- function(net_type, minModuleSize, softPower,
                            merge_modules = TRUE,
                            export = TRUE, data = datExpr,
                            network_visibility = -2){
  #Set wd
  if(export == TRUE){
    wd <- paste(parent.dir, net_type, paste0("minmod",minModuleSize), sep = "/")
    print(wd)
    dir.create(wd)
    setwd(wd)
    wd <<- wd
  }

  #Statistics for weighted correlation network:
  datExpr <- data
  adjacency <- adjacency(datExpr, power = softPower, type = net_type)
  TOM <- TOMsimilarity(adjacency, TOMType = net_type)
  dissTOM <- 1-TOM


  #Hclust of taxa using TOM dissimilarity as distance metric
  geneTree = hclust(as.dist(dissTOM), method = "average")

  if(export == TRUE){
    sizeGrWindow(12,9)
    plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
         labels = FALSE, hang = 0.04)
  }

  #Identify modules from hclust TOM dissimilarity
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  table(dynamicMods) #module asigned to each sample is given as a whole number
  dynamicColors = labels2colors(dynamicMods) # Convert numeric lables into colors
  table(dynamicColors)

  if(export == TRUE){
    # Plot the taxa dendogram with  module colours
    sizeGrWindow(8,6)
    tiff(file = "taxaDendro.tiff", width = 3000, height = 1500, res = 300)
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors")
    dev.off()
  }

  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  head(MEs)

  #Assign variables to global environment
  assign("adjacency", adjacency, envir = .GlobalEnv)
  TOM <<- TOM
  dissTOM <<- dissTOM
  MEs <<- MEs
  dynamicColors <<- dynamicColors
  moduleColors <<- dynamicColors
  geneTree <<- geneTree

  if(merge_modules == TRUE && length(dynamicColors %>% unique) > 1){
    print("Merge modules")
    # Hclust of module eigengenes (determine if any modules are highly correlated)
    MEDiss = 1-cor(MEs) #Dissimilarity of module eigengenes
    METree = hclust(as.dist(MEDiss), method = "average")

    MEDissThres = 0.25
    if(export == TRUE){
      sizeGrWindow(7, 6) # Plot the result
      plot(METree, main = "Clustering of module eigengenes",
           xlab = "", sub = "")
      abline(h=MEDissThres, col = "red") #Plot the cut line into the dendrogram
    }

    # Automatically merge similar modules
    merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    mergedColors = merge$colors # The merged module colors
    mergedMEs = merge$newMEs # Eigengenes of the new merged modules:


    if(export == TRUE){
      #Merged modules dendogram
      sizeGrWindow(12, 9)
      tiff(file = "taxaDendro_merged.tiff", width = 4000, height = 2000, res = 400)
      plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                          c("Dynamic Tree Cut", "Merged dynamic"),
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05)
      dev.off()
    }
    moduleColors = mergedColors # Replace module colors with the newly merged colors


    # Construct numerical labels corresponding to the colors
    colorOrder = c("grey", standardColors(50))
    moduleLabels = match(moduleColors, colorOrder)-1
    MEs = mergedMEs

    #Assign variables to global environment
    moduleColors <<- moduleColors
    MEs <<- MEs

  } #merge redundant modules

  if(export == TRUE){
    #Heatmap with Modules Dendogram
    #Use to optimize minModuleSize
    #Look for "corners" of heat
    #If there are cornerns of heat within modules, the minimum module size might be too low
    plotTOM = dissTOM^(network_visibility) # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
    diag(plotTOM) = NA # Set diagonal to NA for a nicer plot
    sizeGrWindow(9,9)
    tiff("network_heatmap.tiff", heigh = 1500, width = 3000, res = 300)
    TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all taxa")
    dev.off()}

  # Define numbers of genes and samples (required for downstream analysis)
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)

  #Export number of genes and samples
  nGenes <<- nGenes
  nSamples <<- nSamples


}


#Results:
#adjacency #adjacency matrix
#TOM #TOM matrix
#dissTOM #TOM as a dissimilarity metric
#MEs %>% head #module Eigengenes
#dynamicColors %>% unique() #unmerged colors
#moduleColors %>% unique() #merged colors
#plot(geneTree) #geneTree
#wd #working directory created by the parameters
#nGenes
#nSamples
# Function: Clinical trait correlation -------------------------------------
module_traits <- function(x){
  #Set wd
  setwd(wd)


  # Recalculate Module Eigengenes with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)

  #Correlate module eigengenes with clinical traits
  head(MEs)
  head(datTraits)

  moduleTraitCor = cor(MEs, datTraits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

  moduleTraitCor
  moduleTraitPvalue

  sizeGrWindow(10,6)
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3))

  # Module - Clinical Trait Heatmap
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  #Export
  tiff("trait_heatmap_wcgna.tiff", width = 1500, height = 1500, res = 300)
  par(mar = c(6, 6.8, 3, 1))

  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = redWhiteGreen(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.6,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()

  #Assign variables to global environment
  moduleTraitCor <<- moduleTraitCor
}

# Function: Investigate module --------------------------------
investigate_module <- function(clin.trait,#Clinical variables of interest (e.g gad7_total)
                               mod #WCNA module to investigate voi
){
  #Purpose:
  # Investigate the relationship between Module Membership (kME) and the clinical signifcance of each taxa
  #Also identifies hub taxa and most clinically significant taxa within the module
  #Module membership (kME) = measures how connected taxa is to other taxa in the module. Hub taxa have the highest kME

  #Code:
  #=== Prep ===
  for (i in mod){
    #Set module
    print(i)
    module <- i
    modNames = substring(names(MEs), 3) #Remove ME from module names
    column = match(module, modNames)
    moduleTaxa = moduleColors==module

    #Create new directory for MODULE
    tempdir_parent <- paste0(wd, "/", module)
    dir.create(tempdir_parent)
    setwd(tempdir_parent)

    # === Calculate Module Membership for each taxa: ===
    #Module Membership = correlation between Module eigengene and the taxa
    geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
    MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

    names(geneModuleMembership) = paste("MM", modNames, sep="")
    names(MMPvalue) = paste("p.MM", modNames, sep="")

    head(geneModuleMembership)
    head(MMPvalue)

    member_df <- cbind(geneModuleMembership, MMPvalue, moduleColors)
    head(member_df)



    for (j in clin.trait){
      print(j)
      trait <- j
      #Create new directory for TRAIT
      tempdir <- paste0(tempdir_parent, "/", trait)
      dir.create(tempdir)
      setwd(tempdir)

      #Clinical trait data
      voi <- as.data.frame(eval(parse(text = paste0 ("datTraits$", trait))))
      names(voi) <- trait



      #Extract correlation sign (negative or positive)
      cor.sign_df <- moduleTraitCor %>% as.data.frame() %>% rownames_to_column(var = "module.col") %>%
        mutate(module.col = substring(module.col, 3)) %>%
        filter(module.col == module) %>% select(trait)
      cor.sign <- ifelse(cor.sign_df[1,1] < 0, "neg", "pos")
      print("cor.sign")


      # === Taxa - trait significance ===
      #Trait significance = correlation between taxa and variable of interest
      geneTraitSignificance = as.data.frame(cor(datExpr, voi, use = "p"))
      GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
      names(geneTraitSignificance) = paste("GS.", names(voi), sep="")
      names(GSPvalue) = paste("p.GS.", names(voi), sep="")

      traitsig_df <- cbind(geneTraitSignificance, GSPvalue, moduleColors)
      head(traitsig_df)


      #=== Hub taxa ===
      #Top 5% taxa for module membership
      top.5 = (nrow(member_df %>% filter(moduleColors == module)) * 0.05) %>% floor #Goal: identify top 5% connected genes in module. This Line: find how many taxa represent 5% of module.
      hubtaxa <- member_df %>% filter(moduleColors == module) %>% arrange(desc(eval(parse(text = paste0("MM", module))))) %>%
        slice_max(order_by = eval(parse(text = paste0("MM", module))),n = top.5)

      #=== Top significant taxa for clinical trait ===
      #Top 20% of taxa with R^2 > 0.2
      colname <- paste0("GS.", names(voi))
      top.20 <- (nrow(traitsig_df %>% filter(moduleColors == module)) * 0.2) %>% round(., digits = 0) #Goal: identify top 10% connected genes in module. This Line: find how many taxa represent 10% of module.

      if(cor.sign == "neg"){
        top_taxa <- traitsig_df %>% filter(moduleColors == module) %>% arrange((eval(parse(text = colname)))) %>%
          slice_max(order_by = desc((eval(parse(text = colname)))),n = top.20) %>% filter(abs(eval(parse(text = colname))) > 0.2)
        print("cor.sign = negative")
      }
      if(cor.sign == "pos"){
        top_taxa <- traitsig_df %>% filter(moduleColors == module) %>% arrange(desc(eval(parse(text = colname)))) %>%
          slice_max(order_by = eval(parse(text = colname)), n = top.20) %>% filter(abs(eval(parse(text = colname))) > 0.2)
        print("cor.sign = positive")
      }

      #=== Exports ===
      #Correlation indicates that genes highly significantly associated with a trait of often the most important elements of the modules associated with the trait
      sizeGrWindow(7, 7)
      par(mfrow = c(1,1))
      tiff("modulemembership_Abscorr.tiff", width = 1500, height = 1500, res = 300)
      verboseScatterplot(abs(geneModuleMembership[moduleTaxa, column]),
                         abs(geneTraitSignificance[moduleTaxa, 1]),
                         xlab = paste("Module membership in", module, "module"),
                         ylab = paste("Absolute taxa significance for", trait),
                         main = paste("Module membership vs. taxa significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()

      tiff("modulemembership_corr.tiff", width = 1500, height = 1500, res = 300)


      df <- data.frame(abs(geneModuleMembership[moduleTaxa, column]),
                    geneTraitSignificance[moduleTaxa, 1]) %>% mutate_all(as.numeric)
      colnames(df) <- c("taxa_mod_membership", "clinical_correlation")

      p <- ggscatter(x = "taxa_mod_membership",
                y = "clinical_correlation",
                data = df,
                add = "reg.line",  # Add regressin line
                add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson", label.sep = "\n",
                                      label.x = max(df$taxa_mod_membership)* 0.8,
                                      label.y = ifelse(cor.sign == "pos",
                                                       min(df$clinical_correlation) * 0.8,
                                                       max(df$clinical_correlation) * 0.8)
                                      ),
                xlab = paste("Module membership in", module, "module"),
                ylab = paste("Correlation with", trait),
                title = paste("Module membership vs. taxa significance\n")
                ) +
        geom_hline(yintercept=0)
      print(p)
      dev.off()

      #Export clinical trait taxa
      write.csv(file = paste0(module, "_", trait,"_TOPclinicaltaxa.csv"), top_taxa)

      if(cor.sign == "neg"){
        traitsig_df_export <- traitsig_df %>% filter(moduleColors == module) %>% arrange((eval(parse(text = colname))))
      }
      if(cor.sign == "pos"){
        traitsig_df_export <- traitsig_df %>% filter(moduleColors == module) %>% arrange(desc(eval(parse(text = colname))))
      }

      write.csv(file = paste0(module, "_", trait,"_ALLclinicaltaxa.csv"), traitsig_df_export)
    }

    #Assign key variables to global environment
    geneModuleMembership.global <<- geneModuleMembership
    geneTraitSignificance.global <<- geneTraitSignificance
    moduleTaxa.global <<- moduleTaxa

    #Export module-level tables
    setwd(tempdir_parent)
    write.csv(file = paste0(module, "_hubtaxa.csv"), hubtaxa)
    write.csv(file = paste0(module, "taxa_membership.csv"),
              member_df %>% arrange(desc(eval(parse(text = paste0("MM", module))))))

  }
}












# Function: Boostrapping for module stability -----------------------------------------------------
#Purpose: Low module size leads to more clusters that have different correlation and significance with clinical traits
#However their validity as stable clusters is unclear because of small module size.
#Here we test stability by bootstrapping our original dataset and testing whether the modules appear in the boostrapped sample sets with similar hub genes
#Code:
#Save original dataset
datExpr_parent <- datExpr
datTraits_parent <- datTraits
datMeta_parent <- datMeta

#Inspect dataframes
head(datExpr)
head(datTraits)
head(datMeta)


#=== Jaccard Similarity Function ===
#https://www.statology.org/jaccard-similarity-in-r/
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

#=== Module Membership Functions ===
module_membership <- function(data){
  # === Calculate Module Membership for each taxa: ===
  #Module Membership = correlation between Module eigengene and the taxa
  #reqiured to calculate hub taxa

  modNames = substring(names(MEs), 3) #Remove ME from module names

  geneModuleMembership = as.data.frame(cor(data, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")

  member_df <- cbind(geneModuleMembership, MMPvalue, moduleColors)

  #Export to upper-level environment
  return(member_df)
} #calculates module membership for ALL modules
taxa_in_module <- function(top.taxa.prop, module,
                           membership.data #output from module_membership
){
  #=== Hub taxa ===
  #Top 5% taxa for module membership
  top.taxa = (nrow(membership.data %>% filter(moduleColors == module)) * top.taxa.prop) %>% floor
  hubtaxa <- membership.data %>% filter(moduleColors == module) %>% arrange(desc(eval(parse(text = paste0("MM", module))))) %>%
    slice_max(order_by = eval(parse(text = paste0("MM", module))),n = top.taxa)

  #=== All taxa in module ===
  all_members <- membership.data %>% filter(moduleColors == module)

  export <- list(all_members, hubtaxa); names (export) <- c("all_members", "hubtaxa")
  return(export)
} #returns all taxa and hub taxa for SPECIFIC modules


#Goal: Bootstrap 1000 times

#=== Goal: Resample  and report: ===
#1) Module membership
#2) Cluster stability using Jaccard Similarity Index (doi:10.1016/j.csda.2006.11.025)

#Network inputs
module_stability <- function(data.input = datExpr_parent,
                             net,
                             input.minModulesize,
                             input.softpower,
                             nboot = 1000,
                             top_taxa_prop = 0.2,
                             export = TRUE
){

  #Original network
  print("Setup original network)")
  network_modules(net_type = net, minModuleSize = input.minModulesize,
                  softPower = input.softpower, export = FALSE, data = data.input)

  print("Original network - Module member of individual taxa")
  og_member_df <- module_membership(data = data.input)

  module_list <- unique(moduleColors)
  og.modules <- module_list
  print(og.modules)
  #Original Module Taxa
  for(i in 1:length(module_list)){
    module.name <- module_list[i]
    name <- paste0("original_", module.name)
    assign(paste(name), taxa_in_module(top.taxa.prop = 0.05, module = module.name,
                                       membership.data = og_member_df))
  }



  #Bootstrapping loop
  bootstrap.result.members <- data.frame()
  bootstrap.result.hub <- data.frame()
  for (k in 1:nboot){
    print(paste("Resample set", k))
    #Resample with replacement
    set.seed(k)
    resampled <- data.input[sample(nrow(data.input), replace = TRUE), ]
    network_modules(net_type = net, minModuleSize = input.minModulesize,
                    softPower = input.softpower, export = FALSE, data = resampled)
    boot_member_df <- module_membership(data = resampled)

    module_list <- unique(moduleColors)
    boot.modules <- module_list
    #Loop to get taxa list from bootstrapped modules
    for(i in 1:length(module_list)){
      module.name <- module_list[i]
      name <- paste0("bootstrap_", module.name)
      assign(paste(name), taxa_in_module(top.taxa.prop = top_taxa_prop, module = module.name,
                                         membership.data = boot_member_df))
    }

    #=== Module similarity matrix (all bootstrap vs. all original) ===
    #== All taxa ==
    jaccard.matrix <- data.frame()
    for(i in og.modules){ #get list of taxa in original module
      og.module <- i #module name/colour
      og <- eval(parse(text = paste0("original_", i, "$all_members")))
      for(j in boot.modules){ #get list of taxa in bootstrapped module. Compare boostrapped v original via Jaccard similarity.
        boot <- eval(parse(text = paste0("bootstrap_", j, "$all_members")))
        j.stat <- jaccard(og %>% rownames(), boot %>% rownames())
        #Export similarity of bootstrap module with original module
        bootstrap.module <- j
        row <- cbind(og.module, bootstrap.module, j.stat)
        jaccard.matrix <- rbind(jaccard.matrix, row)
      }
    }

    # === Record  bootstrap cluster jaccard statistic ===
    #Goal:
    #Assign the most similar bootstrapped module to each original module for this itration

    # == Modules ==
    #Code:
    #All taxa
    head(jaccard.matrix %>% arrange(desc(j.stat)))
    dumdf <- jaccard.matrix %>% group_by(og.module) %>% slice(which.max(j.stat)) %>% mutate(bootstrap = k) #assigns a bootstrapped module to each original module
    bootstrap.result.members <- rbind(bootstrap.result.members, dumdf) #add assignments to results dataframe


    #==Hub taxa  ==
    #= * WARNING * this feautre is unvalidated and does not output meaningful results =
    jaccard.hub.matrix <- data.frame()
    for(i in og.modules){ #get list of taxa in original module
      og.module <- i #og module colour
      #Find bootstrap module colour assigned to og module colour
      df1 <- dumdf %>% filter(og.module %in% i) %>% select(bootstrap.module)
      bootstrap.module <- df1$bootstrap.module

      #Bootstrap and original hub taxa
      og <- eval(parse(text = paste0("original_", i, "$hubtaxa")))
      boot <- eval(parse(text = paste0("bootstrap_", bootstrap.module, "$hubtaxa")))

      #Intersect og and boot
      #OG hub taxa is top 5% of taxa
      #Boot hub taxa is top 20% of taxa (default) or user defined proportion (top_taxa_prop)
      #Intersect to find whether top 5% hub taxa in original dataset are conserved in boostrapped top 20% hub taxa for the module
      boot.int <- intersect(rownames(og), rownames(boot))


      #Compare hub taxa (Jaccard Similarity)
      j.stat <- jaccard(og %>% rownames(), boot.int)

      #Export similarity of bootstrap module with original module
      row <- cbind(og.module, bootstrap.module, j.stat)
      jaccard.hub.matrix <- rbind(jaccard.hub.matrix, row)}

    # === Record  bootstrap cluster jaccard statistic ===
    #Hub taxa
    head(jaccard.hub.matrix %>% arrange(desc(j.stat)))
    dumdf <- jaccard.hub.matrix %>% group_by(og.module) %>% slice(which.max(j.stat)) %>% mutate(bootstrap = k)
    bootstrap.result.hub <- rbind(bootstrap.result.hub, dumdf)
  }

  member.summary <- bootstrap.result.members %>%
    group_by(og.module) %>% mutate(j.stat = as.numeric(j.stat)) %>% summarise(mean = mean(j.stat))

  hub.summary <- bootstrap.result.hub %>%
    group_by(og.module) %>% mutate(j.stat = as.numeric(j.stat)) %>% summarise(mean = mean(j.stat))

  #Export CSV
  if(export == TRUE){
    wd <- paste(parent.dir, net, paste0("minmod",input.minModulesize), sep = "/")
    setwd(wd)
    write.csv(member.summary, "ValidateModuleMembers_summary.csv")
    #write.csv(hub.summary, "ValidateHubTaxa_summary.csv") #removed from exports because it is unvalidated
    write.csv(bootstrap.result.members, "ValidateModuleMembers.csv")
    write.csv(bootstrap.result.hub, "ValidateHubTaxa.csv")
  }

  return(list(member.summary, hub.summary, bootstrap.module))
}


# Function - Visualize Network - cytoscape -----------------------------------------------
cyto_fun <- function(module,
                     r = 0.4,
                     netname){
  print(cytoscapePing ())
  print(cytoscapeVersionInfo ())
  #Setup datasets (univeral for any module)
  taxa_modkey <- adjacency %>% rownames %>% as.data.frame %>% cbind(moduleColors) #taxa - module key
  cor.matrix <- rcorr(datExpr %>% as.matrix(), type = "pearson")

  #Module network visualization
  taxa_mod <- taxa_modkey %>% filter(moduleColors == module) %>% rename (taxa = ".")
  cor.module <- cor.matrix$r %>% data.frame %>%
    select(taxa_mod$taxa) %>%
    rownames_to_column(var = "taxa") %>% filter(taxa %in% taxa_mod$taxa) %>% column_to_rownames(var = "taxa")

  #== Prune correlations ==
  #Remove correlations below threshold (e.g r > 0.2)
  cor.module[cor.module < r] <- 0

  #Remove taxa that don't have any correlations after hard-thresholding
  rmtaxa <- rowSums(cor.module) %>% data.frame() %>% rename (corsums = ".") %>% filter(corsums == 1) %>% rownames()
  cor.module <- cor.module %>% select(-rmtaxa) %>%
    rownames_to_column(var = "taxa") %>% filter(!(taxa %in% rmtaxa)) %>% column_to_rownames(var = "taxa")

  #== Cytoscape visualization ==
  net <- igraph::graph_from_adjacency_matrix(cor.module %>% as.matrix(), weighted=T, mode="undirected", diag=F)
  createNetworkFromIgraph(net, paste(module, netname, r, sep = "_"))
}









# Run signed network ------------------------------------------------------
select_network <- FALSE
if(select_network == TRUE){
net <- "signed"
network_selection(net_type = net)
softpower <- 3
visibility <- (-10^(-12))
#=== Differing module size ===
#size = 10
network_modules(net_type = net, minModuleSize = 10, softPower = softpower, network_visibility = visibility)
module_traits()
investigate_module(clin.trait = c("phq9_total", "gad7_total", "age", "BMI"), mod = c("green", "brown", "yellow"))
module_stability(net = "signed", input.minModulesize = 10, input.softpower = softpower)

#size = 15
network_modules(net_type = net, minModuleSize = 15, softPower = softpower, network_visibility = visibility)
module_traits()
investigate_module(clin.trait = c("phq9_total", "gad7_total", "age", "BMI"), mod = c("blue", "yellow", "brown"))
module_stability(net = "signed", input.minModulesize = 15, input.softpower = softpower)

#size = 50
network_modules(net_type = net, minModuleSize = 50, softPower = softpower, network_visibility = visibility)
module_traits()
investigate_module(clin.trait = c("phq9_total", "gad7_total", "age", "BMI"), mod = "blue")
module_stability(net = "signed", input.minModulesize = 50, input.softpower = softpower)
}


# Run select network - minMod = 50  ---------------------------------------
getwd()
parent.dir <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/"
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod50")
visibility <- (-10^(-12))
datExpr <- clr %>% select(-metadata) %>%
  as.matrix() %>% as.data.frame() %>% column_to_rownames(var = "Sample") %>% mutate_all(as.numeric)

datTraits <- clr %>% select(numeric_trait) %>%
  as.matrix() %>% as.data.frame() %>% column_to_rownames(var = "Sample") %>% mutate_all(as.numeric)
traitColors <- WGCNA::numbers2colors(datTraits, signed = FALSE)

datMeta <- clr %>% select(metadata) %>%
  as.matrix() %>% as.data.frame() %>% column_to_rownames(var = "Sample")


network_selection(net_type = "signed")
network_modules(net_type = "signed", minModuleSize = 50, softPower = 3, network_visibility = (-10^(-12)))
module_traits()
investigate_module(clin.trait = c("phq9_total", "gad7_total", "age", "BMI"), mod = c("brown", "blue"))

cyto_fun(module = "brown", netname = "signed_mindmod50", r = 0.4)
cyto_fun(module = "brown", netname = "signed_mindmod50", r = 0.2)

setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod50")
write.csv(MEs, "eigentaxa.csv")


# Run Select Network - minMod = 25 ----------------------------------------
getwd()
parent.dir <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/"

softpower <- 3
visibility <- (-10^(-12))
datExpr <- clr %>% select(-metadata) %>%
  as.matrix() %>% as.data.frame() %>% column_to_rownames(var = "Sample") %>% mutate_all(as.numeric)

datTraits <- clr %>% select(numeric_trait) %>%
  as.matrix() %>% as.data.frame() %>% column_to_rownames(var = "Sample") %>% mutate_all(as.numeric)
traitColors <- WGCNA::numbers2colors(datTraits, signed = FALSE)

datMeta <- clr %>% select(metadata) %>%
  as.matrix() %>% as.data.frame() %>% column_to_rownames(var = "Sample")


network_selection(net_type = "signed")
network_modules(net_type = "signed", minModuleSize = 25, softPower = 3, network_visibility = (-10^(-12)))
module_traits()
investigate_module(clin.trait = c("phq9_total", "gad7_total", "age", "BMI"), mod = c("brown", "blue"))

setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod25")
write.csv(MEs, "eigentaxa.csv")

#Run stability 
module_stability(net = "signed", input.minModulesize = 25, input.softpower = softpower)



# Rule-out age correlation with gad/phq --------------------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna")
#Age
cor.test(datTraits$gad7_total, datTraits$age, method = "pearson")
cor.test(datTraits$phq9_total, datTraits$age, method = "pearson")

#GAD/PHQ9
cor.test(datTraits$gad7_total, datTraits$phq9_total, method = "pearson")



#Export
sink(file = "variable_correlation.txt")
print("GAD7 vs. age")
cor.test(datTraits$gad7_total, datTraits$age, method = "pearson")
print("PHQ9 vs. age")
cor.test(datTraits$phq9_total, datTraits$age, method = "pearson")
print("GAD7 vs. PHQ9")
cor.test(datTraits$gad7_total, datTraits$phq9_total, method = "pearson")
print("GAD7 vs. DARS")
cor.test(datTraits$gad7_total, datTraits$dars17_total, method = "pearson")
print("BMI vs. age")
cor.test(datTraits$BMI, datTraits$age, method = "pearson")
print("BMI vs. gad7")
cor.test(datTraits$BMI, datTraits$gad7_total, method = "pearson")
print("BMI vs. phq9")
cor.test(datTraits$BMI, datTraits$phq9_total, method = "pearson")
sink()




# Functions: Top Taxa plots ----------------------------------------------------
important_taxa_scatter <- function(df, #dataframe containing taxa and clinical traits (character)
                                   trait, #name of clinical trait (character)
                                   taxa, #col 1 = taxa names as a character list; col 2 = clinical correlation (numeric). Taxonomic format must be the same as df column names
                                   module, #name (character) for export file name
                                   export_params = c(4000, 4000, 300, 1), #width, height, res, ncol
                                   textsize = 10 #numeric
){
  #Abbreviate taxa names for axis label
  taxa.names <- taxa[1]
  taxa <- taxa %>% mutate(label = sp_genus(dat_glom, taxa.names))

  #Scatter plot loop
  plot.df <- data.frame()
  plotlist <- list()
  for (i in 1:nrow(taxa)){

    #ggpubr plot input parameters
    data <- df
    y.name <- eval(parse(text = paste0("taxa[", i, ",1]")))
    x <- trait

    ylabel <- eval(parse(text = paste0("taxa[", i, ",3]")))
    x_rcoef <- max(eval(parse(text = paste(data, x, sep = "$"))), na.rm=T) * 0.25
    y_rcoef <- max(eval(parse(text = paste(data, y.name, sep = "$"))), na.rm=T) * 1.2


    print(x_rcoef)
    print(y_rcoef)
    #scatter plot
    p <- ggscatter(x = paste(x),
                   y = paste(y.name),
                   data = get(data),
                   add = "reg.line",  # Add regression line
                   conf.int = TRUE, # Add confidence interval
                   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                   cor.coeff.args = list(method = "pearson",
                                         label.x = x_rcoef,
                                         label.y = y_rcoef,
                                         p.accuracy = 0.001
                                         ),
                   cor.coef.size = 3,
                   add.params = list(color = "black", fill = "lightgray"),
                   ylab = ylabel,
                   color = "white",
                   fill = "white") %>% ggadd("point", alpha = 0.4, size = 1.8, color = "gray26") %>% ggpar(font.y = (textsize), font.x = textsize)

    #plot list for loop
    plotname <- paste("scatter", i, sep  = ".")
    assign(paste(plotname), p)
    plot.df <- rbind(plot.df, plotname)
    plotlist[[i]] <- p

  }

  p <- cowplot::plot_grid(plotlist = plotlist, scale = 0.8, ncol = export_params[4])
  print(p)
  tiff(paste0("scatter_wcna_", module, "_", trait, ".tiff"), width = export_params[1], height = export_params[2], res = export_params[3])
  print(p)
  dev.off()

}

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

gdens <- function(df, #name of dataframe
                  x, #variable to plot
                  tsize = 24,
                  colour_hex
                  )
{ #input(object, text, text, colour/hex)
  ggplot(df, aes_string(x = x)) +
    geom_density(alpha=0.5, fill = colour_hex, color = colour_hex) +
    theme_minimal_hgrid(12) +
    theme(text = element_text(size=tsize),
          plot.margin = unit(c(0,1,0,0), "cm"))

}



hub_taxa_density <- function(df, #dataframe containing taxa and clinical traits (character)
                             taxa, #col 1 = taxa names as a character list; Taxonomic format must be the same as df column names
                             module, #name (character) for export file name
                             export_params = c(4000, 4000, 300, 1), #width, height, res, ncol
                             textsize = 10,
                             hex
                             )
  {
  #Abbreviate taxa names for axis label
  taxa.names <- taxa[1]
  taxa <- taxa %>% mutate(label = sp_genus(dat_glom, taxa.names))

  #Scatter plot loop
  plot.df <- data.frame()
  plotlist <- list()
  for (i in 1:nrow(taxa)){
    #ggpubr plot input parameters
    data <- df
    x.name <- eval(parse(text = paste0("taxa[", i, ",1]")))
    xlabel <- eval(parse(text = paste0("taxa[", i, ",2]")))

    #scatter plot
    p <- gdens(x = paste(x.name),
               df = data,
               tsize = textsize,
               colour_hex = hex[i])
    print(p)

    #plot list for loop
    plotname <- paste("plot", i, sep  = ".")
    assign(paste(plotname), p)
    plot.df <- rbind(plot.df, plotname)
    plotlist[[i]] <- p

  }

  p <- cowplot::plot_grid(plotlist = plotlist, scale = 0.8, ncol = export_params[4])

  print(p)
  tiff(paste0("density_wcna_", module, ".tiff"), width = export_params[1], height = export_params[2], res = export_params[3])
  print(p)
  dev.off()

}



# Hub taxa colours (hex codes) ------------------------------------------------------------
hubhex <- c("#6DA34D", "#D74E09", "#0E6BA8", "#933982")

# Top taxa - brown module -----------------------------------------------------------
#Import clincally important taxa list
age.important <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod50/brown/age/brown_age_TOPclinicaltaxa.csv") %>% data.frame %>% rename(taxa = X)
gad.important <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod15/brown/gad7_total/brown_gad7_total_TOPclinicaltaxa.csv") %>% data.frame %>% rename(taxa = X)
phq.important <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod50/brown/phq9_total/brown_phq9_total_TOPclinicaltaxa.csv") %>% data.frame %>% rename(taxa = X)

#Import brown module hub taxa
brown.hub <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod50/brown/brown_hubtaxa.csv") %>% data.frame %>% rename(taxa = X)

#CLR dataset
clr[1:10, 1:20] #df

#Set wd for exports
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod50/brown")


#Plot hub taxa (density plots)
hub_taxa_density(df = rel, taxa = brown.hub %>% select(taxa),
                 module = "brown_rel", export_params = c(7000, #width
                                                     1000, #height
                                                     300, #res
                                                     4), #ncol

                 textsize = 10,
                 hex = hubhex
                 )

hub_taxa_density(df = clr, taxa = brown.hub %>% select(taxa),
                 module = "brown", export_params = c(2300, #width
                                                         1600, #height
                                                         300, #res
                                                         2), #ncol
                 hex = hubhex,
                 textsize = 10
)


# Import stability data ---------------------------------------------------------
stabimport <- function(df){
  df <- df %>% select(-X) %>% rename(module = og.module)
  return(df)
}
mod10.stability <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod10/ValidateModuleMembers_summary.csv") %>% stabimport
mod15.stability <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod15/ValidateModuleMembers_summary.csv") %>% stabimport
mod50.stability <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod50/ValidateModuleMembers_summary.csv") %>% stabimport
mod25.stability <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod25/ValidateModuleMembers_summary.csv") %>% stabimport

# Function: Stability plots ---------------------------------------------------------



stabplot <- function(data, #stability dataframe
                     alpha = 1 #geom_rectangle alpha (background)
                     )
  {
ggdotchart(data, x = "module", y = "mean",
                sorting = "ascending",
                ggtheme = theme_pubr(),
                add = "segment",
                add.params=list(color = "black"),
                ylab = "Module Stability",
                xlab = "",
                dot.size = 4,
                shape = 21) +
  ylim(0,1) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
  geom_point(aes(colour = module), size = 3.5) + scale_color_manual(values = palette) + theme(legend.position = "none") +
  geom_rect(ymin = 0, ymax = 0.5, xmin = -Inf, xmax = Inf,
            fill = "red", alpha = 0.02*alpha) +
  geom_hline(yintercept = 0.5, linetype = "dotted", size = 0.5) +
  geom_rect(ymin = 0.5, ymax = 0.6, xmin = -Inf, xmax = Inf,
            fill = "orange", alpha = 0.02*alpha) +
  geom_hline(yintercept = 0.6, linetype = "dotted", size = 0.5) +
  geom_rect(ymin = 0.6, ymax = 0.75, xmin = -Inf, xmax = Inf,
            fill = "palegreen", alpha = 0.02*alpha) +
  geom_hline(yintercept = 0.75, linetype = "dotted", size = 0.5) +
  geom_rect(ymin = 0.75, ymax = 1, xmin = -Inf, xmax = Inf,
            fill = "green", alpha = 0.03*alpha)
  }


# Export Stability Plots --------------------------------------------------

#Set wd for exports
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed")

#Export plots
w = 1200
h = 1200
res = 300

tiff("mod10_stability.tiff", width = w, height = h, res = res)
palette <- mod10.stability %>% arrange(mean) %>% pull(module) #in order of most to least stable
mod10.stability %>% stabplot(alpha = 1.2)
dev.off()

tiff("mod15_stability.tiff", width = w, height = h, res = res)
palette <- mod15.stability %>% arrange(mean) %>% pull(module) #in order of most to least stable
mod15.stability %>% stabplot(alpha = 1.5)
dev.off()

tiff("mod50_stability.tiff", width = w, height = h, res = res)
palette <- mod50.stability %>% arrange(mean) %>% pull(module) #in order of most to least stable
mod50.stability %>% stabplot(alpha = 2)
dev.off()

tiff("mod25_stability.tiff", width = w, height = h, res = res)
palette <- mod25.stability %>% arrange(mean) %>% pull(module) #in order of most to least stable
mod25.stability %>% stabplot(alpha = 2)
dev.off()


# Taxa significance plots - Edit - Highlight Hub Taxa --------------------------------------------

#Import brown module hub taxa
brown.hub <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod50/brown/brown_hubtaxa.csv") %>% data.frame %>% rename(taxa = X)

#Import minmod 50 module membership
mod50_mm <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed/minmod50/brown/browntaxa_membership.csv") %>% rename(taxa = X)


#Ensure required df in global environment
head(datExpr)
head(datTraits)
head(brown.hub)




#Plot  function
scatter_taxasig_hub <- function(taxa_mm, #taxa module membership dataframe
                            module, #Moduel name (character)
                            trait, #Clinical trait (character; must match name in datTraits)
                            cor.sign, #"pos" or "neg"
                            hex = hubhex, #hubhex variable stored in - Functions: Top Taxa plots - section
                            hubtax,
                            scale.label.y = 0.8, #how far along y-axis R-squared and p-value reported.
                            scale.label.x = 0.8 #how far along x-axis R-squared and p-value reported.
                            )
  #Plots module membership vs. taxa clinical significance (i.e corrleation with trait)
  #Adds colour highlights to hub taxa
  {
  #Get clinical significance dataframe
  voi <- as.data.frame(eval(parse(text = paste0 ("datTraits$", trait))))
  taxa_clinsig = as.data.frame(cor(datExpr, voi, use = "p")) %>% rename(clinical_correlation = `eval(parse(text = paste0("datTraits$", trait)))`) %>%
    rownames_to_column(var = "taxa")
  print(taxa_clinsig)
  #Setup plotting dataframe
  modNames = substring(names(MEs), 3) #Remove ME from module names
  column = match(module, modNames) #Locate column belonging to module

  #df <- data.frame(abs(taxa_mm[moduleTaxa.global, column]), taxa_clinsig[moduleTaxa.global, 1]) %>% mutate_all(as.numeric) %>% mutate(taxa = taxa_mm[moduleTaxa.global,] %>% rownames)
  df <- taxa_mm %>% filter(moduleColors == paste(module)) %>% select(taxa, paste0("MM",module)) %>% left_join(taxa_clinsig, by = "taxa")
  colnames(df) <- c("taxa", "taxa_mod_membership", "clinical_correlation")
  #Setup point highlight dataframe
  df_highlight <- df %>% filter(taxa %in% hubtax)


  #GGplot
  p <- ggscatter(x = "taxa_mod_membership",
                 y = "clinical_correlation",
                 data = df,
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE, # Add confidence interval
                 cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                 cor.coeff.args = list(method = "pearson", label.sep = "\n",
                                       label.x = max(df$taxa_mod_membership)* scale.label.x,
                                       label.y = ifelse(cor.sign == "pos",
                                                        min(df$clinical_correlation) * scale.label.y,
                                                        max(df$clinical_correlation) * scale.label.y)
                 ),
                 xlab = paste("Module membership in", module, "module"),
                 ylab = paste("Correlation with", trait),
                 title = " ") +
    theme(legend.position = "right",
          legend.title = element_blank()) +
    geom_hline(yintercept=0) +
    geom_point(data = df_highlight,
               aes(x = taxa_mod_membership, y = clinical_correlation, color = taxa, fill = taxa),
               size = 2.5, shape = 16) +
    scale_color_manual(values=c(hex)) +
    scale_fill_manual(values = c(hex))
  print(p)
}


#Set wd
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/signed")

#Set hubhex
hubhex_reorder <- c("#6DA34D", "#933982","#D74E09", "#0E6BA8")

#Plots
p1 <- scatter_taxasig_hub(module = "brown", trait = "gad7_total", cor.sign = "neg",
                hubtax = brown.hub %>% pull(taxa),
                hex = hubhex_reorder,
                taxa_mm = mod50_mm,
                scale.label.y = 0.7,
                scale.label.x = 0.7)
p2 <- scatter_taxasig_hub(module = "brown", trait = "phq9", cor.sign = "neg",
                          hubtax = brown.hub %>% pull(taxa),
                          hex = hubhex_reorder,
                          taxa_mm = mod50_mm,
                          scale.label.y = 0.7,
                          scale.label.x = 0.7)
p3 <- scatter_taxasig_hub(module = "brown", trait = "age", cor.sign = "pos",
                          hubtax = brown.hub %>% pull(taxa),
                          hex = hubhex_reorder,
                          taxa_mm = mod50_mm,
                          scale.label.x = 0.7)
p4 <- scatter_taxasig_hub(module = "brown", trait = "BMI", cor.sign = "neg",
                          hubtax = brown.hub %>% pull(taxa),
                          hex = hubhex_reorder,
                          taxa_mm = mod50_mm,
                          scale.label.x = 0.7,
                          scale.label.y = 0.7)


tiff("brown_taxasig.tiff", width = 2400, height = 4000, res = 300)
cowplot::plot_grid(p1 + theme(legend.position = "none") + expand_limits(x = 0.85),
                   p2 + theme(legend.position = "none") + expand_limits(y = 0.05) + expand_limits(x = 0.85),
                   p3 + theme(legend.position = "none") + expand_limits(x = 0.85),
                   p4 + theme(legend.position = "none") + expand_limits(x = 0.85),
                   get_legend(p1) %>% as_ggplot(),
                   scale = 0.9, ncol = 2)
dev.off()



# Run Age and BMI corrected Network ------------------------------------------------------
# === Age and BMI adjusted data ===
#set wd
parent.dir <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected"

setwd(parent.dir)
#Remove sample with age = NA
covariate.filter <- function(df){
  out <- df %>% filter(!(is.na(age)) & !(is.na(BMI)))
  return(out)
}

datExpr_age <- clr %>% covariate.filter %>% select(-metadata) %>%
  as.matrix() %>% as.data.frame() %>% column_to_rownames(var = "Sample") %>% mutate_all(as.numeric)

datTraits_age <- clr %>% select(numeric_trait) %>% covariate.filter %>%
  as.matrix() %>% as.data.frame() %>% column_to_rownames(var = "Sample") %>% mutate_all(as.numeric)

#number of samples
  #n = 163; 22aug22
datExpr_age %>% nrow
datTraits_age %>% nrow


#qc
(rownames(datExpr_age) == rownames(datTraits_age)) %>% all #qc data is aligned

#covariate correction
datAdj <- empiricalBayesLM(data = datExpr_age,
                           removedCovariates = datTraits_age %>% select(age, BMI), verbose = 10)

#Input data for WCNA
datExpr <- datAdj$adjustedData
datTraits <- datTraits_age %>% select(-c(age, BMI))
traitColors <- WGCNA::numbers2colors(datTraits, signed = FALSE)


# === Run WCNA ===
run_covar_nonselect <- FALSE
if(run_covar_nonselect == TRUE){
  net <- "signed"
  network_selection(net_type = net)

  softpower <- 3

  network_modules(net_type = net, minModuleSize = 10, softPower = softpower)
  module_traits()
  investigate_module(clin.trait = c("phq9_total", "gad7_total"), mod = c("blue", "yellow", "brown", "red", "turquoise"))
  module_stability(net = "signed", input.minModulesize = 10, input.softpower = softpower)


  network_modules(net_type = net, minModuleSize = 15, softPower = softpower)
  module_traits()
  investigate_module(clin.trait = c("phq9_total", "gad7_total"), mod = c("brown", "yellow"))
  module_stability(net = "signed", input.minModulesize = 15, input.softpower = softpower)
}
 # == Select Model ==
run_covar_select <- TRUE
run.mod.stability <- FALSE
run.cytoscape <- FALSE
if(run_covar_select == TRUE){
  net <- "signed"
  network_selection(net_type = net)
  softpower <- 3
  network_modules(net_type = net, minModuleSize = 25, softPower = softpower)
  module_traits()
  investigate_module(clin.trait = c("phq9_total", "gad7_total", "dars17_total"), mod = c("brown"))
if(run.mod.stability == TRUE){
    module_stability(net = "signed", input.minModulesize = 25, input.softpower = softpower)
  }
}

# === Export Module Eigengenes === 
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25")
write.csv(MEs, "eigentaxa.csv")



  # = Export network vis to cytoscape =
run_covar_cytoscape = FALSE
if(run_covar_cytoscape == TRUE){
cyto_fun(module = "brown", netname = "signed_adjusted_mindmod25", r = 0.4)
cyto_fun(module = "brown", netname = "signed_adjusted_mindmod25", r = 0.35)
cyto_fun(module = "brown", netname = "signed_adjusted_mindmod25", r = 0.2)
}

# === Plot stability graphs ===
setwd(parent.dir)

w = 1200
h = 1200
res = 300

mod10.stability <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod10/ValidateModuleMembers_summary.csv") %>% rename(module = og.module)
mod15.stability <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod15/ValidateModuleMembers_summary.csv") %>% rename(module = og.module)
mod25.stability <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/ValidateModuleMembers_summary.csv") %>% rename(module = og.module)


tiff("mod10_stability.tiff", width = w, height = h, res = res)
palette <- mod10.stability %>% arrange(mean) %>% pull(module)
mod10.stability %>% stabplot(alpha = 1.5)
dev.off()

tiff("mod15_stability.tiff", width = w, height = h, res = res)
palette <- mod15.stability %>% arrange(mean) %>% pull(module)
mod15.stability %>% stabplot(alpha = 2)
dev.off()

tiff("mod25_stability.tiff", width = w, height = h, res = res)
palette <- mod25.stability %>% arrange(mean) %>% pull(module)
mod25.stability %>% stabplot(alpha = 3)
dev.off()




# Top taxa - Adjusted - Brown -----------------------------------------------------------
#Import clincally important taxa list
gad.important <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/brown/gad7_total/brown_gad7_total_TOPclinicaltaxa.csv") %>% data.frame %>% rename(taxa = X)
phq.important <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/brown/phq9_total/brown_phq9_total_TOPclinicaltaxa.csv") %>% data.frame %>% rename(taxa = X)
dars.important <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/brown/dars17_total/brown_dars17_total_TOPclinicaltaxa.csv") %>% data.frame %>% rename(taxa = X)


#Import brown module hub taxa
brown.hub <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/brown/brown_hubtaxa.csv") %>% data.frame %>% rename(taxa = X)

#Dataset; adjusted CLR taxa abundance and clinical trait
clr_adj <- clr %>% select(Sample, dars17_total, gad7_total, phq9_total) %>% left_join(datAdj$adjustedData %>% data.frame %>% rownames_to_column(var = "Sample"), by = "Sample")
head(clr_adj)

#Set wd for exports
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/brown")

#Plot hub taxa (density plots)
hub_taxa_density(df = clr_adj, taxa = brown.hub %>% select(taxa),
                 module = "brown", export_params = c(4000, #width
                                                     1000, #height
                                                     300, #res
                                                     3), #ncol
                 hex = hubhex,
                 textsize = 11
)




# Taxa significance plots - Adjusted Brown--------------------------------------------
#Import brown module hub taxa
brown.hub <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/brown/brown_hubtaxa.csv") %>% data.frame %>% rename(taxa = X)


#Import minmod 25 module membership
mod50_mm <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected/signed/minmod25/brown/browntaxa_membership.csv") %>% rename(taxa = X)


#Ensure required df in global environment
head(datExpr)
head(datTraits)
head(brown.hub)




#Plot  function
scatter_taxasig_hub <- function(taxa_mm, #taxa module membership dataframe
                                module, #Moduel name (character)
                                trait, #Clinical trait (character; must match name in datTraits)
                                cor.sign, #"pos" or "neg"
                                hex = hubhex, #hubhex variable stored in - Functions: Top Taxa plots - section
                                hubtax,
                                scale.label.y = 0.8, #how far along y-axis R-squared and p-value reported.
                                scale.label.x = 0.8 #how far along x-axis R-squared and p-value reported.
)
  #Plots module membership vs. taxa clinical significance (i.e corrleation with trait)
  #Adds colour highlights to hub taxa
{
  #Get clinical significance dataframe
  voi <- as.data.frame(eval(parse(text = paste0 ("datTraits$", trait))))
  taxa_clinsig = as.data.frame(cor(datExpr, voi, use = "p")) %>% rename(clinical_correlation = `eval(parse(text = paste0("datTraits$", trait)))`) %>%
    rownames_to_column(var = "taxa")
  print(taxa_clinsig)
  #Setup plotting dataframe
  modNames = substring(names(MEs), 3) #Remove ME from module names
  column = match(module, modNames) #Locate column belonging to module

  #df <- data.frame(abs(taxa_mm[moduleTaxa.global, column]), taxa_clinsig[moduleTaxa.global, 1]) %>% mutate_all(as.numeric) %>% mutate(taxa = taxa_mm[moduleTaxa.global,] %>% rownames)
  df <- taxa_mm %>% filter(moduleColors == paste(module)) %>% select(taxa, paste0("MM",module)) %>% left_join(taxa_clinsig, by = "taxa")
  colnames(df) <- c("taxa", "taxa_mod_membership", "clinical_correlation")
  #Setup point highlight dataframe
  df_highlight <- df %>% filter(taxa %in% hubtax)


  #GGplot
  p <- ggscatter(x = "taxa_mod_membership",
                 y = "clinical_correlation",
                 data = df,
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE, # Add confidence interval
                 cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                 cor.coeff.args = list(method = "pearson", label.sep = "\n",
                                       label.x = max(df$taxa_mod_membership)* scale.label.x,
                                       label.y = ifelse(cor.sign == "pos",
                                                        min(df$clinical_correlation) * scale.label.y,
                                                        max(df$clinical_correlation) * scale.label.y)
                 ),
                 xlab = paste("Module membership in", module, "module"),
                 ylab = paste("Correlation with", trait),
                 title = " ") +
    theme(legend.position = "right",
          legend.title = element_blank()) +
    geom_hline(yintercept=0) +
    geom_point(data = df_highlight,
               aes(x = taxa_mod_membership, y = clinical_correlation, color = taxa, fill = taxa),
               size = 2.5, shape = 16) +
    scale_color_manual(values=c(hex)) +
    scale_fill_manual(values = c(hex))
  print(p)
}


#Set wd
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV5/wcna/covariate-corrected")

#Set hubhex
hubhex_reorder <- c("#6DA34D", "#D74E09", "#0E6BA8")

#Plots
p1 <- scatter_taxasig_hub(module = "brown", trait = "gad7", cor.sign = "neg",
                          hubtax = brown.hub %>% pull(taxa),
                          hex = hubhex_reorder,
                          taxa_mm = mod50_mm,
                          scale.label.y = 0.7,
                          scale.label.x = 0.7)
p2 <- scatter_taxasig_hub(module = "brown", trait = "phq9", cor.sign = "neg",
                          hubtax = brown.hub %>% pull(taxa),
                          hex = hubhex_reorder,
                          taxa_mm = mod50_mm,
                          scale.label.y = 0.7,
                          scale.label.x = 0.7)
p3 <- scatter_taxasig_hub(module = "brown", trait = "dars17", cor.sign = "pos",
                          hubtax = brown.hub %>% pull(taxa),
                          hex = hubhex_reorder,
                          taxa_mm = mod50_mm,
                          scale.label.x = 0.7)


tiff("brown_taxasig.tiff", width = 2800, height = 3000, res = 300)
cowplot::plot_grid(p1 + theme(legend.position = "none") + expand_limits(x = 0.85),
                   p2 + theme(legend.position = "none") + expand_limits(y = 0.05) + expand_limits(x = 0.85),
                   p3 + theme(legend.position = "none") + expand_limits(x = 0.85),
                   get_legend(p1) %>% as_ggplot(),
                   scale = 0.85, ncol = 2)
dev.off()





