#Statistical tools        Primarily PERMANOVA, alpha diversity and the CLR transformation.
library(vegan)            #install.packages("vegan")
library(iNEXT)            #install.packages("iNEXT")
library(Tjazi)            #devtools::install_github("thomazbastiaanssen/Tjazi")
library(lme4)
library(lmerTest)

#Data Wrangling
library(tidyverse)        #install.packages("tidyverse")
library(knitr)            #install.packages("knitr")
library(waldo)            #install.packages("waldo")

#Plotting
library(ggplot2)          #install.packages("ggplot2")
library(ggforce)          #install.packages("ggforce")
library(patchwork)        #install.packages("patchwork")
library(ggbeeswarm)       #install.packages("ggbeeswarm")
library(metafolio)        #install.packages("metafolio")

#Disable strings automatically being read in as factors to avoid unintuitive behaviour.
options(stringsAsFactors = F)

#Set contrasts for type-3 sum of squares 
options(contrasts = c("contr.treatment","contr.poly"))


#Set a seed for the purposes of reproducibility in this document.
set.seed(1)


#Load in the genus level count table and the metadata file. 
counts   <- read.delim("raw/genus_table_from_dada2.csv", sep = ",", row.names = 1, header = T)
metadata <- read.delim("raw/metadata_DNBS.csv", sep = ",") %>% 
  mutate(Treatment = str_replace(Treatment, pattern = "x \\+ F", replacement = "x+F"))

#Fix metadata in accordance with genus names
metadata$ID <- gsub(metadata$ID, pattern = "-", replacement = ".")
metadata = metadata[!metadata$ID %in% c("DNBS.17", "DNBS.29", "DNBS.38", "DNBS.41", "DNBS.8"),]

waldo::compare(sort(colnames(counts)), sort(metadata$ID))

#metadata$Timepoint = factor(metadata$Timepoint, levels = c("FMT", "D0", "D7", "D42"))
metadata = metadata[metadata$Timepoint != "D7",]
metadata$Timepoint = factor(metadata$Timepoint, levels = c("FMT", "D0", "D32"))
#metadata$Timepoint = factor(metadata$Timepoint, levels = c("FMT", "D0", "D7", "D32"))
metadata = metadata[metadata$Experiment == "Pathological",]

counts = counts[,metadata$ID]

#Fork off your count data so that you always have an untouched version handy.
genus   <- counts

#make sure our count data is all numbers
genus   <- apply(genus,c(1,2),function(x) as.numeric(as.character(x)))

#Remove features with prevalence < 10% in two steps:
#First, determine how often every feature is absent in a sample
n_zeroes <- rowSums(genus == 0)

#Then, remove features that are absent in more than your threshold (90% in this case).
genus    <- genus[n_zeroes <= round(ncol(genus) * 0.90),]

#Perform a CLR transformation
genus.exp <- clr_c(genus)



# Alpha Diversity ---------------------------------------------------------
#Compute alpha diversity using a wrapper around the iNEXT library

alpha_diversity = get_asymptotic_alpha(species = counts, verbose = FALSE) 

#Add metadata for plotting and stats. Make sure the count table and metadata match up!
alpha_diversity$Timepoint  = metadata$Timepoint
alpha_diversity$ID         = metadata$animal_ID
alpha_diversity$Experiment = metadata$Experiment
alpha_diversity$type       = metadata$sample_type
alpha_diversity$Treatment  = metadata$Treatment


#Plot alpha diversity all at once using pipes
p_alpha <- alpha_diversity %>%
  
  #Wrangle the data to long format for easy plotting
  pivot_longer(!c(Timepoint, ID, Experiment, type, Treatment)) %>%

  #Pipe it all directly into ggplot2
  ggplot(aes(x     = Treatment,
             y     = value, 
             fill  = interaction(Treatment, Timepoint), 
             group = Treatment)) + 
  geom_boxplot(alpha = 1/2, coef = 100, position = position_dodge(), show.legend = FALSE) + 
  geom_beeswarm(size = 4, cex = 3, shape = 21, show.legend = FALSE) + 
 # scale_fill_manual(values = c("CTRL.FMT" = "#08306b", "abx+FMTCTR.D0" = "#9ecae1", "abx+FMTCTR.D7" = "#4292c6", "abx+FMTCTR.D32" = "#08519c", 
 #                              "DNBS.FMT" = "#67000d", "abx+FMTDNBS.D0" = "#fc9272", "abx+FMTDNBS.D7" = "#ef3b2c", "abx+FMTDNBS.D32" = "#a50f15"), 
 #                                        "Legend") +
  scale_fill_manual(values = c("CTRL.FMT" = "#08306b", "abx+FMTCTR.D0" = "#9ecae1", "abx+FMTCTR.D32" = "#08519c", 
                               "DNBS.FMT" = "#67000d", "abx+FMTDNBS.D0" = "#fc9272", "abx+FMTDNBS.D32" = "#a50f15"), 
                    "Legend") +

  guides(shape = "none", fill = guide_legend(override.aes = list(shape = c(23, 23, 21, 21, 21, 21)))) +
  ggh4x::facet_nested(name~Timepoint*type, scales = "free", strip = ggh4x::strip_nested(bleed = TRUE)) + theme_bw()  +

  ylab("Alpha diversity index") + xlab(NULL) + 
  theme(text = element_text(size = 14), axis.text.x = element_text(size = 11))


#Alpha Div stats Chao1
alpha_diversity %>%
  
  #Wrangle the data to long format for easy plotting
  pivot_longer(!c(Timepoint, ID, Experiment, type, Treatment)) %>%
  filter(name == "Chao1") %>%
  filter(type == "fecal sample") %>%
  filter(Timepoint != "D7") %>%
  lmerTest::lmer(value ~ Treatment * Timepoint + (1|ID), data = .) %>%
  car::Anova(., type = "III") %>% 
  
  capture.output(.,file = "stats_pathologic_arm/alpha_div_chao1.txt")



#Alpha Div stats Simpson Index
alpha_diversity %>%
  
  #Wrangle the data to long format for easy plotting
  pivot_longer(!c(Timepoint, ID, Experiment, type, Treatment)) %>%
  filter(name == "Simpson Index") %>%
  filter(type == "fecal sample") %>%
  filter(Timepoint != "D7") %>%
  lmerTest::lmer(value ~ Treatment * Timepoint + (1|ID), data = .)  %>%
  car::Anova(., type = "III") %>% 
  
  capture.output(.,file = "stats_pathologic_arm/alpha_div_simps.txt")


#Alpha Div stats shannon Entropy
alpha_diversity %>%
  
  #Wrangle the data to long format for easy plotting
  pivot_longer(!c(Timepoint, ID, Experiment, type, Treatment)) %>%
  filter(name == "Shannon Entropy") %>%
  filter(type == "fecal sample") %>%
  filter(Timepoint != "D7") %>%
  lmerTest::lmer(value ~ Treatment * Timepoint + (1|ID), data = .)%>%
  car::Anova(., type = "III") %>% 
  
  capture.output(.,file = "stats_pathologic_arm/alpha_div_shan.txt")


# Beta Diversity ----------------------------------------------------------


#Apply the base R principal component analysis function on our CLR-transformed data.
data.a.pca  <- prcomp(t(genus.exp))

#Extract the amount of variance the first four components explain for plotting. 
pc1 <- round(data.a.pca$sdev[1]^2/sum(data.a.pca$sdev^2),4) * 100
pc2 <- round(data.a.pca$sdev[2]^2/sum(data.a.pca$sdev^2),4) * 100
pc3 <- round(data.a.pca$sdev[3]^2/sum(data.a.pca$sdev^2),4) * 100
pc4 <- round(data.a.pca$sdev[4]^2/sum(data.a.pca$sdev^2),4) * 100

#Extract the scores for every sample for the first four components for plotting. 
pca  = data.frame(PC1 = data.a.pca$x[,1], 
                  PC2 = data.a.pca$x[,2], 
                  PC3 = data.a.pca$x[,3], 
                  PC4 = data.a.pca$x[,4])

#Add relevant information from the metadata
pca$ID                  = metadata$animal_ID
pca$Treatment           = metadata$Treatment
pca$`Sample Type`       = metadata$sample_type
pca$Timepoint           = metadata$Timepoint


#Compute euclidean distance over CLR-transformed values (i.e. Aitchison distance).
meta_perm  = metadata[metadata$Timepoint %in% c("D0", "D32"),]
genus_perm = genus.exp[,meta_perm$ID]
dis_ait = dist(t(genus_perm), method = "euclidean")
#Perform a PERMANOVA (PERmutational Multivariate ANalysis Of VAriance) test.
beta_div_path <- adonis2(dis_ait ~  Timepoint *Treatment,  
                         data = meta_perm, method = "euclidean", permutations = 10000, by = "terms")

beta_div_path %>%
  capture.output(.,file = "stats_pathologic_arm/beta_div_permanova.txt")


#First, the main plot. Plot the first two components of the PCA
p_beta <- ggplot(pca %>% mutate(type = "\nMicrobial composition over time (PCA)\n"), 
            aes(x       = PC1, 
                y       = PC2, 
                fill    = interaction(Treatment, Timepoint),
                colour  = Treatment,
                group   = ID, 
                shape = `Sample Type`, 
                alpha = `Sample Type`)) +  
  
  #Create the points and ellipses
  geom_path() +
  geom_point(size=5, col = "black", alpha = 1) + 
  stat_ellipse(geom = "polygon", alpha = 1/4, show.legend = FALSE) +
  annotate(geom = "text", 
           label = 
"Timepoint:
Treatment:" , 
           x = I(0.625), y = I(0.175), colour = "black",  size = 6, hjust = 0) +
  annotate(geom = "text", 
           label = 
"⁂
🞱" , 
           x = I(0.925), y = I(0.175), colour = "black", size = 6, hjust = 1/2) +
  
  
  #Adjust appearance
  guides(shape = "none", fill = guide_legend(override.aes = list(shape = c(23, 23, 21, 21, 21, 21)))) +
  scale_shape_manual(values = c("FMT Pool" = 23, "fecal sample" = 21)) +
  # scale_fill_manual(values = c("CTRL.FMT" = "#08306b", "abx+FMTCTR.D0" = "#9ecae1", "abx+FMTCTR.D7" = "#4292c6", "abx+FMTCTR.D32" = "#08519c", 
  #                              "DNBS.FMT" = "#67000d", "abx+FMTDNBS.D0" = "#fc9272", "abx+FMTDNBS.D7" = "#ef3b2c", "abx+FMTDNBS.D32" = "#a50f15"), 
  #                                        "Legend") +
  scale_fill_manual(values = c("CTRL.FMT" = "#08306b", "abx+FMTCTR.D0" = "#9ecae1", "abx+FMTCTR.D32" = "#08519c", 
                               "DNBS.FMT" = "#67000d", "abx+FMTDNBS.D0" = "#fc9272", "abx+FMTDNBS.D32" = "#a50f15"), 
                    "Legend") +
  scale_alpha_manual(values = c("FMT Pool" = 0, "fecal sample" = 1)) + 
  scale_colour_manual(values = c("abx+FMTCTR" = "#3690c0", "abx+FMTDNBS" = "#ef6548")) +
  facet_wrap(~type) +
  #Adjust labels
  xlab(paste("PC1: ", pc1,  "%", sep="")) + 
  ylab(paste("PC2: ", pc2,  "%", sep="")) + 
  theme_bw() + guides(colour = "none", alpha = "none") + 
  theme(text = element_text(size = 14))



# Differential abundance analysis -----------------------------------------
genus.glm =   fw_glmer(x             = genus_perm,
                       f             = ~ Treatment * Timepoint + (1|animal_ID), 
                       metadata      = meta_perm, 
                       adjust.method = "BH", order = "ac")

write.csv(genus.glm, file = "stats_pathologic_arm/genus_lv_lmer_stats.csv")

#hist(genus.glm$`anovas.Treatment:Timepoint Pr(>F).BH`)


genBH <-  genus.exp[genus.glm[genus.glm$`anovas.Treatment:Timepoint Pr(>F).BH` < 0.2,"feature"],]

#Plot the features that show a group effect at q < 0.2

p_genus <- genBH %>%
  t() %>%
  as.data.frame() %>%
  add_column(Timepoint = metadata$Timepoint, 
             Treatment = metadata$Treatment)  %>%
  pivot_longer(!c("Timepoint", "Treatment"))  %>%
  
  left_join(., genus.glm %>% dplyr::select("feature", "anovas.Treatment:Timepoint Pr(>F).BH",
                                           "anovas.Treatment:Timepoint Pr(>F)"), 
            by = c("name" = "feature")) %>% 
  
  mutate(stars = case_when(`anovas.Treatment:Timepoint Pr(>F)` < 0.001  ~ "⁂",
                           `anovas.Treatment:Timepoint Pr(>F)` < 0.01  ~ "🞱🞱",
                           `anovas.Treatment:Timepoint Pr(>F)` < 0.05 ~ "🞱", .default = ""),  
           name = str_replace(name, ".*ales_", "")) %>%  

  ggplot(aes(x     = Timepoint, 
             y     = value, 
             fill  = interaction(Treatment, Timepoint))) + 
  geom_boxplot(alpha = 1/2, coef = 100, show.legend = FALSE) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.75), show.legend = FALSE) + 
  
  geom_text(data = . %>% 
              dplyr::select(name, `anovas.Treatment:Timepoint Pr(>F)`, stars) %>% 
              distinct(name, .keep_all = TRUE), 
            inherit.aes = FALSE, show.legend = FALSE,
            aes(label = stars,  alpha = `anovas.Treatment:Timepoint Pr(>F)` < 0.05),
            x = 2.5, y = Inf, vjust = 2, color = "black", size = 4) +
  geom_vline(xintercept = 1.5, colour = "darkgray", linetype = "dashed", linewidth = 1) +
  
  facet_wrap(~name, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("CTRL.FMT" = "#08306b", "abx+FMTCTR.D0"  = "#9ecae1", "abx+FMTCTR.D7"  = "#4292c6", "abx+FMTCTR.D32"  = "#08519c", 
                               "DNBS.FMT" = "#67000d", "abx+FMTDNBS.D0" = "#fc9272", "abx+FMTDNBS.D7" = "#ef3b2c", "abx+FMTDNBS.D32" = "#a50f15"), "Legend") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) +
    guides(shape = "none", fill = guide_legend(override.aes = list(shape = c(23, 23, 21, 21, 21, 21)))) +
  
  ylab("Genus-level abundance (CLR)") + xlab(NULL) + theme_bw() + 
  theme(text = element_text(size = 14))

# GBMs --------------------------------------------------------------------

#Load GBMS like we did with the genus-level counts and metadata above. 
GBMs   <- read.delim("raw/GBMs_siobhain.csv", sep = ",", row.names = 3, header = T)

#Reorder the raw feature table to match the metadata
GBMs <- GBMs[,metadata$ID]

#Make sure our count data is all numbers
GBMs   <- apply(GBMs,c(1,2),function(x) as.numeric(as.character(x)))

#Remove features with prevalence < 10% in two steps:
#First, determine how often every feature is absent in a sample
n_zeroes_GBMs <- rowSums(GBMs == 0)

#Then, remove features that are absent in more than your threshold (90% in this case).
GBMs    <- GBMs[n_zeroes_GBMs <= round(ncol(GBMs) * 0.90),]  

#Perform a CLR transformation
GBMs.exp <- clr_c(GBMs)

GBMs_perm = GBMs.exp[,meta_perm$ID]


#This function fits the equivalent of lmer(feature ~ Treatment * Timepoint + (1|animal_ID)) for each feature.
#It also performs an appropriate Benjamini-Hochberg correction on the p-values. 
GBMs.glm =   fw_glmer(x             = GBMs_perm,
                      f             = ~ Treatment * Timepoint + (1|animal_ID), 
                      metadata      = meta_perm, 
                      adjust.method = "BH", order = "ac")

write.csv(GBMs.glm, file = "stats_pathologic_arm/GBMs_lv_lmer_stats.csv")

#hist(GBMs.glm$`anovas.Treatment:Timepoint Pr(>F).BH`)

#~  significant interactions in the GBM-level data

GBMBH <-  GBMs.exp[GBMs.glm[GBMs.glm$`anovas.Treatment:Timepoint Pr(>F).BH` < 0.1,"feature"],]

#Plot the features that show a group effect at q < 0.1 (none)

# GBMBH %>%
#   t() %>%
#   as.data.frame() %>%
#   add_column(Timepoint = metadata$Timepoint, 
#              Treatment = metadata$Treatment)  %>%
#   pivot_longer(!c("Timepoint", "Treatment"))  %>%
#   mutate(name = str_replace(name, ".*ales_", "")) %>% 
#   ggplot(aes(x     = Timepoint, 
#              y     = value, 
#              fill  = interaction(Treatment, Timepoint))) + 
#   geom_boxplot(alpha = 1/2, coef = 100) +
#   geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.75)) + 
#   
#   facet_wrap(~name, scales = "free_y", ncol = 4) +
#   scale_fill_manual(values = c("CTRL.FMT" = "#08306b", "abx+FMTCTR.D0"  = "#9ecae1", "abx+FMTCTR.D7"  = "#4292c6", "abx+FMTCTR.D32"  = "#08519c", 
#                                "DNBS.FMT" = "#67000d", "abx+FMTDNBS.D0" = "#fc9272", "abx+FMTDNBS.D7" = "#ef3b2c", "abx+FMTDNBS.D32" = "#a50f15"), "Legend") +
#   ylab(NULL) + xlab(NULL) + theme_bw() + theme(text = element_text(size = 14))

# GMMs --------------------------------------------------------------------

#Load GMMS like we did with the genus-level counts and metadata above. 
GMMs   <- read.delim("raw/GMMs_siobhain.csv", sep = ",", row.names = 3, header = T)

#Reorder the raw feature table to match the metadata
GMMs <- GMMs[,metadata$ID]

#Make sure our count data is all numbers
GMMs   <- apply(GMMs,c(1,2),function(x) as.numeric(as.character(x)))

#Remove features with prevalence < 10% in two steps:
#First, determine how often every feature is absent in a sample
n_zeroes_GMMs <- rowSums(GMMs == 0)

#Then, remove features that are absent in more than your threshold (90% in this case).
GMMs    <- GMMs[n_zeroes_GMMs <= round(ncol(GMMs) * 0.90),]  

#Perform a CLR transformation
GMMs.exp <- clr_c(GMMs)

GMMs_perm = GMMs.exp[,meta_perm$ID]


#This function fits the equivalent of lmer(feature ~ TypeCS * Time + (1|Participant)) for each feature.
#It also performs an appropriate Benjamini-Hochberg correction on the p-values. 
GMMs.glm =   fw_glmer(x             = GMMs_perm,
                      f             = ~ Treatment * Timepoint + (1|animal_ID), 
                      metadata      = meta_perm, 
                      adjust.method = "BH", order = "ac")
write.csv(GMMs.glm, file = "stats_pathologic_arm/GMMs_lv_lmer_stats.csv")

#hist(GMMs.glm$`anovas.Treatment:Timepoint Pr(>F).BH`)

#~ No interactions detected 

GMMBH <-  GMMs.exp[GMMs.glm[GMMs.glm$`anovas.Treatment:Timepoint Pr(>F).BH` < 0.1,"feature"],]

#Plot the features that show a group effect at q < 0.1

# GMMBH %>%
#   t() %>%
#   as.data.frame() %>%
#   add_column(Timepoint = metadata$Timepoint, 
#              Treatment = metadata$Treatment)  %>%
#   pivot_longer(!c("Timepoint", "Treatment"))  %>%
#   mutate(name = str_replace(name, ".*ales_", "")) %>% 
#   ggplot(aes(x     = Timepoint, 
#              y     = value, 
#              fill  = interaction(Treatment, Timepoint))) + 
#   geom_boxplot(alpha = 1/2, coef = 100) +
#   geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.75)) + 
#   
#   facet_wrap(~name, scales = "free_y", ncol = 4) +
#   scale_fill_manual(values = c("CTRL.FMT" = "#08306b", "abx+FMTCTR.D0"  = "#9ecae1", "abx+FMTCTR.D7"  = "#4292c6", "abx+FMTCTR.D32"  = "#08519c", 
#                                "DNBS.FMT" = "#67000d", "abx+FMTDNBS.D0" = "#fc9272", "abx+FMTDNBS.D7" = "#ef3b2c", "abx+FMTDNBS.D32" = "#a50f15"), "Legend") +
#   ylab(NULL) + xlab(NULL) + theme_bw() + theme(text = element_text(size = 14))




# anansi ------------------
# library(anansi)
# 
# metab    <- read.delim("raw/Pathological_FMT_for_anansi.csv", sep = ",", row.names = 1, header = T)
# KOs      <- read.delim("raw/pred_metagenome_unstrat.tsv", sep = "\t", row.names = 1, header = T)
# metadata <- read.delim("raw/metadata_DNBS.csv", sep = ",")
# 
# 
# #Fix metadata in accordance with genus names
# metadata$ID <- gsub(metadata$ID, pattern = "-", replacement = ".")
# metadata = metadata[!metadata$ID %in% c("DNBS.17", "DNBS.29", "DNBS.38", "DNBS.41", "DNBS.8"),]
# 
# metadata$Timepoint = factor(metadata$Timepoint, levels = c("FMT", "D0", "D7", "D32"))
# metadata = metadata[metadata$Experiment == "Pathological",]
# 
# metadata = metadata[metadata$Timepoint == "D32",]
# 
# metadata = metadata[!metadata$ID %in% c("DNBS.4", "DNBS.15", "DNBS.28", "DNBS.33", "DNBS.43", "DNBS.100", "DNBS.104"),]
# 
# KOs      = KOs[,metadata$ID]
# 
# KOs.exp  = clr_c(KOs)
# 
# t2       = t(KOs.exp)
# 
# metab    = metab[,metadata$ID]
# 
# metab.exp = clr_c(metab) 
# 
# t1       = t(metab.exp)
# 
# data(dictionary)
# anansi_dic = anansi_dic
# web <- weaveWebFromTables(tableY = t1, tableX = t2, dictionary = anansi_dic)
# 
# metadata$Treatment_fix <- str_replace(metadata$Treatment,pattern =  " \\+ ", replacement = "_")
# 
# anansi_out <- anansi(web      = web,             #Generated above
#                      formula  = ~ Treatment_fix, #Define formula to be fitted 
#                      metadata = metadata,        #where is the metadata
#                      verbose  = T                #To let you know what's happening
# )
# 
# anansiLong <- spinToLong(anansi_output = anansi_out, translate = T, 
#                          Y_translation = anansi::cpd_translation, 
#                          X_translation = anansi::KO_translation)  
# #Now it's ready to be plugged into ggplot2, though let's clean up a bit more. 
# 
# #Only consider interactions where the entire model fits well enough. 
# anansiLong <- anansiLong[anansiLong$model_full_q.values < 0.1,]
# 
# p_anansi_corplot <- ggplot(data = anansiLong, 
#                          aes(x      = r.values, 
#                              y      = feature_X, 
#                              fill   = type, 
#                              alpha  = model_disjointed_Treatment_p.values < 0.05)) + 
#   
#   #Make a vertical dashed red line at x = 0
#   geom_vline(xintercept = 0, linetype = "dashed", colour = "red")+
#   
#   #Points show  raw correlation coefficients
#   geom_point(shape = 21, size = 3) + 
#   
#   #facet per compound
#   ggforce::facet_col(~feature_Y, space = "free", scales = "free_y") + 
#   
#   #fix the scales, labels, theme and other layout
#   scale_y_discrete(limits = rev, position = "right") +
#   scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 1/3)) +
#   scale_fill_manual(values = c("abx+FMTCTR" = "#2166ac", 
#                                "abx+FMTDNBS"  = "#b2182b", 
#                                "All"        = "gray"))+
#   theme_bw() + 
#   ylab(NULL) + 
#   xlab("Pearson's rho")
# 
# 
# 
# outPlots = spinToPlots(anansi_out,
#                        target = anansi_out@input@web@dictionary &
#                          anansi_out@output@model_results$modelfit@q.values < 0.1, 
#                        Y_translation = anansi::cpd_translation, 
#                        X_translation = anansi::KO_translation, translate = T )
# 
# #load ggplot2 and patchwork for plotting
# 
# library(ggplot2)
# library(patchwork)
# 
# plotted = lapply(outPlots, FUN = function(p){
#   
#   #Main ggplot call
#   ggplot(data = p$data, aes(x = X, y = Y, fill = groups)) +
#     
#     #Establish geoms:
#     geom_point(shape = 21) +
#     geom_smooth(method = "lm") +
#     theme_bw() +
#     
#     #Improve annotation:
#     scale_fill_manual(values = c("abx+FMTCTR"  = "#2166ac", 
#                                  "abx+FMTDNBS" = "#b2182b"))+
#     ylab(p$name[1]) +
#     xlab(p$name[2]) +
#     ggtitle(paste(p$name[1], "vs", p$name[2]))
#   
# })
# 
# #Call patchwork to unify and arrange the first plot
# 
# p_anansi_dotplot <- wrap_plots(plotted) + plot_layout(guides = 'collect')
# 
