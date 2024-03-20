#.rs.restartR()
# script for generating a PCA with profiles generated using enviMass
rm(list=ls(all=TRUE)) 
setwd("C:/Users/slr257/Documents/Project 3/All_objectives_P3WWMetabo")

library(enviMass)
library(psych)
library(devtools)
library(FactoMineR)
library(factoextra)
library(reshape2)
library(dplyr)
library(pheatmap)
library(tidyr)
library(ggpubr)
# install_github("vqv/ggbiplot")

#### LOAD DATA ####

# project_name <- "P3_WW_Metabo_230827"
# project_name <- "P3_WW_Metabo_230907" # this project has target data
# project_name <- "P3_WW_Metabo_230917"
# project_name <- "P3_WW_Metabo_230918_Copy"
# project_name <- "P3_WW_Metabo_231001_istdnormchange"

project_name <- "P3_WW_Metabo_231008"

# project_name <- "P3_WW_Metabo_231008_Copy_ArchivedonOct15"

# original results, componentized with "reduce_comp = "greedy"" and enviMass project "P3_WW_Metabo_230918"
results_folder <- "results_new_workflow_231013"

# # try with noncomponentized peaks "reduce_comp = FALSE" and enviMass project "P3_WW_Metabo_230918"
# results_folder <- "results_noncomponentized_232001"

# # # try with componentized peaks "reduce_comp = "greedy"" and enviMass project "P3_WW_Metabo_231001" where ISTD normalization threshold
# # # is changed to 8 samples....
# results_folder <- "results_istd_norm_231001"


load_profileList(file = file.path("C:/Users/slr257/Documents/Project 3/enviMass/",project_name,"results/profileList_pos"))
load(file = file.path("C:/Users/slr257/Documents/Project 3/enviMass/",project_name,"results/links_profiles_pos"))
measurements<-read.csv(file=file.path("C:/Users/slr257/Documents/Project 3/enviMass/",project_name,"/dataframes/measurements"), colClasses = "character")
# write.csv(measurements, "real_measurements.csv")

# this matrix contains all features as columns, samples as rows, and intensity as entries
# function description for "profiles_to_matrix" can be found here: https://www.envibee.ch/eng/enviMass/faqs/faqs_3_3.php
# "reduce_comp" should either be "greedy" or "balanced" for componentized profiles
matrix <- as.data.frame(profiles_to_matrix(profileList_pos, links_profiles_pos, sort_by = c("number_peaks_sample", "mean_int"), sort_decreasing = TRUE, 
                             reduce_comp = "greedy", n_profiles = NULL, only_sample_peaks = F, n_latest_peaks = NULL, mean_above_blind = NULL, normalize = F))

# extracted peak information from profileList_pos
base_data <- as.data.frame(profileList_pos$peaks)
base_data <- base_data[base_data$profileIDs %in% colnames(matrix),] #includes only componentized profIDs
n_all_profiles <- length(unique(base_data$profileIDs))

# sampleID, date, and location information
sample_info <- measurements[,c(1,2,5,6)]
names(sample_info)[names(sample_info) == 'ID'] <- 'sampleIDs'

# associates peakIDs, and mz values with sample information (date, time, location)
base_data <- merge(base_data, sample_info, by = 'sampleIDs')
colnames(base_data)[2] <- c("mz")

#remove triplicates from base_data by taking average values for m/z, intensity, and RT for each sample date and location
base_data_notrip <- base_data %>%
  group_by(profileIDs, Place, Date) %>%
  summarise(mean_mz = mean(mz), mean_RT = mean(RT), mean_int = mean(intensity))

# extracts profile information (profile IDs) from profileList_pos
profile_summary <- profileList_pos$index_pro

colnames(profile_summary)[4] <- "profileIDs"
# information for each sample and profile ID along with date, mean_mz_location...
working_data <- as.data.frame(merge(base_data, profile_summary, by = "profileIDs")[,c(1,2,3,4,5,11,12,15,25,26,27)])

# only keep distinct IDs in a new working data dataframe
wd_distinctID <- working_data %>% distinct(profileIDs, .keep_all = TRUE)

# subset data based on location and only include unique profileIDs
INF_profiles <- working_data[working_data$Place == "INF",] %>% distinct(profileIDs, .keep_all = TRUE)
EFF_profiles <- working_data[working_data$Place == "EFF",] %>% distinct(profileIDs, .keep_all = TRUE)


##### GROUP BY NUMBER OF DETECTIONS AND LOCATION ####
matrix.df.count <- matrix
matrix.df.count$sampleIDs <- rownames(matrix.df.count)
matrix.df.count <- as.data.frame(merge(matrix.df.count, sample_info, by = "sampleIDs"))

# NEED TO CHANGE THIS WHEN ENVIMASS IS DONE PROCESSING NEW DATA
# take averages of triplicates for INF and EFF locations
ncol_1 <- nrow(wd_distinctID)+1


mat.df.INF <- pivot_longer(matrix.df.count[matrix.df.count$Place == "INF",], cols = 2:ncol_1, 
                                  names_to = "profile_ID",values_to = "intensity")

mat.df.INF <- mat.df.INF %>%
  group_by(Date, profile_ID, Place) %>%
  summarise(mean_int = mean(intensity))

mat.df.EFF <- pivot_longer(matrix.df.count[matrix.df.count$Place == "EFF",], cols = 2:ncol_1,
                                  names_to = "profile_ID",values_to = "intensity")

mat.df.EFF <- mat.df.EFF %>%
  group_by(Date, profile_ID, Place) %>%
  summarise(mean_int = mean(intensity))


all.intensities <- rbind(mat.df.INF, mat.df.EFF)

library(stats)
daily.intensities <- all.intensities %>%
  group_by(Date, Place) %>%
  summarise(mean_int2 = mean(mean_int), sd_int = sd(mean_int))

daily.intensities$Place <- factor(daily.intensities$Place, levels = c("INF","EFF"))
ggplot(data=daily.intensities, aes(x=Date, y = log10(mean_int2), fill = Place, color = Place)) +
  geom_point()

# relationship between intensities in INF and EFF? - pearson shows if there is a linear relationship
cor.test(daily.intensities$mean_int2[daily.intensities$Place == "INF"], daily.intensities$mean_int2[daily.intensities$Place == "EFF"],
         # alternative = c("two.sided", "less", "greater"),
         method = c("pearson"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)





# define number of detections
mat.EFF.detects <- mat.df.EFF %>%
  group_by(profile_ID) %>%
  summarise(EFF = sum(mean_int != 0))

# define number of detections
mat.INF.detects <- mat.df.INF %>%
  group_by(profile_ID) %>%
  summarise(INF = sum(mean_int != 0))

all.detects <- melt(merge(mat.INF.detects, mat.EFF.detects, by = "profile_ID"))
colnames(all.detects) <- c("profile_ID","Location","Count")
ggplot(data=all.detects, aes(x=Count, fill = Location)) +
  geom_bar(position = "dodge")
  

# define types in new "all.detects" dataframe
all.detects.type <- merge(mat.INF.detects, mat.EFF.detects, by = "profile_ID")

# careful with the order of ifelse statements
all.detects.type$type <- with(all.detects.type, ifelse(INF == 0 & EFF != 0, 'EFF Only',
                              ifelse(EFF == 0 & INF != 0, 'INF Only',
                                     # ifelse(INF == 14 & EFF != 14, "INF14",
                                     #        ifelse(INF != 14 & EFF == 14, "EFF14",
                                              ifelse(INF == 14 & EFF == 14, 'Core', "Intermittent"))))

type.date <- merge(rbind(mat.df.INF,mat.df.EFF),all.detects.type, by = "profile_ID")
colnames(base_data_notrip)[1] <- "profile_ID"

base_data_notrip_type <- merge(base_data_notrip, type.date[,c(1,7)], by = "profile_ID")


# add intensity classification based on individual meeting on Oct 4th 2023

summ_peakarea <- base_data_notrip_type %>%
  group_by(profile_ID, Place,Date, type)%>%
  summarise(mean_peakarea = mean(mean_int))


# plan - separate summ_peakarea into INF and EFF dataframes
# use for-loop to compare peak areas of same profile IDs in INF and EFF - IF they exist in both INF and EFF
# if profile IDs do not exist in INF and EFF, new category == old cateory (probably INF Only or EFF Only)

# proIDs <- unique(summ_peakarea$profile_ID)
# summ_peakarea$type_new <- rep(NA,nrow(summ_peakarea))

# for(ID in proIDs){
#   # ID = 2
#   INF_dat <- summ_peakarea[summ_peakarea$Place == "INF" & summ_peakarea$profile_ID == ID,]
#   EFF_dat <- summ_peakarea[summ_peakarea$Place == "EFF" & summ_peakarea$profile_ID == ID,]
# 
#   if(nrow(INF_dat) > 0 & nrow(EFF_dat) > 0){
#     if(mean(INF_dat$mean_peakarea) < 0.9*mean((EFF_dat$mean_peakarea))){
#       summ_peakarea$type_new[summ_peakarea$profile_ID == ID] <- "Form"
# 
#     }else if(mean(INF_dat$mean_peakarea) > 1.1*mean((EFF_dat$mean_peakarea))){
#       summ_peakarea$type_new[summ_peakarea$profile_ID == ID] <- "Biotrans"
# 
#     }else if(mean(INF_dat$mean_peakarea) < 1.1*mean((EFF_dat$mean_peakarea)) & mean(INF_dat$mean_peakarea) > 0.9*mean((EFF_dat$mean_peakarea))){
#       summ_peakarea$type_new[summ_peakarea$profile_ID == ID] <- "No Change"
#     }
#   }else{
#     summ_peakarea$type_new[summ_peakarea$profile_ID == ID] <- summ_peakarea$type[summ_peakarea$profile_ID == ID]
#   }
#   print(ID)
# }


#save(summ_peakarea, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/summ_peakarea.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/summ_peakarea.RData"))


ggplot(data=summ_peakarea, aes(x=Date, fill = type_new)) +
  geom_histogram(stat="count")+
  facet_wrap(~Place)

summ_peakarea$Place <- factor(summ_peakarea$Place, levels = c("INF","EFF"))
summ_peakarea$type_new <- factor(summ_peakarea$type_new,levels = c("INF Only","EFF Only","Form","Biotrans","No Change"))
summ_peakarea$Date <- with(summ_peakarea, ifelse(Date == "2020-11-18","Nov18",
                                         ifelse(Date == "2020-11-19","Nov19",
                                                ifelse(Date == "2020-11-20","Nov20",
                                                       ifelse(Date == "2020-11-21","Nov21",
                                                              ifelse(Date == "2020-11-22","Nov22",
                                                                     ifelse(Date == "2020-11-23","Nov23",
                                                                            ifelse(Date == "2020-11-24","Nov24",
                                                                                   ifelse(Date == "2020-11-25","Nov25",
                                                                                          ifelse(Date == "2020-11-26","Nov26",
                                                                                                 ifelse(Date == "2020-11-27","Nov27",
                                                                                                        ifelse(Date == "2020-11-28","Nov28",
                                                                                                               ifelse(Date == "2020-11-29","Nov29",
                                                                                                                      ifelse(Date == "2020-11-30","Nov30",
                                                                                                                             ifelse(Date == "2020-12-01","Dec1","x")))))))))))))))
summ_peakarea$Date <- factor(summ_peakarea$Date, levels = c("Nov18","Nov19","Nov20","Nov21","Nov22",
                                                    "Nov23","Nov24","Nov25","Nov26","Nov27",
                                                    "Nov28","Nov29","Nov30","Dec1"))

summ_peakarea$Place <- factor(summ_peakarea$Place, levels = c("INF","EFF"))
summ_peakarea$type_new <- factor(summ_peakarea$type_new,levels = c("INF Only","EFF Only","Form","Biotrans","No Change"))

# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/featurecountsbytypebyday_FIG1_231008.tiff"), units="in", width=14, height=8, res=300)
ggplot(data=summ_peakarea, aes(x=Date, fill = type_new)) +
  geom_histogram(stat="count")+
  facet_wrap(~Place)+
  ylab("Count")+
  scale_fill_manual(values = c("#FF6666","#00BFC4","darkgreen","purple","orange"))+
  theme(text = element_text(size=18),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.size = unit(0.32, 'cm'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray97"),
        panel.grid.minor = element_line(color = "gray97"),
        panel.border = element_rect(color = "grey", fill = NA))
dev.off()

#### FIG1 data exported in a table format ####
type_counts <- summ_peakarea %>%
  group_by(Date, Place, type_new) %>%
  summarise(n = n())

type_counts_w <- spread(type_counts, type_new,n)
write.csv(type_counts_w,"Fig1_typecounts_table.csv")

# try relative abundance plot for features intensities using based on new categories

# not sure if this is the right approach
relabund_int <- summ_peakarea %>%
  group_by(Date, Place) %>%
  summarise(total = sum(mean_peakarea))

relabund_int2 <- merge(summ_peakarea,relabund_int, by = c("Date","Place")) 
relabund_int2$relabund <- (relabund_int2$mean_peakarea/relabund_int2$total)*100

#subset date, place, ID plus new type and relabund
x5 <- relabund_int2[,c(1:3,6,8)]

relabund_int2$profile_ID <- as.factor(relabund_int2$profile_ID)

# relabund.m <- relabund_int2 %>% pivot_longer(cols = c("Date","Place","profile_ID","relabund"))


ggplot(data=x5, aes(x=Date, fill = profile_ID))+
  geom_bar()

x5$profile_ID <- as.factor(x5$profile_ID)

x5$type_new <- factor(x5$type_new,levels = c("INF Only","EFF Only","Form","Biotrans","No Change"))
x5$Place <- factor(x5$Place, levels = c("INF","EFF"))
colnames(x5) <- c("Date","Place","profile_ID","Type","Relative_Abundance")

x5$Date <- factor(x5$Date, levels = c("Nov18","Nov19","Nov20","Nov21","Nov22",
                                      "Nov23","Nov24","Nov25","Nov26","Nov27",
                                      "Nov28","Nov29","Nov30","Dec1"))

#### PLOT WITH RELATIVE ABUNDANCE BY NEW CATEGORIES - MAYBE ADD CORE VS NON CORE COMPARISON?
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/relative_abundance_FIG2_231008.tiff"), units="in", width=14, height=8, res=300)
# ggplot(x5, aes(x = Date, y = Relative_Abundance, fill = Type))+
#   geom_bar(stat = "identity")+
#   facet_wrap(~Place)+
#   ylab("Relative Abundance (%)")+
#   scale_fill_manual(values = c("#FF6666","#00BFC4","darkgreen","purple","orange"))+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.32, 'cm'),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()



#subset date, place, ID plus new type and relabund
x6 <- relabund_int2[,c(1:3,4,8)]


x6$profile_ID <- as.factor(x6$profile_ID)

x6$type <- factor(x6$type,levels = c("INF Only","EFF Only","Intermittent","Core"))
x6$Place <- factor(x6$Place, levels = c("INF","EFF"))
colnames(x6) <- c("Date","Place","profile_ID","Type","Relative_Abundance")

x6$Date <- factor(x6$Date, levels = c("Nov18","Nov19","Nov20","Nov21","Nov22",
                                      "Nov23","Nov24","Nov25","Nov26","Nov27",
                                      "Nov28","Nov29","Nov30","Dec1"))




# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/relative_abundance_corenoncore_231008.tiff"), units="in", width=14, height=8, res=300)
# ggplot(x6, aes(x = Date, y = Relative_Abundance, fill = Type))+
#   geom_bar(stat = "identity")+
#   facet_wrap(~Place)+
#   ylab("Relative Abundance (%)")+
#   scale_fill_manual(values = c("#FF6666","#00BFC4","purple","orange"))+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.32, 'cm'),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()

# TRY COMPARING RELABUND OF CORE VS NON CORE

#subset date, place, ID plus new type and relabund
x7 <- relabund_int2[,c(1:3,4,6,8)]

library(ggpattern) # maybe the pattern thing will work one day?
library(ggplot2)
#[compare_bars$Category == "type_new",]

x7$type[x7$type == "INF Only"] <- "INF_Only"
x7$type[x7$type == "EFF Only"] <- "EFF_Only"

x7$type <- factor(x7$type, levels = c("INF_Only","EFF_Only","Intermittent","Core"))

#### THIS IS THE PLOT TO KEEP! ####
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/relative_abundance_all_info_231015.tiff"), units="in", width=14, height=8, res=300)
ggplot(x7, aes(x = Date, y = relabund, fill = type_new, alpha = type))+
  geom_bar(stat = "identity")+
  facet_wrap(~Place)+
  ylab("Relative Abundance (%)")+
  scale_alpha_manual(values=c(1, 1, 0.6, 1), guide = 'none') +
  scale_fill_manual(values = c("#FF6666","#00BFC4","darkgreen","purple","orange"))+
  theme(text = element_text(size=18),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.size = unit(0.32, 'cm'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray97"),
        panel.grid.minor = element_line(color = "gray97"),
        panel.border = element_rect(color = "grey", fill = NA))
# dev.off()
####


#### FIG2 data exported in a table format ####
type_relabundFIG2 <- x7 %>%
  group_by(Date, Place, type_new) %>%
  summarise(sum_relabund = sum(relabund))

type_relabundFIG2_w <- spread(type_relabundFIG2, type_new,sum_relabund)
write.csv(type_relabundFIG2_w,"Fig2_relabund_table.csv")

####

n_core <- length(unique(x7$profile_ID[x7$type == "Core"]))
n_INFOnly <- length(unique(x7$profile_ID[x7$type == "INF_Only"]))
n_EFFOnly <- length(unique(x7$profile_ID[x7$type == "EFF_Only"]))
n_Inter <- length(unique(x7$profile_ID[x7$type == "Intermittent"]))
n_Form <- length(unique(x7$profile_ID[x7$type_new == "Form"]))
n_Biotrans <- length(unique(x7$profile_ID[x7$type_new == "Biotrans"]))
n_nochange <- length(unique(x7$profile_ID[x7$type_new == "No Change"]))

# ggplot(x7, aes(x = Date, y = relabund, fill = type_new, pattern = type))+
#   geom_bar_pattern(inherit.aes = TRUE,stat = "identity",position = "stack",
#                    color = "black", 
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.2,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 0.6) +
#   facet_wrap(~Place)+
#   ylab("Relative Abundance (%)")+
#   scale_pattern_manual(values = c(INF_Only = "none", EFF_Only = "none",Intermittent = "none",Core = "stripe")) +
#   scale_fill_manual(values = c("#FF6666","#00BFC4","darkgreen","purple","orange"))+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.32, 'cm'),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))



# 5. run bar plot workflow
summ_type <- base_data_notrip_type %>%
  group_by(profile_ID, Date, Place, type)%>%
  summarise(n = n())
summ_type$Place <- factor(summ_type$Place, levels = c("INF","EFF"))
summ_type$type <- factor(summ_type$type,levels = c("INF Only","EFF Only","Intermittent","Core"))

summ_type$Date <- with(summ_type, ifelse(Date == "2020-11-18","Nov18",
                                         ifelse(Date == "2020-11-19","Nov19",
                                            ifelse(Date == "2020-11-20","Nov20",
                                                ifelse(Date == "2020-11-21","Nov21",
                                                      ifelse(Date == "2020-11-22","Nov22",
                                                            ifelse(Date == "2020-11-23","Nov23",
                                                                  ifelse(Date == "2020-11-24","Nov24",
                                                                          ifelse(Date == "2020-11-25","Nov25",
                                                                              ifelse(Date == "2020-11-26","Nov26",
                                                                                    ifelse(Date == "2020-11-27","Nov27",
                                                                                            ifelse(Date == "2020-11-28","Nov28",
                                                                                                  ifelse(Date == "2020-11-29","Nov29",
                                                                                                          ifelse(Date == "2020-11-30","Nov30",
                                                                                                                ifelse(Date == "2020-12-01","Dec1","x")))))))))))))))
summ_type$Date <- factor(summ_type$Date, levels = c("Nov18","Nov19","Nov20","Nov21","Nov22",
                                                    "Nov23","Nov24","Nov25","Nov26","Nov27",
                                                    "Nov28","Nov29","Nov30","Dec1"))

# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/featurecountsbytypebyday_FIG1_231015.tiff"), units="in", width=14, height=8, res=300)
# ggplot(data=summ_type, aes(x=Date, fill = type)) +
#   geom_histogram(stat="count")+
#   facet_wrap(~Place)+
#   ylab("Count")+
#   scale_fill_manual(values = c("#FF6666","#00BFC4","purple","orange"))+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.32, 'cm'),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()


#### test for correlation between total number of features in INF and EFF over time
overall_feat_counts_by_date <- summ_type %>%
  group_by(Place, Date) %>%
  summarize(n = n())

feat.counts.inf <- overall_feat_counts_by_date[overall_feat_counts_by_date$Place == "INF",]
feat.counts.eff <- overall_feat_counts_by_date[overall_feat_counts_by_date$Place == "EFF",]

# NO CORRELATION
cor.test(feat.counts.inf$n, feat.counts.eff$n,
         # alternative = c("two.sided", "less", "greater"),
         method = c("pearson"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)

sd.counts.inf <- sd(feat.counts.inf$n)
mean.count.inf <- mean(feat.counts.inf$n)
CV.counts.inf <- sd.counts.inf/mean.count.inf

sd.counts.eff <- sd(feat.counts.eff$n)
mean.count.eff <- mean(feat.counts.eff$n)
CV.counts.eff <- sd.counts.eff/mean.count.eff

# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/total_feature_counts_bylocation_231001.tiff"), 
#      units="in", width=12, height=8, res=300)
# ggplot(data=summ_type, aes(x=Date, fill = Place)) +
#   geom_histogram(stat="count", position = "dodge")+
#   # facet_wrap(~Place)+
#   ylab("Count")+
#   scale_fill_manual(values = c("#FF6666","#00BFC4"))+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.5, 'cm'),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()


p <- ggplot(data=summ_type, aes(x=Date, fill = Place)) +
  geom_histogram(stat="count", position = "dodge")
pg <- ggplot_build(p)
x<- pg$data[[1]]

# #values for writing
# n_core <- all.detects.type[all.detects.type$type == "Core",]
# n_INFOnly <- all.detects.type[all.detects.type$type == "INF Only",]
# n_EFFOnly <- all.detects.type[all.detects.type$type == "EFF Only",]
# n_Inter <- all.detects.type[all.detects.type$type == "Intermittent",]
  

# look into componentization - you probably need this for the final figure - greedy is okay



#### GENERATE HISTOGRAMS FOR MZ INTENSITY AND RT AT INF AND EFF LOCATIONS SEPARATELY #### ---------------
# get average mz, RT, and intensity for each profile_ID by location
sum_base_data <- base_data_notrip %>%
  group_by(profile_ID,Place) %>%
  summarise(avg_mz = mean(mean_mz),avg_RT = mean(mean_RT),avg_int = mean(mean_int))

sum_base_data$Place <- factor(sum_base_data$Place, levels = c("INF","EFF"))

# use "sum_base_data" to generate histograms of mz, RT, and int for INF and EFF
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/mz_histogram_231016.tiff"), units="in", width=10, height=8, res=300)
# ggplot(sum_base_data, aes(x=avg_mz, fill = Place))+ 
#   geom_histogram(position = "identity",alpha = 0.3)+
#   # facet_wrap(~Place)+
#   xlab("Average m/z by profile ID")+
#   ylab("Count")+
#   ylim(0,1000)+
#   labs(fill = "Location")+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         legend.key.size = unit(0.32, 'cm'),
#         axis.text.x = element_text(vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_text(size = 16),
#         # axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()
# 
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/RT_histogram_231016.tiff"), units="in", width=10, height=8, res=300)
# ggplot(sum_base_data, aes(x=avg_RT/60, fill = Place))+ 
#   geom_histogram(position = "identity",alpha = 0.3)+
#   # facet_wrap(~Place)+
#   xlab("Average RT (min) by profile ID")+
#   ylab("Count")+
#   ylim(0,1000)+
#   labs(fill = "Location")+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         legend.key.size = unit(0.32, 'cm'),
#         axis.text.x = element_text(vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_text(size = 16),
#         # axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()
# 
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/PeakArea_histogram_231016.tiff"), units="in", width=10, height=8, res=300)
# ggplot(sum_base_data, aes(x=log10(avg_int), fill = Place))+ 
#   geom_histogram(position = "identity",alpha = 0.3)+
#   # facet_wrap(~Place)+
#   xlab("Average log Peak Area by chemical feature ID")+
#   ylab("Count")+
#   ylim(0,2000)+
#   labs(fill = "Location")+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         legend.key.size = unit(0.32, 'cm'),
#         axis.text.x = element_text(vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_text(size = 16),
#         # axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()

# T-TESTS FOR MZ RT AND PEAK AREAS IN INF AND EFF
sum_base_data_inf <- sum_base_data[sum_base_data$Place == "INF",]
sum_base_data_eff <- sum_base_data[sum_base_data$Place == "EFF",]

t.test(sum_base_data_inf$avg_mz, y = sum_base_data_eff$avg_mz,
       alternative = c("two.sided"),
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95)

t.test(sum_base_data_inf$avg_RT, y = sum_base_data_eff$avg_RT,
       alternative = c("two.sided"),
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95)

t.test(sum_base_data_inf$avg_int, y = sum_base_data_eff$avg_int,
       alternative = c("two.sided"),
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95)



# just to double check this is correct
length(unique(sum_base_data$profile_ID))



ggplot(base_data_notrip[base_data_notrip$Place == "INF",], aes(x=log10(mean_int), fill = Date))+ 
  geom_histogram(position = "identity",alpha = 0.3)+
  # facet_wrap(~Place)+
  xlab("Average log Peak Area by chemical feature ID")+
  ylab("Count")+
  ylim(0,700)+
  labs(fill = "Location")+
  theme(text = element_text(size=18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.32, 'cm'),
        axis.text.x = element_text(vjust = 1, hjust=1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        # axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray97"),
        panel.grid.minor = element_line(color = "gray97"),
        panel.border = element_rect(color = "grey", fill = NA))

ggplot(base_data_notrip[base_data_notrip$Place == "EFF",], aes(x=log10(mean_int), fill = Date))+ 
  geom_histogram(position = "identity",alpha = 0.3)+
  # facet_wrap(~Place)+
  xlab("Average log Peak Area by chemical feature ID")+
  ylab("Count")+
  ylim(0,700)+
  labs(fill = "Location")+
  theme(text = element_text(size=18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.32, 'cm'),
        axis.text.x = element_text(vjust = 1, hjust=1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        # axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray97"),
        panel.grid.minor = element_line(color = "gray97"),
        panel.border = element_rect(color = "grey", fill = NA))

ggplot(base_data_notrip[base_data_notrip$Place == "EFF",], aes(x = Date, y=log10(mean_int), fill = Date))+ 
  geom_boxplot(position = "identity")+
  # facet_wrap(~Place)+
  xlab(element_blank())+
  ylab("Count")+
  # ylim(0,700)+
  labs(fill = "Location")+
  theme(text = element_text(size=18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.32, 'cm'),
        axis.text.x = element_text(vjust = 1, hjust=1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        # axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray97"),
        panel.grid.minor = element_line(color = "gray97"),
        panel.border = element_rect(color = "grey", fill = NA))




#### TEST PEAK AREA CORRELATION ####

# group data by place and date, take average of peak area
new.peakarea <- base_data_notrip %>%
  group_by(Place,Date) %>%
  summarise(avg_PA = mean(mean_int))

bd_notrip_eff <- new.peakarea[new.peakarea$Place == "EFF",]
bd_notrip_inf <- new.peakarea[new.peakarea$Place == "INF",]

cor.test(bd_notrip_eff$avg_PA, bd_notrip_inf$avg_PA,
         # alternative = c("two.sided", "less", "greater"),
         method = c("pearson"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)


#### TEST CORRELATION FOR NUMBER OF CHEMICAL FEATURES PER DAY AND MEAN PEAK AREA
# cor.test(bd_notrip_eff$avg_PA, bd_notrip_inf$avg_PA,
#          # alternative = c("two.sided", "less", "greater"),
#          method = c("pearson"),
#          exact = NULL, conf.level = 0.95, continuity = FALSE)
# 
# cor.test(feat.counts.inf$n, feat.counts.eff$n,
#          # alternative = c("two.sided", "less", "greater"),
#          method = c("pearson"),
#          exact = NULL, conf.level = 0.95, continuity = FALSE)


cor.test(feat.counts.inf$n, bd_notrip_inf$avg_PA,
         # alternative = c("two.sided", "less", "greater"),
         method = c("spearman"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)

cor.test(feat.counts.eff$n, bd_notrip_eff$avg_PA,
         # alternative = c("two.sided", "less", "greater"),
         method = c("spearman"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)

plot(feat.counts.eff$n,bd_notrip_eff$avg_PA)
lm(feat.counts.eff$n ~ bd_notrip_eff$avg_PA)

plot(feat.counts.inf$n,bd_notrip_inf$avg_PA)
lm(feat.counts.inf$n ~ bd_notrip_inf$avg_PA)


new.peakarea$Place <- factor(new.peakarea$Place, levels = c("INF","EFF"))

new.peakarea$Date <- with(new.peakarea, ifelse(Date == "2020-11-18","Nov18",
                                         ifelse(Date == "2020-11-19","Nov19",
                                                ifelse(Date == "2020-11-20","Nov20",
                                                       ifelse(Date == "2020-11-21","Nov21",
                                                              ifelse(Date == "2020-11-22","Nov22",
                                                                     ifelse(Date == "2020-11-23","Nov23",
                                                                            ifelse(Date == "2020-11-24","Nov24",
                                                                                   ifelse(Date == "2020-11-25","Nov25",
                                                                                          ifelse(Date == "2020-11-26","Nov26",
                                                                                                 ifelse(Date == "2020-11-27","Nov27",
                                                                                                        ifelse(Date == "2020-11-28","Nov28",
                                                                                                               ifelse(Date == "2020-11-29","Nov29",
                                                                                                                      ifelse(Date == "2020-11-30","Nov30",
                                                                                                                             ifelse(Date == "2020-12-01","Dec1","x")))))))))))))))
new.peakarea$Date <- factor(new.peakarea$Date, levels = c("Nov18","Nov19","Nov20","Nov21","Nov22",
                                                    "Nov23","Nov24","Nov25","Nov26","Nov27",
                                                    "Nov28","Nov29","Nov30","Dec1"))
dates <- c("Nov18","Nov19","Nov20","Nov21","Nov22",
           "Nov23","Nov24","Nov25","Nov26","Nov27",
           "Nov28","Nov29","Nov30","Dec1")


# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/PeakArea_barchart_231016.tiff"), units="in", width=10, height=8, res=300)
# ggplot(data=new.peakarea, aes(x=Date, y = avg_PA, fill = Place)) +
#   geom_bar(stat = "identity", position = "dodge")+
#   # facet_wrap(~Place)+
#   ylab("Mean Peak Area")+
#   scale_fill_manual(values = c("#FF6666","#00BFC4"))+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.5, 'cm'),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()

ggbarplot(new.peakarea, x = "Date", y = "avg_PA",
          color = "Place", 
          position = position_dodge(0.8))+
  stat_compare_means(aes(group = Place), label = "p.signif", label.y = 29)

base_data_notrip1 <- base_data_notrip

base_data_notrip1$Date <- with(base_data_notrip1, ifelse(Date == "2020-11-18","Nov18",
                                               ifelse(Date == "2020-11-19","Nov19",
                                                      ifelse(Date == "2020-11-20","Nov20",
                                                             ifelse(Date == "2020-11-21","Nov21",
                                                                    ifelse(Date == "2020-11-22","Nov22",
                                                                           ifelse(Date == "2020-11-23","Nov23",
                                                                                  ifelse(Date == "2020-11-24","Nov24",
                                                                                         ifelse(Date == "2020-11-25","Nov25",
                                                                                                ifelse(Date == "2020-11-26","Nov26",
                                                                                                       ifelse(Date == "2020-11-27","Nov27",
                                                                                                              ifelse(Date == "2020-11-28","Nov28",
                                                                                                                     ifelse(Date == "2020-11-29","Nov29",
                                                                                                                            ifelse(Date == "2020-11-30","Nov30",
                                                                                                                                   ifelse(Date == "2020-12-01","Dec1","x")))))))))))))))
base_data_notrip1$Date <- factor(base_data_notrip1$Date, levels = c("Nov18","Nov19","Nov20","Nov21","Nov22",
                                                          "Nov23","Nov24","Nov25","Nov26","Nov27",
                                                          "Nov28","Nov29","Nov30","Dec1"))


## THIS IS ALSO A GOOD PLOT
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/PeakArea_PSIG_231016.tiff"), units="in", width=10, height=8, res=300)
# ggbarplot(base_data_notrip1, x = "Date", y = "mean_int",add = "mean_se",
#           add.params = list(color = "#4C4E52"),
#           color = "darkgrey", fill = "Place", position = position_dodge(0.8))+
#   stat_compare_means(aes(group = Place), label = "p.signif", label.y = 3.5*10^7,method = "t.test", hide.ns = T)+
#   ylab("Mean Peak Area")+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.5, 'cm'),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()

res <- desc_statby(base_data_notrip, measure.var = "mean_int",
                   grps = c("Date", "Place"))

res2 <- desc_statby(base_data_notrip, measure.var = "mean_int",
                   grps = c("Place"))

# get difference plot
diff <- bd_notrip_inf$avg_PA - bd_notrip_eff$avg_PA
diff.alldat <- as.data.frame(cbind(dates, diff))
diff.alldat$color <- ifelse(diff.alldat$diff < 0, "Negative","Positive")
diff.alldat$hjust <- ifelse(diff.alldat$diff > 0, 1.3, -0.3)

diff.alldat$dates <- factor(diff.alldat$dates, levels = c("Nov18","Nov19","Nov20","Nov21","Nov22",
                                                                    "Nov23","Nov24","Nov25","Nov26","Nov27",
                                                                    "Nov28","Nov29","Nov30","Dec1"))

# this is the plot to keep!!!
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/PeakArea_DIFF_231021.tiff"), units="in", width=10, height=8, res=300)
# ggplot(diff.alldat,aes(dates,as.numeric(diff),label="",hjust=hjust))+
#   geom_text(aes(y=0,color=color))+
#   geom_bar(stat="identity",position="identity",aes(fill = color))+
#   ylab("Difference in Mean Peak Area (INF - EFF)")+
#   xlab(element_blank())+
#   scale_fill_manual(values=c(Positive="red",Negative="blue"))+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.5, 'cm'),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()

# this analysis should actually just be boxplots showing distibutions of peak areas on every day of the campaign across INF and EFF

base_data_notrip$Place <- factor(base_data_notrip$Place, levels = c("INF","EFF"))


# PLOT TO KEEP
ggplot(data=base_data_notrip, aes(x=Date, y = log10(mean_int), fill = Place)) +
  geom_boxplot()+
  # geom_jitter()+
  # facet_wrap(~Place)+
  ylab("Mean Peak Area")+
  scale_fill_manual(values = c("#FF6666","#00BFC4"))+
  theme(text = element_text(size=18),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray97"),
        panel.grid.minor = element_line(color = "gray97"),
        panel.border = element_rect(color = "grey", fill = NA))+
  stat_compare_means(aes(group = Date), label = "p.signif", method = "t.test")


# compares mean peak area across INF and EFF on every date
# https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
compare_means(mean_int ~ Place, data = base_data_notrip, 
              group.by = "Date", method = "t.test")

ggboxplot(base_data_notrip, x = "Date", y = "mean_int", color = "Place")+
  stat_compare_means(aes(group = Place), label = "p.signif", method = "t.test", hide.ns = T)


ggplot(data=base_data_notrip, aes(x=Date, y = mean_int, fill = Place)) +
  geom_bar(stat = "identity",position = "dodge")+
  # facet_wrap(~Place)+
  ylab("Mean Peak Area")+
  scale_fill_manual(values = c("#FF6666","#00BFC4"))+
  theme(text = element_text(size=18),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray97"),
        panel.grid.minor = element_line(color = "gray97"),
        panel.border = element_rect(color = "grey", fill = NA))







#### HEATMAP ANALYSIS ####
core_proIDs <- as.data.frame(unique(x7$profile_ID[x7$type == "Core"]))
colnames(core_proIDs) <- c("profile_ID")
# turn matrix into dataframe
# matrix.df <- as.data.frame(matrix)
# matrix.df.comp <- as.data.frame(matrix)
matrix.df.greedy <- as.data.frame(matrix)

# transpose matrix to subset by column values (profileIDs are now factors)
greedy.t <- as.data.frame(t(matrix.df.greedy))
greedy.t$profileID <- rownames(greedy.t)
matrix.core <- greedy.t[greedy.t$profileID %in% core_proIDs$profile_ID,]

#check there are no all zero rows
checkzeros <- sum(rowSums(matrix.core[,c(1:84)]) == 0) # equals zero - good


matrix.core.t <- as.data.frame(t(matrix.core[,c(1:84)]))
matrix.core.t$sampleIDs <- rownames(matrix.core.t)
#only include certain columns from sample info dataframe
matrix.core.t <- merge(matrix.core.t, sample_info[,c(1,3,4)], by = "sampleIDs")

library(reshape2)

m <- melt(matrix.core.t)
colnames(m) <- c("ID","Place","Date","profileID","Intensity")

library(dplyr)
# one date and one location for each feature ID
core.grouped <- m %>%
  group_by(Place, Date, profileID) %>%
  summarise(mean_int = mean(Intensity))

core.grouped$Date <- with(core.grouped, ifelse(Date == "2020-11-18","Nov18",
                                         ifelse(Date == "2020-11-19","Nov19",
                                                ifelse(Date == "2020-11-20","Nov20",
                                                       ifelse(Date == "2020-11-21","Nov21",
                                                              ifelse(Date == "2020-11-22","Nov22",
                                                                     ifelse(Date == "2020-11-23","Nov23",
                                                                            ifelse(Date == "2020-11-24","Nov24",
                                                                                   ifelse(Date == "2020-11-25","Nov25",
                                                                                          ifelse(Date == "2020-11-26","Nov26",
                                                                                                 ifelse(Date == "2020-11-27","Nov27",
                                                                                                        ifelse(Date == "2020-11-28","Nov28",
                                                                                                               ifelse(Date == "2020-11-29","Nov29",
                                                                                                                      ifelse(Date == "2020-11-30","Nov30",
                                                                                                                             ifelse(Date == "2020-12-01","Dec1","x")))))))))))))))
core.grouped$Date <- factor(core.grouped$Date, levels = c("Nov18","Nov19","Nov20","Nov21","Nov22",
                                                    "Nov23","Nov24","Nov25","Nov26","Nov27",
                                                    "Nov28","Nov29","Nov30","Dec1"))
library(tidyr)
# now average intensities for triplicate measurements at each location and date
new.plot <- as.data.frame(spread(core.grouped, key = profileID, value = mean_int))
new.plot <- as.data.frame(new.plot %>% arrange(factor(Place, levels = c("INF","EFF"))))

# if(results_folder == "results_noncomponentized_232001"){
#   ncol_2 <- nrow(n_core)+2
# }else if(results_folder == "results_comp_nonorm_enviMass_230918"){
#   ncol_2 <- 360
# }

ncol_2 <- n_core+2
# pheatmap(t(new.plot[,c(3:ncol_2)]), scale = "row", labels_col = t(new.plot$Date))

annot_col2 <- data.frame(row.names = rownames(new.plot), Location = new.plot[,1])
rownames(new.plot) <- rownames(annot_col2)

# pheatmap(t(new.plot[,c(3:ncol_2)]), scale = "row",cluster_cols = F, labels_col = t(new.plot$Date), show_rownames = F,
#          annotation_col = annot_col2,clustering_method = "ward.D2")

# changing colors of heatmap with diverging scale
ncol <- 300

## Make a vector with n colors
# cols <- RColorBrewer:::brewer.pal(11,"PRGn")  # OR c("purple","white","orange")
cols <- c("red","white","#0012b0")
rampcols <- colorRampPalette(colors = rev(cols), space="Lab")(ncol)

# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/HCA_core_profiles_231015.tiff"), units="in", width=8, height=10, res=300)
# pheatmap(t(new.plot[,c(3:ncol_2)]), scale = "row",cluster_cols = F, labels_col = t(new.plot$Date), show_rownames = F,
#          annotation_col = annot_col2, cutree_rows = 3,clustering_method = "ward.D2", color = rampcols)
# dev.off()

# get optimal number of clusters here: https://www.statology.org/hierarchical-clustering-in-r/

library(cluster)
#calculate gap statistic for each number of clusters (up to 10 clusters)
gap_stat <- clusGap(t(new.plot[,c(3:ncol_2)]), FUN = hcut, nstart = 25, K.max = 10, B = 50)

# #produce plot of clusters vs. gap statistic
fviz_gap_stat(gap_stat)

# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/HCA_optimal_clust_231015.tiff"), units="in", width=10, height=8, res=300)
# fviz_nbclust(t(new.plot[,c(3:ncol_2)]),FUNcluster = hcut, method = "wss")
# dev.off()

# means of each profile
means <- rowMeans(matrix.core[,c(1:84)])


res <- pheatmap(t(new.plot[,c(3:ncol_2)]), scale = "row",cluster_cols = F, labels_col = t(new.plot$Date), show_rownames = F,
         annotation_col = annot_col2, cutree_rows = 3,clustering_method = "ward.D2", color = rampcols)

# get boxplots of intensity distributions for each cluster per day
core.clust <- as.data.frame(cbind(t(new.plot[,c(3:ncol_2)]), cluster = cutree(res$tree_row, k = 3)))

# featureIDs and cluster numbers
featID.clustnum <- as.data.frame(cbind(rownames(core.clust),core.clust$cluster))
colnames(featID.clustnum) <- c("profileID","cluster")

core.grouped.clust <- merge(core.grouped, featID.clustnum, by = "profileID")

core.grouped.clust$Place <- factor(core.grouped.clust$Place, levels = c("INF","EFF"))

# Figure out how clusters are labeled
ggplot(core.grouped.clust, aes(x=Date, y=log10(mean_int), fill = Place)) + 
  geom_boxplot()+
  facet_wrap(~cluster)

# cluster 1 here is same as in plot
# cluster 2 here is cluster 3 in plot
# cluster 3 here is cluster 2 in plot

core.grouped.clust$cluster <- factor(core.grouped.clust$cluster, levels = c("1", "2", "3"), 
                  labels = c("Cluster 1","Cluster 3","Cluster 2"))

core.grouped.clust$cluster <- factor(core.grouped.clust$cluster, levels = c("Cluster 1","Cluster 2","Cluster 3"))

# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/cluster_intensity_grid_231017.tiff"), units="in", width=12, height=8, res=300)
# ggplot(core.grouped.clust, aes(x=Date, y=log10(mean_int), fill = Place)) + 
#   geom_boxplot()+
#   facet_wrap(~cluster)+
#   ylab("Log (Peak Area)")+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         legend.key.size = unit(0.32, 'cm'),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()

num_feat_percluster <- core.grouped.clust %>%
  group_by(Place, Date, cluster) %>%
  summarise(n = n())

## find how many features in each cluster are associated with Form, Biotrans, and No Change groups

res <- pheatmap(t(new.plot[,c(3:ncol_2)]), scale = "row",cluster_cols = F, labels_col = t(new.plot$Date), show_rownames = F,
                annotation_col = annot_col2, cutree_rows = 3,clustering_method = "ward.D2", color = rampcols)
core.clust <- as.data.frame(cbind(t(new.plot[,c(3:ncol_2)]), cluster = cutree(res$tree_row, k = 3)))

# cluster numbers are reordered to match figure labels
prof_ID_clust1 <- as.numeric(rownames(core.clust[core.clust$cluster == 1,]))
prof_ID_clust3 <- as.numeric(rownames(core.clust[core.clust$cluster == 2,]))
prof_ID_clust2 <- as.numeric(rownames(core.clust[core.clust$cluster == 3,]))

proID_Biotrans <- as.numeric(unique(summ_peakarea$profile_ID[summ_peakarea$type_new == "Biotrans"]))
proID_Form <- as.numeric(unique(summ_peakarea$profile_ID[summ_peakarea$type_new == "Form"]))
proID_NoChange <- as.numeric(unique(summ_peakarea$profile_ID[summ_peakarea$type_new == "No Change"]))

# cluster 1 intersects
c1_bt <- prof_ID_clust1[prof_ID_clust1 %in% proID_Biotrans] # 109 out of 111 in Biotrans
c1_f <- prof_ID_clust1[prof_ID_clust1 %in% proID_Form] # 0 in Form
c1_nc <- prof_ID_clust1[prof_ID_clust1 %in% proID_NoChange] # 2 in No Change

# cluster 2 intersects
c2_bt <- prof_ID_clust1[prof_ID_clust2 %in% proID_Biotrans] # 3 in biotransform
c2_f <- prof_ID_clust1[prof_ID_clust2 %in% proID_Form] # 108 out of 116 in Form
c2_nc <- prof_ID_clust1[prof_ID_clust2 %in% proID_NoChange] # 5 in No Change

# cluster 3 intersects
c3_bt <- prof_ID_clust1[prof_ID_clust3 %in% proID_Biotrans] # 41 out of 125 in Biotrans
c3_f <- prof_ID_clust1[prof_ID_clust3 %in% proID_Form] # 34 out of 125 in Form
c3_nc <- prof_ID_clust1[prof_ID_clust3 %in% proID_NoChange] # 50 out of 125 in No Change

sum_base_data2 <- sum_base_data
sum_base_data2$clust <- rep(NA, nrow(sum_base_data2))
sum_base_data2$clust[sum_base_data2$profile_ID %in% prof_ID_clust1] <- "Cluster 1"
sum_base_data2$clust[sum_base_data2$profile_ID %in% prof_ID_clust2] <- "Cluster 2"
sum_base_data2$clust[sum_base_data2$profile_ID %in% prof_ID_clust3] <- "Cluster 3"
sum_base_data2 <- na.omit(sum_base_data2)

tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/cluster_mz_density_231025.tiff"), units="in", width=10, height=8, res=300)
ggplot(sum_base_data2, aes(x=avg_mz, fill = clust))+
  # geom_histogram(position = "identity",alpha = 0.3)+
  geom_density(alpha = 0.6)+
  # facet_wrap(~ Place)+
  xlab("Average m/z by profile ID")+
  ylab("Density")+
  ylim(0,0.005)+
  labs(fill = "Cluster")+
  theme(text = element_text(size=18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.32, 'cm'),
        axis.text.x = element_text(vjust = 1, hjust=1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        # axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray97"),
        panel.grid.minor = element_line(color = "gray97"),
        panel.border = element_rect(color = "grey", fill = NA))
dev.off()

tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/cluster_RT_density_231025.tiff"), units="in", width=10, height=8, res=300)
ggplot(sum_base_data2, aes(x=avg_RT/60, fill = clust))+
  # geom_histogram(position = "identity",alpha = 0.3)+
  geom_density(alpha = 0.6)+
  # facet_wrap(~ Place)+
  xlab("Average RT (min)")+
  ylab("Density")+
  ylim(0,0.15)+
  labs(fill = "Cluster")+
  theme(text = element_text(size=18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.32, 'cm'),
        axis.text.x = element_text(vjust = 1, hjust=1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        # axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray97"),
        panel.grid.minor = element_line(color = "gray97"),
        panel.border = element_rect(color = "grey", fill = NA))
dev.off()

tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/cluster_PA_density_231025.tiff"), units="in", width=10, height=8, res=300)
ggplot(sum_base_data2, aes(x=log10(avg_int), fill = clust))+
  # geom_histogram(position = "identity",alpha = 0.3)+
  geom_density(alpha = 0.6)+
  # facet_wrap(~ Place)+
  xlab("Average Log Peak Area")+
  ylab("Density")+
  ylim(0,1)+
  labs(fill = "Cluster")+
  theme(text = element_text(size=18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.32, 'cm'),
        axis.text.x = element_text(vjust = 1, hjust=1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        # axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray97"),
        panel.grid.minor = element_line(color = "gray97"),
        panel.border = element_rect(color = "grey", fill = NA))
dev.off()

#### END HEATMAP ANALYSIS #####



#### REACTION COUNTING ####
# Count the number of reactions observed for only persistent dataset
# Import functions defined in other R script
source("enviMass_reaction_counting_FUNCTIONS.R")
btrules <- read.csv("scholee_gulde_biotransformations.csv", header = TRUE)

# step 1: get INF profile IDs for each desired dataset
# set one: core higher in INF than EFF (cluster 1)
res <- pheatmap(t(new.plot[,c(3:ncol_2)]), scale = "row",cluster_cols = F, labels_col = t(new.plot$Date), show_rownames = F,
         annotation_col = annot_col2, cutree_rows = 3,clustering_method = "ward.D2", color = rampcols)
core.clust <- as.data.frame(cbind(t(new.plot[,c(3:ncol_2)]), cluster = cutree(res$tree_row, k = 3)))

colnames(INF_profiles)[1] <- c("profile_ID")
colnames(EFF_profiles)[1] <- c("profile_ID")
colnames(base_data)[8] <- c("profile_ID")
colnames(working_data)[1] <- c("profile_ID")
colnames(base_data_notrip)[1] <- c("profile_ID")

#### RUN ANALYSIS ON ALL PROFILES ####

# step 1 - predict
INF_pred <- mz_search(INF_profiles$profile_ID, INF_profiles$mean_mz, btrules)
colnames(INF_pred) <- c("prof_ID_INF", "prof_mass_INF",colnames(INF_pred)[3:21])

# step 2 - match
#ALL_ranged_matches <- mz_ranged_mass_function(EFF_profiles, INF_pred, 1)
#save(ALL_ranged_matches, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/ALL_ranged_matches.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/ALL_ranged_matches.RData"))

# write.csv(INF_pred, "INF_pred.csv")
# write.csv(EFF_profiles, "EFF_profiles.csv")
# write.csv(base_data_notrip, "base_data_notrip.csv")


# step 3 - count
# ALL_daily_counts <- daily_counts_function_v2(ALL_ranged_matches, base_data_notrip, btrules)
#save(ALL_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/ALL_daily_counts.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/ALL_daily_counts.RData"))

# # testing using base_data_notrip - same results...
# ALL_daily_counts_HUD2 <- daily_counts_function_v2(ALL_ranged_matches_HUD, base_data_notrip, btrules)
# daily_barplot_function(ALL_daily_counts_HUD2, yaxis = 2500)
# boxplot_function(ALL_daily_counts_HUD2)

# step 4 - plot
tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/ALL_daily_counts_231014.tiff"), units="in", width=10, height=8, res=300)
daily_barplot_function(ALL_daily_counts, yaxis = 10000)
dev.off()

tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/ALL_daily_counts_boxplot_231014.tiff"), units="in", width=10, height=8, res=300)
boxplot_function(ALL_daily_counts)
dev.off()


# step 5 - apply RT filter and sort matches accordingly
# ALL_RT_filter <- RT_detect_criteria(ALL_ranged_matches,working_data)
#save(ALL_RT_filter, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/ALL_RT_filter.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/ALL_RT_filter.RData"))

#ALL_RT_daily_counts <- daily_counts_function_v2(ALL_RT_filter[ALL_RT_filter$RT_match == TRUE,], base_data_notrip, btrules)
#save(ALL_RT_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/ALL_RT_daily_counts.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/ALL_RT_daily_counts.RData"))

tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/ALL_RT_daily_counts_RT_231010.tiff"), units="in", width=10, height=8, res=300)
daily_barplot_function(ALL_RT_daily_counts, yaxis = 10000)
dev.off()

tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/ALL_daily_counts_231010_RT_boxplot_231010.tiff"), units="in", width=10, height=8, res=300)
boxplot_function(ALL_RT_daily_counts, 1000)
dev.off()

# # write.csv(ALL_ranged_matches, "real_matches_231008.csv")
# ALL_RT_daily_counts_ratio1 <- daily_counts_function_v2(ALL_RT_filter[ALL_RT_filter$RT_match == TRUE & ALL_RT_filter$detect_ratio == 1,], base_data_notrip, btrules)
# save(ALL_RT_daily_counts_ratio1, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/ALL_RT_daily_counts_ratio1.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/ALL_RT_daily_counts_ratio1.RData"))
# 
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/ALL_RT_daily_counts_ratio1_231014.tiff"), units="in", width=10, height=8, res=300)
# daily_barplot_function(ALL_RT_daily_counts_ratio1, yaxis = 1000)
# dev.off()
# 
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/ALL_RT_daily_counts_ratio1_boxplot_231014.tiff"), units="in", width=10, height=8, res=300)
# boxplot_function(ALL_RT_daily_counts_ratio1)
# dev.off()


#--------------------------------------------------------------------------------------------------------------------------

#### try rarefaction - 10/20/2023 ####

# 1. randomly select 4,000 samples from predicted biotransformations
# set.seed makes sure you select the sample random samples every time you run the sample() function
set.seed(2347723)
# pick number of rows to randomly sample
n_sample.inf <- 8000
n_sample.eff <- 6000
# sample dataset
predict_INF_rare <- INF_pred[sample(1:n_sample.inf),]
search_EFF_rare <- EFF_profiles[sample(1:n_sample.eff),]

# run searching workflow

rare_ranged_matches2 <- mz_ranged_mass_function(search_EFF_rare, predict_INF_rare, 1)
save(rare_ranged_matches2, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/rare_ranged_matches2.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/rare_ranged_matches2.RData"))


colnames(base_data)[8] <- c("profile_ID")
rare_daily_counts2 <- daily_counts_function_v2(rare_ranged_matches2, base_data_notrip, btrules)
save(rare_daily_counts2, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/rare_daily_counts2.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/rare_daily_counts2.RData"))

colnames(working_data)[1] <- c("profile_ID")
rare_RT_filter2 <- RT_detect_criteria(rare_ranged_matches2,working_data)
save(rare_RT_filter2, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/rare_RT_filter2.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/rare_RT_filter2.RData"))

rare_RT_daily_counts2 <- daily_counts_function_v2(rare_RT_filter2[rare_RT_filter2$RT_match == TRUE,], base_data, btrules)
save(rare_RT_daily_counts2, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/rare_RT_daily_counts2.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/rare_RT_daily_counts2.RData"))

# try boxplot
tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/RARE_RT_dcbox_ALL_231020.tiff"), units="in", width=10, height=8, res=300)
boxplot_function(rare_RT_daily_counts2, 1000)
dev.off()


#### New analysis from Oct 8th 2023-------------------- changing INF predict and EFF search
n_core <- length(unique(x7$profile_ID[x7$type == "Core"]))
n_INFOnly <- length(unique(x7$profile_ID[x7$type == "INF_Only"]))
n_EFFOnly <- length(unique(x7$profile_ID[x7$type == "EFF_Only"]))
n_Inter <- length(unique(x7$profile_ID[x7$type == "Intermittent"]))
n_Form <- length(unique(x7$profile_ID[x7$type_new == "Form"]))
n_Biotrans <- length(unique(x7$profile_ID[x7$type_new == "Biotrans"]))
n_nochange <- length(unique(x7$profile_ID[x7$type_new == "No Change"]))

INFONly <- as.data.frame(unique(x7$profile_ID[x7$type == "INF_Only"]))
colnames(INFONly) <- c("profile_ID")

# step 1: get profile IDs for chemical features in INF_Only and Biotrans datasets
proID_INFonly <- as.numeric(unique(summ_peakarea$profile_ID[summ_peakarea$type_new == "INF Only"]))
proID_Biotrans <- as.numeric(unique(summ_peakarea$profile_ID[summ_peakarea$type_new == "Biotrans"]))
# proID_nochange <- as.numeric(unique(summ_peakarea$profile_ID[summ_peakarea$type_new == "No Change"]))

predictIDs <- c(proID_INFonly,proID_Biotrans)

# subset INF profiles based on categories to be used for prediction
predict_subprofs <- INF_profiles[INF_profiles$profile_ID %in% predictIDs,]
colnames(predict_subprofs)[1] <- c("profile_ID")
save(predict_subprofs, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/predict_subprofs.RData"))

# step 2: get profile IDs for chemical features in EFF only and Form datasets
proID_EFFonly <- as.numeric(unique(summ_peakarea$profile_ID[summ_peakarea$type_new == "EFF Only"]))
proID_Form <- as.numeric(unique(summ_peakarea$profile_ID[summ_peakarea$type_new == "Form"]))

searchIDs <- c(proID_EFFonly,proID_Form)

# subset EFF profiles based on categories to be searched for biotransformation products
search_subprofs <- EFF_profiles[EFF_profiles$profile_ID %in% searchIDs,]
colnames(search_subprofs)[1] <- c("profile_ID")
save(search_subprofs, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/search_subprofs.RData"))

# predict biotransformations
predict_overall <- mz_search(predict_subprofs$profile_ID, predict_subprofs$mean_mz, btrules)
colnames(predict_overall) <- c("prof_ID_INF", "prof_mass_INF",colnames(predict_overall)[3:21])
save(predict_overall, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/predict_overall.RData"))


# 
# #### analysis on "predict overall" predicted biotransformation masses - 10/21/2023 ####

# define a function within the apply function to count values meeting the criteria (in scan range) among each column in the dataframe
y3 <- as.data.frame(apply(predict_overall[,c(3:21)], 2, function(x) sum(x > 100.0000 & x < 999.9999)))
sum(y3)
# total possible is (7126*19 = 135394) but only 132394 meet the criteria here

# 
# #### END predict overall analysis ####


#### real analysis #### INF Only + Biotrans with EFF Only and Form ### IN MANUSCRIPT ####
#overall_ranged_matches <- mz_ranged_mass_function(search_subprofs, predict_overall, 1)
#save(overall_ranged_matches, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/overall_ranged_matches.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/overall_ranged_matches.RData"))

colnames(base_data)[8] <- c("profile_ID")
overall_daily_counts <- daily_counts_function_v2(overall_ranged_matches, base_data_notrip, btrules)
save(overall_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/overall_daily_counts.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/overall_daily_counts.RData"))

tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/INFbiotrans_EFFform_daily_counts_231010.tiff"), units="in", width=10, height=8, res=300)
daily_barplot_function(overall_daily_counts, yaxis = 2500)
dev.off()

# try boxplot
tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/INFbiotrans_EFFform_daily_counts_boxplot_231010.tiff"), units="in", width=10, height=8, res=300)
boxplot_function(overall_daily_counts)
dev.off()


#colnames(working_data)[1] <- c("profile_ID")
#overall_RT_filter <- RT_detect_criteria(overall_ranged_matches,working_data)
#save(overall_RT_filter, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/overall_RT_filter.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/overall_RT_filter.RData"))

#overall_RT_daily_counts <- daily_counts_function_v2(overall_RT_filter[overall_RT_filter$RT_match == TRUE,], base_data, btrules)
#save(overall_RT_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/overall_RT_daily_counts.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/overall_RT_daily_counts.RData"))

tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/INFbiotrans_EFFform_daily_counts_RT_231008.tiff"), units="in", width=10, height=8, res=300)
tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/INFbiotrans_EFFform_daily_counts_RT_231116.tiff"), units="in", width=10, height=8, res=300)
daily_barplot_function(overall_RT_daily_counts, yaxis = 2000)
dev.off()

# try boxplot
tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/INFbiotrans_EFFform_daily_counts_RT_boxplot_231118.tiff"), units="in", width=10, height=8, res=300)
boxplot_function(overall_RT_daily_counts,225)
dev.off()

# get total unique matches
total_unique_matches <- overall_RT_filter[overall_RT_filter$RT_match == TRUE,]

require(gridExtra)
plot2 <- daily_barplot_function(overall_RT_daily_counts, yaxis = 2000)
plot1 <- boxplot_function(overall_RT_daily_counts,225)

# BROAD REACTION COUNTS USING CATEGORIES - NEW FIGURE 3 AS OF OCT 17TH
tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/INFOnlyBiotransEFFForm_RT_rxncounts_231017.tiff"),
     units="in", width=16, height=8, res=300)
grid.arrange(plot1, plot2, ncol=2, widths = c(4.2/10,5.8/10))
dev.off()

# INFORMATION FOR TABLE S10 IN THE SI
stats <- as.data.frame(cbind(mean = rowMeans(overall_RT_daily_counts),sd = apply(overall_RT_daily_counts, 1, sd), median = apply(overall_RT_daily_counts, 1, median)))
stats$cov <- stats$sd/stats$mean

write.csv(stats, "overall_stats_231024.csv")

sum_bt_counts <- as.data.frame(cbind(sum = colSums(overall_RT_daily_counts), mean = colMeans(overall_RT_daily_counts),
                                     sd = apply(overall_RT_daily_counts, 2, sd), median = apply(overall_RT_daily_counts, 2, median)))
sum_bt_counts$Date <- c("Nov18","Nov19","Nov20","Nov21","Nov22",
                        "Nov23","Nov24","Nov25","Nov26","Nov27",
                        "Nov28","Nov29","Nov30","Dec1")


sum_bt_counts <- merge(sum_bt_counts, feat.counts.eff, by = "Date")
sum_bt_counts$Date <- factor(sum_bt_counts$Date, levels = c("Nov18","Nov19","Nov20","Nov21","Nov22",
                                                            "Nov23","Nov24","Nov25","Nov26","Nov27",
                                                            "Nov28","Nov29","Nov30","Dec1"))

colnames(sum_bt_counts)[7] <-c("n_feat_EFF")
sum_bt_counts$perc_bts <- sum_bt_counts$sum*100/sum_bt_counts$n_feat_EFF
write.csv(sum_bt_counts, "overall_stats_daily_231025.csv")

cor.test(sum_bt_counts$sum,sum_bt_counts$n_feat_EFF, method = "pearson")
lm(sum_bt_counts$sum ~ sum_bt_counts$n_feat_EFF)


# try using only cluster one for prediction biotransformations
prof_ID_clust1 <- as.numeric(rownames(core.clust[core.clust$cluster == 1,]))
clust1_subprofs <- INF_profiles[INF_profiles$profile_ID %in% prof_ID_clust1,]
colnames(clust1_subprofs)[1] <- c("profile_ID")

clust1_pred <- mz_search(clust1_subprofs$profile_ID, clust1_subprofs$mean_mz, btrules)
colnames(clust1_pred) <- c("prof_ID_INF", "prof_mass_INF",colnames(clust1_pred)[3:21])

#### analyze clust1_pred #### - 2005 total predicted biotransformatio products in analytical window
y4 <- as.data.frame(apply(clust1_pred[,c(3:21)], 2, function(x) sum(x > 100.0000 & x < 999.9999)))
sum(y4)


# clust1_ranged_matches <- mz_ranged_mass_function(search_subprofs, clust1_pred, 1)
# save(clust1_ranged_matches, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_ranged_matches.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_ranged_matches.RData"))

# clust1_daily_counts <- daily_counts_function_v2(clust1_ranged_matches, base_data, btrules)
# save(clust1_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_daily_counts.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_daily_counts.RData"))

tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/clust1_daily_counts_231008.tiff"), units="in", width=10, height=8, res=300)
daily_barplot_function(clust1_daily_counts, yaxis = 250)
dev.off()

tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/clust1_daily_counts_boxplot_231008.tiff"), units="in", width=10, height=8, res=300)
boxplot_function(clust1_daily_counts, yaxis = 50)
dev.off()

# clust1_RT_filter <- RT_detect_criteria(clust1_ranged_matches,working_data)
# save(clust1_RT_filter, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_RT_filter.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_RT_filter.RData"))

# clust1_RT_daily_counts <- daily_counts_function_v2(clust1_RT_filter[clust1_RT_filter$RT_match == TRUE,], base_data, btrules)
# save(clust1_RT_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_RT_daily_counts.RData"))
load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_RT_daily_counts.RData"))

tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/clust1_daily_counts_RT_231008.tiff"), units="in", width=10, height=8, res=300)
daily_barplot_function(clust1_RT_daily_counts, yaxis = 250)
dev.off()

tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/clust1_daily_counts_boxplot_RT_231008.tiff"), units="in", width=10, height=8, res=300)
boxplot_function(clust1_RT_daily_counts,30)
dev.off()


require(gridExtra)
plot4 <- daily_barplot_function(clust1_RT_daily_counts, yaxis = 200)
plot3 <- boxplot_function(clust1_RT_daily_counts,30)

# BROAD REACTION COUNTS USING CATEGORIES - NEW FIGURE 3 AS OF OCT 17TH
tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Clust1_EFFForm_RT_rxncounts_231024.tiff"),
     units="in", width=16, height=8, res=300)
grid.arrange(plot3, plot4, ncol=2, widths = c(4.2/10,5.8/10))
dev.off()

# data for table S11 in SI
stats_clust1 <- as.data.frame(cbind(mean = rowMeans(clust1_RT_daily_counts),sd = apply(clust1_RT_daily_counts, 1, sd), median = apply(clust1_RT_daily_counts, 1, median)))
stats_clust1$cov <- stats_clust1$sd/stats_clust1$mean
write.csv(stats_clust1, "clust1_rxn_stats_231024.csv")

clust1_bt_counts <- as.data.frame(cbind(sum = colSums(clust1_RT_daily_counts), mean = colMeans(clust1_RT_daily_counts),
                                     sd = apply(clust1_RT_daily_counts, 2, sd), median = apply(clust1_RT_daily_counts, 2, median)))
clust1_bt_counts$Date <- c("Nov18","Nov19","Nov20","Nov21","Nov22",
                        "Nov23","Nov24","Nov25","Nov26","Nov27",
                        "Nov28","Nov29","Nov30","Dec1")


clust1_bt_counts <- merge(clust1_bt_counts, feat.counts.eff, by = "Date")
clust1_bt_counts$Date <- factor(clust1_bt_counts$Date, levels = c("Nov18","Nov19","Nov20","Nov21","Nov22",
                                                            "Nov23","Nov24","Nov25","Nov26","Nov27",
                                                            "Nov28","Nov29","Nov30","Dec1"))

colnames(clust1_bt_counts)[7] <-c("n_feat_EFF")
clust1_bt_counts$perc_bts <- clust1_bt_counts$sum*100/clust1_bt_counts$n_feat_EFF
write.csv(clust1_bt_counts, "clust1_stats_daily_231025.csv")





# 
# 
# # # plot a profile?
# # profile1 <- as.data.frame(cbind(as.numeric(core.clust6[rownames(core.clust6) == 8362,])[1:28],date_nums))
# # profile1$feature <- rep(c("profile1"),28)
# # profile2 <- as.data.frame(cbind(as.numeric(core.clust6[rownames(core.clust6) == 6530,])[1:28],date_nums))
# # profile2$feature <- rep(c("profile2"),28)
# # date_nums <- rep(c("Nov18","Nov19","Nov20","Nov21","Nov22","Nov23","Nov24","Nov25","Nov26","Nov27","Nov28","Nov29","Nov30","Dec1"),2)
# # 
# # profile_data <- rbind(profile1, profile2)
# # colnames(profile_data) <- c("Intensity","Date","feature")
# # profile_data$Intensity <- as.numeric(profile_data$Intensity)
# # 
# # ggplot(profile_data, aes(x=Date, y=log10(Intensity), color=feature)) + 
# #   geom_point(size=6)
# 
# 
# # step 1: collect profile IDs of the groups you are interested in counting reactions for
# # cluster numbers are adjusted for figure naming
# prof_ID_clust1 <- as.numeric(rownames(core.clust[core.clust$cluster == 1,]))
# prof_ID_clust2 <- as.numeric(rownames(core.clust[core.clust$cluster == 4,]))
# prof_ID_clust3 <- as.numeric(rownames(core.clust[core.clust$cluster == 2,]))
# prof_ID_clust4 <- as.numeric(rownames(core.clust[core.clust$cluster == 3,]))
# prof_ID_infonly <- summ_type$profile_ID[summ_type$type == "INF Only"]
# 
# # step 2: subset INF profiles to only include features in desired dataset for prediction
# core_subprofs <- INF_profiles[INF_profiles$profileIDs %in% rownames(core.clust),]
# clust1_subprofs <- INF_profiles[INF_profiles$profileIDs %in% prof_ID_clust1,]
# clust2_subprofs <- INF_profiles[INF_profiles$profileIDs %in% prof_ID_clust2,]
# clust3_subprofs <- INF_profiles[INF_profiles$profileIDs %in% prof_ID_clust3,]
# clust4_subprofs <- INF_profiles[INF_profiles$profileIDs %in% prof_ID_clust4,]
# infonly_subprofs <- INF_profiles[INF_profiles$profileIDs %in% prof_ID_infonly,]
# 
# # step 3: make sure dataframes have correct column names
# colnames(core_subprofs)[1] <- c("profile_ID")
# colnames(clust1_subprofs)[1] <- c("profile_ID")
# colnames(clust2_subprofs)[1] <- c("profile_ID")
# colnames(clust3_subprofs)[1] <- c("profile_ID")
# colnames(clust4_subprofs)[1] <- c("profile_ID")
# colnames(infonly_subprofs)[1] <- c("profile_ID")
# colnames(EFF_profiles)[1] <- c("profile_ID")
# colnames(base_data)[8] <- c("profile_ID")
# 
# 
# 
# # step 4: apply function to generate predicted values - make sure profileIDs column name is profile_ID
# core_pred <- mz_search(core_subprofs$profile_ID, core_subprofs$mean_mz, btrules)
# colnames(core_pred) <- c("prof_ID_INF", "prof_mass_INF",colnames(core_pred)[3:21])
# 
# clust1_pred <- mz_search(clust1_subprofs$profile_ID, clust1_subprofs$mean_mz, btrules)
# colnames(clust1_pred) <- c("prof_ID_INF", "prof_mass_INF",colnames(clust1_pred)[3:21])
# 
# clust2_pred <- mz_search(clust2_subprofs$profile_ID, clust2_subprofs$mean_mz, btrules)
# colnames(clust2_pred) <- c("prof_ID_INF", "prof_mass_INF",colnames(clust2_pred)[3:21])
# 
# clust3_pred <- mz_search(clust3_subprofs$profile_ID, clust3_subprofs$mean_mz, btrules)
# colnames(clust3_pred) <- c("prof_ID_INF", "prof_mass_INF",colnames(clust3_pred)[3:21])
# 
# clust4_pred <- mz_search(clust4_subprofs$profile_ID, clust4_subprofs$mean_mz, btrules)
# colnames(clust4_pred) <- c("prof_ID_INF", "prof_mass_INF",colnames(clust4_pred)[3:21])
# 
# infonly_pred <- mz_search(infonly_subprofs$profile_ID, infonly_subprofs$mean_mz, btrules)
# colnames(infonly_pred) <- c("prof_ID_INF", "prof_mass_INF",colnames(infonly_pred)[3:21])
# 
# 
# # step 5: find matched with 1 ppm range
# core_ranged_matches <- mz_ranged_mass_function(EFF_profiles, core_pred,1)
# clust1_ranged_matches <- mz_ranged_mass_function(EFF_profiles, clust1_pred, 1)
# clust2_ranged_matches <- mz_ranged_mass_function(EFF_profiles, clust2_pred, 1)
# clust3_ranged_matches <- mz_ranged_mass_function(EFF_profiles, clust3_pred, 1)
# clust4_ranged_matches <- mz_ranged_mass_function(EFF_profiles, clust4_pred, 1)
# 
# 
# # infonly_ranged_matches <- mz_ranged_mass_function(EFF_profiles, infonly_pred, 1)
# 
# # step 6: count daily matches
# core_daily_counts <- daily_counts_function_v2(core_ranged_matches, base_data, btrules) 
# clust1_daily_counts <- daily_counts_function_v2(clust1_ranged_matches, base_data, btrules)
# clust2_daily_counts <- daily_counts_function_v2(clust2_ranged_matches, base_data, btrules)
# clust3_daily_counts <- daily_counts_function_v2(clust3_ranged_matches, base_data, btrules)
# clust4_daily_counts <- daily_counts_function_v2(clust4_ranged_matches, base_data, btrules)
# 
# # INFonly_daily_counts <- daily_counts_function_v2(infonly_ranged_matches, base_data, btrules)
# 
# #save all generated dataframes for easy loading
# # matched masses
# save(core_ranged_matches, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/core_ranged_matches.RData"))
# save(clust1_ranged_matches, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_ranged_matches.RData"))
# save(clust2_ranged_matches, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust2_ranged_matches.RData"))
# save(clust3_ranged_matches, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust3_ranged_matches.RData"))
# save(clust4_ranged_matches, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust4_ranged_matches.RData"))
# # save(infonly_ranged_matches, file = "~/Project 3/All_objectives_P3WWMetabo/Robjects_toload/infonly_ranged_matches.RData")
# 
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/core_ranged_matches.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_ranged_matches.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust2_ranged_matches.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust3_ranged_matches.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust4_ranged_matches.RData"))
# # load("~/Project 3/All_objectives_P3WWMetabo/Robjects_toload/infonly_ranged_matches.RData")
# 
# 
# # daily count dataframes
# save(core_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/core_daily_counts.RData"))
# save(clust1_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_daily_counts.RData"))
# save(clust2_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust2_daily_counts.RData"))
# save(clust3_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust3_daily_counts.RData"))
# save(clust4_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust4_daily_counts.RData"))
# # save(INFonly_daily_counts, file = "~/Project 3/All_objectives_P3WWMetabo/Robjects_toload/INFonly_daily_counts.RData")
# 
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/core_daily_counts.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_daily_counts.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust2_daily_counts.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust3_daily_counts.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust4_daily_counts.RData"))
# # load("~/Project 3/All_objectives_P3WWMetabo/Robjects_toload/INFonly_daily_counts.RData")
# 
# dev.off()
# # step 7: plot stuff
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/core_rxns_230930.tiff"), units="in", width=10, height=8, res=300)
# daily_barplot_function(core_daily_counts, yaxis = 2000)
# dev.off()
# 
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/clust1_rxns_230930.tiff"), units="in", width=10, height=8, res=300)
# daily_barplot_function(clust1_daily_counts, yaxis = 500)
# dev.off()
# 
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/clust2_rxns_230930.tiff"), units="in", width=10, height=8, res=300)
# daily_barplot_function(clust2_daily_counts, yaxis = 500)
# dev.off()
# 
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/clust3_rxns_230930.tiff"), units="in", width=10, height=8, res=300)
# daily_barplot_function(clust3_daily_counts, yaxis = 500)
# dev.off()
# 
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/clust4_rxns_230930.tiff"), units="in", width=10, height=8, res=300)
# daily_barplot_function(clust4_daily_counts, yaxis = 500)
# dev.off()
# 
# # tiff("INFonly_rxns_230929.tiff", units="in", width=10, height=8, res=300)
# # daily_barplot_function(INFonly_daily_counts, yaxis = 5000)
# # dev.off()
# 
# clust1_daily_counts_sumstats <- clust1_daily_counts
# clust1_daily_counts_sumstats$mean <- rowMeans(clust1_daily_counts)
# clust1_daily_counts_sumstats$sd <- apply(clust1_daily_counts, 1, sd)
# clust1_daily_counts_sumstats$CV <- clust1_daily_counts_sumstats$sd/clust1_daily_counts_sumstats$mean
# 
# clust2_daily_counts_sumstats <- clust2_daily_counts
# clust2_daily_counts_sumstats$mean <- rowMeans(clust2_daily_counts)
# clust2_daily_counts_sumstats$sd <- apply(clust2_daily_counts, 1, sd)
# clust2_daily_counts_sumstats$CV <- clust2_daily_counts_sumstats$sd/clust2_daily_counts_sumstats$mean
# 
# clust3_daily_counts_sumstats <- clust3_daily_counts
# clust3_daily_counts_sumstats$mean <- rowMeans(clust3_daily_counts)
# clust3_daily_counts_sumstats$sd <- apply(clust3_daily_counts, 1, sd)
# clust3_daily_counts_sumstats$CV <- clust3_daily_counts_sumstats$sd/clust3_daily_counts_sumstats$mean
# 
# clust4_daily_counts_sumstats <- clust4_daily_counts
# clust4_daily_counts_sumstats$mean <- rowMeans(clust4_daily_counts)
# clust4_daily_counts_sumstats$sd <- apply(clust4_daily_counts, 1, sd)
# clust4_daily_counts_sumstats$CV <- clust4_daily_counts_sumstats$sd/clust4_daily_counts_sumstats$mean
# 
# #### STARTED HERE ON SATURDAY 9/30/2023
# 
# 
# # ADD RT FILTER FOR REACTION COUNTING
# colnames(working_data)[1] <- c("profile_ID")
# core_RT_filter <- RT_detect_criteria(core_ranged_matches, working_data)
# clust1_RT_filter <- RT_detect_criteria(clust1_ranged_matches,working_data)
# clust2_RT_filter <- RT_detect_criteria(clust2_ranged_matches,working_data)
# clust3_RT_filter <- RT_detect_criteria(clust3_ranged_matches,working_data)
# clust4_RT_filter <- RT_detect_criteria(clust4_ranged_matches,working_data)
# 
# core_RT_daily_counts <- daily_counts_function_v2(core_RT_filter[core_RT_filter$RT_match == TRUE,], base_data, btrules)
# clust1_RT_daily_counts <- daily_counts_function_v2(clust1_RT_filter[clust1_RT_filter$RT_match == TRUE,], base_data, btrules)
# clust2_RT_daily_counts <- daily_counts_function_v2(clust2_RT_filter[clust2_RT_filter$RT_match == TRUE,], base_data, btrules)
# clust3_RT_daily_counts <- daily_counts_function_v2(clust3_RT_filter[clust3_RT_filter$RT_match == TRUE,], base_data, btrules)
# clust4_RT_daily_counts <- daily_counts_function_v2(clust4_RT_filter[clust4_RT_filter$RT_match == TRUE,], base_data, btrules)
# 
# 
# # Load data if this script has already been run before
# save(core_RT_filter, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/core_RT_filter.RData"))
# save(clust1_RT_filter, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_RT_filter.RData"))
# save(clust2_RT_filter, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust2_RT_filter.RData"))
# save(clust3_RT_filter, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust3_RT_filter.RData"))
# save(clust4_RT_filter, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust4_RT_filter.RData"))
# 
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/core_RT_filter.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_RT_filter.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust2_RT_filter.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust3_RT_filter.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust4_RT_filter.RData"))
# 
# save(core_RT_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/core_RT_daily_counts.RData"))
# save(clust1_RT_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_RT_daily_counts.RData"))
# save(clust2_RT_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust2_RT_daily_counts.RData"))
# save(clust3_RT_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust3_RT_daily_counts.RData"))
# save(clust4_RT_daily_counts, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust4_RT_daily_counts.RData"))
# 
# 
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/core_RT_daily_counts.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust1_RT_daily_counts.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust2_RT_daily_counts.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust3_RT_daily_counts.RData"))
# load(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/Robjects_toload/clust4_RT_daily_counts.RData"))
# 
# 
# daily_barplot_function(core_RT_daily_counts, yaxis = 2000)
# daily_barplot_function(clust1_RT_daily_counts, yaxis = 500)
# daily_barplot_function(clust2_RT_daily_counts, yaxis = 500)
# daily_barplot_function(clust3_RT_daily_counts, yaxis = 500)
# daily_barplot_function(clust4_RT_daily_counts, yaxis = 500)
# 
# require(gridExtra)
# plot1 <- daily_barplot_function_left_top(clust1_RT_daily_counts, yaxis = 500)
# plot2 <- daily_barplot_function_right_top(clust2_RT_daily_counts, yaxis = 500)
# plot3 <- daily_barplot_function_left(clust3_RT_daily_counts, yaxis = 500)
# plot4 <- daily_barplot_function_right(clust4_RT_daily_counts, yaxis = 500)
# 
# #clust1 and clust2 reaction counts
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/reaction_counts_4clusters_RT.tiff"), units="in", width=10, height=8, res=300)
# grid.arrange(plot1, plot2,plot3, plot4, ncol=2, widths = c(4.5/10,5.5/10))
# dev.off()
# 
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/reaction_counts_core_230930_RT.tiff"), units="in", width=10, height=8, res=300)
# boxplot_function(core_RT_daily_counts)
# dev.off()
# 
# 
# core_RT_daily_counts_sumstats <- core_RT_daily_counts
# core_RT_daily_counts_sumstats$mean <- rowMeans(core_RT_daily_counts)
# core_RT_daily_counts_sumstats$sd <- apply(core_RT_daily_counts, 1, sd)
# core_RT_daily_counts_sumstats$CV <- core_RT_daily_counts_sumstats$sd/core_RT_daily_counts_sumstats$mean
# 
# write.csv(core_RT_daily_counts_sumstats, "core_RT_sumstats.csv")
# 
# clust1_RT_daily_counts_sumstats <- clust1_RT_daily_counts
# clust1_RT_daily_counts_sumstats$mean <- rowMeans(clust1_RT_daily_counts)
# clust1_RT_daily_counts_sumstats$sd <- apply(clust1_RT_daily_counts, 1, sd)
# clust1_RT_daily_counts_sumstats$CV <- clust1_RT_daily_counts_sumstats$sd/clust1_RT_daily_counts_sumstats$mean
# 
# clust2_RT_daily_counts_sumstats <- clust2_RT_daily_counts
# clust2_RT_daily_counts_sumstats$mean <- rowMeans(clust2_RT_daily_counts)
# clust2_RT_daily_counts_sumstats$sd <- apply(clust2_RT_daily_counts, 1, sd)
# clust2_RT_daily_counts_sumstats$CV <- clust2_RT_daily_counts_sumstats$sd/clust2_RT_daily_counts_sumstats$mean
# 
# clust3_RT_daily_counts_sumstats <- clust3_RT_daily_counts
# clust3_RT_daily_counts_sumstats$mean <- rowMeans(clust3_RT_daily_counts)
# clust3_RT_daily_counts_sumstats$sd <- apply(clust3_RT_daily_counts, 1, sd)
# clust3_RT_daily_counts_sumstats$CV <- clust3_RT_daily_counts_sumstats$sd/clust3_RT_daily_counts_sumstats$mean
# 
# clust4_RT_daily_counts_sumstats <- clust4_RT_daily_counts
# clust4_RT_daily_counts_sumstats$mean <- rowMeans(clust4_RT_daily_counts)
# clust4_RT_daily_counts_sumstats$sd <- apply(clust4_RT_daily_counts, 1, sd)
# clust4_RT_daily_counts_sumstats$CV <- clust4_RT_daily_counts_sumstats$sd/clust4_RT_daily_counts_sumstats$mean
# 
# write.csv(clust1_RT_daily_counts_sumstats, "clust1_RT_daily_counts_sumstats.csv")
# write.csv(clust2_RT_daily_counts_sumstats, "clust2_RT_daily_counts_sumstats.csv")
# write.csv(clust3_RT_daily_counts_sumstats, "clust3_RT_daily_counts_sumstats.csv")
# write.csv(clust4_RT_daily_counts_sumstats, "clust4_RT_daily_counts_sumstats.csv")
# 
# plot4 <- boxplot_function(clust1_RT_daily_counts)
# plot5 <- boxplot_function(clust2_RT_daily_counts)
# plot6 <- boxplot_function(clust3_RT_daily_counts)
# plot7 <- boxplot_function(clust4_RT_daily_counts)
# 
# grid.arrange(plot4, plot5,plot6, plot7, ncol=2, widths = c(4.5/10,5.5/10))
# dev.off()
# 
# m1 <-melt(cbind(clust1_RT_daily_counts, rownames(clust1_RT_daily_counts)))
# colnames(m1) <-c("Reaction","Date","Count")
# m1$Reaction <- factor(m1$Reaction, levels=c("acetylation","hydroxylation","dihydroxylation","demethylation","deethylation","dehydrogenation","hydrogenation","dehydration","chlorine reduction",
#                                                             "deacetylation","glucuronidation","deglucuronidation","sulfonation","desulfonation","formylation","acylation","succinylation","fumarylation","malonylation"))
# m1$Cluster <- rep("Cluster 1", nrow(m1))
# 
# m2 <-melt(cbind(clust2_RT_daily_counts, rownames(clust2_RT_daily_counts)))
# colnames(m2) <-c("Reaction","Date","Count")
# m2$Reaction <- factor(m2$Reaction, levels=c("acetylation","hydroxylation","dihydroxylation","demethylation","deethylation","dehydrogenation","hydrogenation","dehydration","chlorine reduction",
#                                             "deacetylation","glucuronidation","deglucuronidation","sulfonation","desulfonation","formylation","acylation","succinylation","fumarylation","malonylation"))
# m2$Cluster <- rep("Cluster 2", nrow(m2))
# 
# m3 <-melt(cbind(clust3_RT_daily_counts, rownames(clust3_RT_daily_counts)))
# colnames(m3) <-c("Reaction","Date","Count")
# m3$Reaction <- factor(m3$Reaction, levels=c("acetylation","hydroxylation","dihydroxylation","demethylation","deethylation","dehydrogenation","hydrogenation","dehydration","chlorine reduction",
#                                             "deacetylation","glucuronidation","deglucuronidation","sulfonation","desulfonation","formylation","acylation","succinylation","fumarylation","malonylation"))
# m3$Cluster <- rep("Cluster 3", nrow(m3))
# 
# m4 <-melt(cbind(clust4_RT_daily_counts, rownames(clust4_RT_daily_counts)))
# colnames(m4) <-c("Reaction","Date","Count")
# m4$Reaction <- factor(m4$Reaction, levels=c("acetylation","hydroxylation","dihydroxylation","demethylation","deethylation","dehydrogenation","hydrogenation","dehydration","chlorine reduction",
#                                             "deacetylation","glucuronidation","deglucuronidation","sulfonation","desulfonation","formylation","acylation","succinylation","fumarylation","malonylation"))
# m4$Cluster <- rep("Cluster 4", nrow(m4))
# 
# all.clust.long <- rbind(m1,m2,m3,m4)
# 
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/boxplots_rxns_FIG3_230930_RT.tiff"), units="in", width=12, height=10, res=300)
# ggplot(all.clust.long, aes(x=reorder(Reaction, -Count), y=Count, fill = Reaction)) + 
#   geom_boxplot()+
#   facet_wrap(~Cluster)+
#   scale_fill_manual(values = c("#45D6FF","#FFC000","#C00000","#92D050","#FF83FA","#9575cd","#ff9797",
#                                "#ba632f","#2e7d32","#1c3ffd","#020873","#fa5b0f","#EE1289","#551A8B",
#                                "#99b3ff","#990099","#47d1d1","#e6e600","#c2c2a3"))+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         legend.key.size = unit(0.95, 'cm'),
#         # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()
# 
# # try changing variables in boxplot
# 
# # ggplot(all.clust.long, aes(x=reorder(Reaction, -Count), y=Count, fill = Reaction, color = Cluster)) + 
# #   geom_boxplot()+
# #   # facet_wrap(~Cluster)+
# #   scale_fill_manual(values = c("#45D6FF","#FFC000","#C00000","#92D050","#FF83FA","#9575cd","#ff9797",
# #                                "#ba632f","#2e7d32","#1c3ffd","#020873","#fa5b0f","#EE1289","#551A8B",
# #                                "#99b3ff","#990099","#47d1d1","#e6e600","#c2c2a3"))+
# #   theme(text = element_text(size=18),
# #         legend.text = element_text(size = 12),
# #         legend.title = element_text(size = 14),
# #         legend.key.size = unit(0.95, 'cm'),
# #         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
# #         # axis.text.x = element_blank(),
# #         axis.ticks.x = element_blank(),
# #         axis.text.y = element_text(size = 16),
# #         axis.title.x = element_blank(),
# #         panel.background = element_blank(),
# #         panel.grid.major = element_line(color = "gray97"),
# #         panel.grid.minor = element_line(color = "gray97"),
# #         panel.border = element_rect(color = "grey", fill = NA))
# 
# 
# all.clust.long <- rbind(m1,m2,m3,m4)
# 
# all.clust.long$Date <- with(all.clust.long, ifelse(Date == "2020-11-18","Nov18",
#                                                ifelse(Date == "2020-11-19","Nov19",
#                                                       ifelse(Date == "2020-11-20","Nov20",
#                                                              ifelse(Date == "2020-11-21","Nov21",
#                                                                     ifelse(Date == "2020-11-22","Nov22",
#                                                                            ifelse(Date == "2020-11-23","Nov23",
#                                                                                   ifelse(Date == "2020-11-24","Nov24",
#                                                                                          ifelse(Date == "2020-11-25","Nov25",
#                                                                                                 ifelse(Date == "2020-11-26","Nov26",
#                                                                                                        ifelse(Date == "2020-11-27","Nov27",
#                                                                                                               ifelse(Date == "2020-11-28","Nov28",
#                                                                                                                      ifelse(Date == "2020-11-29","Nov29",
#                                                                                                                             ifelse(Date == "2020-11-30","Nov30",
#                                                                                                                                    ifelse(Date == "2020-12-01","Dec1","x")))))))))))))))
# all.clust.long$Date <- factor(all.clust.long$Date, levels = c("Nov18","Nov19","Nov20","Nov21","Nov22",
#                                                           "Nov23","Nov24","Nov25","Nov26","Nov27",
#                                                           "Nov28","Nov29","Nov30","Dec1"))
# 
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/barplot_rxns_FIG4_230930_RT.tiff"), units="in", width=12, height=10, res=300)
# ggplot(data=all.clust.long, aes(x=Date, y=Count, fill = Reaction))+
#   geom_bar(stat="identity", size = 10)+
#   facet_wrap(~Cluster)+
#   scale_fill_manual(values = c("#45D6FF","#FFC000","#C00000","#92D050","#FF83FA","#9575cd","#ff9797",
#                                "#ba632f","#2e7d32","#1c3ffd","#020873","#fa5b0f","#EE1289","#551A8B",
#                                "#99b3ff","#990099","#47d1d1","#e6e600","#c2c2a3"))+
#   ylab("Count")+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         legend.key.size = unit(0.9, 'cm'),
#         # legend.position = "none",
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 18),
#         axis.title.y = element_text(size = 16),
#         axis.ticks.y = element_blank(),
#         axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()
# 
# 
# 
# 
# # try relative abundance plot for features intensities using based_data_notrip_type == core - or base_data_notrip
# 
# relabund_int <- base_data_notrip_type %>%
#   group_by(Date, Place) %>%
#   summarise(total = sum(mean_int))
# 
# relabund_int2 <- merge(base_data_notrip_type,relabund_int, by = c("Date","Place")) 
# relabund_int2$relabund <- (relabund_int2$mean_int/relabund_int2$total)*100
# 
# x5 <- relabund_int2[,c(1:3,7,9)]
# 
# relabund_int2$profile_ID <- as.factor(relabund_int2$profile_ID)
# 
# # relabund.m <- relabund_int2 %>% pivot_longer(cols = c("Date","Place","profile_ID","relabund"))
# 
# 
# ggplot(data=x5, aes(x=Date, fill = profile_ID))+
#   geom_bar()
# 
# x5$profile_ID <- as.factor(x5$profile_ID)
# 
# # SHOW THE RELATIVE ABUNDANCE OF EACH FEATURE BY TYPE IN % ##############################################
# 
# x5$type <- factor(x5$type,levels = c("INF Only","EFF Only","Intermittant","Core"))
# x5$Place <- factor(x5$Place, levels = c("INF","EFF"))
# colnames(x5) <- c("Date","Place","profile_ID","Type","Relative_Abundance")
# 
# x5$Date <- with(x5, ifelse(Date == "2020-11-18","Nov18",
#                                             ifelse(Date == "2020-11-19","Nov19",
#                                                   ifelse(Date == "2020-11-20","Nov20",
#                                                             ifelse(Date == "2020-11-21","Nov21",
#                                                                     ifelse(Date == "2020-11-22","Nov22",
#                                                                             ifelse(Date == "2020-11-23","Nov23",
#                                                                                   ifelse(Date == "2020-11-24","Nov24",
#                                                                                           ifelse(Date == "2020-11-25","Nov25",
#                                                                                                   ifelse(Date == "2020-11-26","Nov26",
#                                                                                                            ifelse(Date == "2020-11-27","Nov27",
#                                                                                                                   ifelse(Date == "2020-11-28","Nov28",
#                                                                                                                          ifelse(Date == "2020-11-29","Nov29",
#                                                                                                                                 ifelse(Date == "2020-11-30","Nov30",
#                                                                                                                                        ifelse(Date == "2020-12-01","Dec1","x")))))))))))))))
# x5$Date <- factor(x5$Date, levels = c("Nov18","Nov19","Nov20","Nov21","Nov22",
#                                                               "Nov23","Nov24","Nov25","Nov26","Nov27",
#                                                               "Nov28","Nov29","Nov30","Dec1"))
# 
# tiff(file.path("~/Project 3/All_objectives_P3WWMetabo/",results_folder,"/featureABUNDANCEbytypebyday_FIG1.1_230930.tiff"), units="in", width=14, height=8, res=300)
# ggplot(x5, aes(x = Date, y = Relative_Abundance, fill = Type))+ 
#   geom_bar(stat = "identity")+
#   facet_wrap(~Place)+
#   ylab("Relative Abundance (%)")+
#   scale_fill_manual(values = c("#FF6666","#00BFC4","#A568D2","orange"))+
#   theme(text = element_text(size=18),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.32, 'cm'),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "gray97"),
#         panel.grid.minor = element_line(color = "gray97"),
#         panel.border = element_rect(color = "grey", fill = NA))
# dev.off()
# #########################################################################################################
# 
# # # try with core dataset
# # base_data_core <- base_data_notrip_type[base_data_notrip_type$type == "Core",]
# # 
# # relabund_core <- base_data_core %>%
# #   group_by(Date, Place) %>%
# #   summarise(total = sum(mean_int))
# # 
# # base_data_core  <- merge(base_data_core,relabund_core, by = c("Date","Place")) 
# # base_data_core$relabund <- (base_data_core$mean_int/base_data_core$total)*100
# # 
# # xcore <- base_data_core[,c(1:3,7,9)]
# # 
# # base_data_core$profile_ID <- as.factor(base_data_core$profile_ID)
# # 
# # ggplot(xcore, aes(x = Date, y = relabund, fill = profile_ID))+ 
# #   geom_bar(stat = "identity")+
# #   facet_wrap(~Place)
# # 
# # 
# # xcore$profile_ID <- as.factor(xcore$profile_ID)
# # 
# # ggplot(xcore, aes(x = Date, y = relabund, fill = type))+ 
# #   geom_bar(stat = "identity")+
#   # facet_wrap(~Place)
# 
