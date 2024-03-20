library(knitr)
library(ggplot2)
library(devtools)
library(tibble)
library(xlsx)
library(plotly)
library(rsconnect)
library(webshot)
library(xml2)
library(enviMass)
library(enviPat)
library(ggplot2)
library(dplyr)
library(rcdk)
library(RMassBank)
library(reshape2)

##### START OF FUNCTION DEFINITONS ####--------------------------------------------------------------------------------

# function for generating predicted masses from INF profiles
# input: 
# profile_ID_column - column of profile IDs from INF profiles
# mean_mz_column - column of mean mz values from INF profiles
# bts - dataframe of listed btrules from literature

# example:
# INF_predicted <- mz_search(INF_search$profile_ID, INF_search$mean_mz, btrules)

mz_search <- function(profile_ID_column, mean_mz_column, bts){
  df <- as.data.frame(cbind(profile_ID_column, mean_mz_column))
  colnames(df) <-c("prof_ID","prof_mass")
  df$hydrox <- df$prof_mass + bts$mass_diff[bts$biotrans == "hydroxylation"]
  df$demeth <- df$prof_mass - bts$mass_diff[bts$biotrans == "demethylation"]
  df$deeth <- df$prof_mass - bts$mass_diff[bts$biotrans == "deethylation"]
  df$dehydro <- df$prof_mass - bts$mass_diff[bts$biotrans == "dehydrogenation"]
  df$hydro <- df$prof_mass + bts$mass_diff[bts$biotrans == "hydrogenation"]
  df$dehydra <- df$prof_mass - bts$mass_diff[bts$biotrans == "dehydration"]
  df$cl_red <- df$prof_mass - bts$mass_diff[bts$biotrans == "chlorine reduction"]
  df$acetyl <- df$prof_mass + bts$mass_diff[bts$biotrans == "acetylation"]
  df$deacetyl <- df$prof_mass - bts$mass_diff[bts$biotrans == "deacetylation"]
  df$gluco <- df$prof_mass + bts$mass_diff[bts$biotrans == "glucuronidation"]
  df$degluco <- df$prof_mass - bts$mass_diff[bts$biotrans == "deglucuronidation"]
  df$sulfo <- df$prof_mass + bts$mass_diff[bts$biotrans == "sulfonation"]
  df$desulfo <- df$prof_mass - bts$mass_diff[bts$biotrans == "desulfonation"]
  df$dihydrox <- df$prof_mass + bts$mass_diff[bts$biotrans == "dihydroxylation"]
  df$formyl <- df$prof_mass + bts$mass_diff[bts$biotrans == "formylation"]
  df$acyl <- df$prof_mass + bts$mass_diff[bts$biotrans == "acylation"]
  df$succ <- df$prof_mass + bts$mass_diff[bts$biotrans == "succinylation"]
  df$fumar <- df$prof_mass + bts$mass_diff[bts$biotrans == "fumarylation"]
  df$malo <- df$prof_mass + bts$mass_diff[bts$biotrans == "malonylation"]
  return(df)
}



##### Matching functions ####

# 1. reaction searching using exact masses - mz_exact_mass_function()
# input:
# measured - a dataframe containing unique Effluent profile IDs as a subset of working_data generated from profileList_pos enviMass output
# prediction - a datframe containing predicted masses of bt products using the mass difference approach and INF profile lists

mz_exact_mass_function <- function(measured, predicted){
  dfList <- list()  ## create empty list
  
  for(i in 1:nrow(measured)){
    EFF_ID <- measured$profile_ID[i]
    mz <- round(measured$mean_mz[i],4)
    matches <- round(predicted,4) %>% filter_all(any_vars(. %in% c(mz)))
    
    # only keep entries where mz value matches EFF mz that we are searching for
    matches[,c(3:21)][matches[,c(3:21)] != mz] <- NA
    matches <- matches[rowSums(matches[,c(3:21)], na.rm = T) != 0,] # remove rows that contain all NA values (no matches)
    
    EFF_Data <- cbind(EFF_ID = rep(EFF_ID,nrow(matches)), mz_EFF = rep(mz, nrow(matches)))
    dfList[[i]] <- cbind(matches, EFF_Data)
    print(i)
  }
  
  #---turn matching results into a dataframe---------
  total <- as.data.frame(do.call(rbind, dfList))
  
  # remove mz values that were matched because of matching profile IDs
  total <- subset(total, prof_ID_INF != EFF_ID)
  
  return(total)
}

# 2. reaction searching using a range of masses
# same as exact mass, but now with "mass_acc" option to set for adding a mass range for detected bt products
mz_ranged_mass_function <- function(measured, predicted, mass_acc){
  # measured <- EFF_profiles
  # predicted <- INF_predicted
  # mass_acc <- 1
  
  dfList <- list()  ## create empty list
  # mass_acc <- 5 # ppm
  
  for(i in 1:nrow(measured)){
    # i = 1
    EFF_ID <- measured$profile_ID[i]
    mz <- round(measured$mean_mz[i],4)
    mz_l <- round(mz - ppm(mz,mass_acc,p = TRUE),4)
    mz_h <- round(mz + ppm(mz,mass_acc,p = TRUE),4)
    range <- seq(mz_l, mz_h, by = 0.0001)
    matches <- round(predicted,4) %>% filter_all(any_vars(. %in% range))
    
    # only keep entries where mz value matches EFF mz that we are searching for
    matches[,c(3:21)][matches[,c(3:21)] < mz_l | matches[,c(3:21)] > mz_h] <- NA
    matches <- matches[rowSums(matches[,c(3:21)], na.rm = T) != 0,] # remove rows that contain all NA values (no matches)
      
    EFF_Data <- cbind(EFF_ID = rep(EFF_ID,nrow(matches)), mz_EFF = rep(mz, nrow(matches)))
    dfList[[i]] <- cbind(matches, EFF_Data)
    print(i)
    
  }
  #---turn matching results into a dataframe---------
  total <- as.data.frame(do.call(rbind, dfList))
  
  # remove mz values that were matched because of matching profile IDs
  total <- subset(total, prof_ID_INF != EFF_ID)
  
  return(total)
  
}

##### Counting function ####

# Description
# input:
# matches - output from mz_XX_mass_function containing mz values where there are matched INF/EFF IDs
# sampleinfo - called "base_data" in other scripts - dataframe containing information on the SPECIFIC DAYS on which each profile ID was detected
# bts - a dataframe containing 19 biotransformation rules from literature

daily_counts_function <- function(matches, sampleinfo, bts){
  
  # here the variable "sampleinfo" is the "working_data" datframe with different dates that profiles are detected
  
  # Count the number of matches per day
  Date <-c("2020-11-18","2020-11-19","2020-11-20","2020-11-21",
           "2020-11-22","2020-11-23","2020-11-24","2020-11-25","2020-11-26","2020-11-27",
           "2020-11-28","2020-11-29","2020-11-30","2020-12-01")
  
  # dataframe to store the results of match counting per day
  daily_counts <- data.frame(matrix(nrow = 19, ncol = 14))
  rownames(daily_counts) <- bts$biotrans
  colnames(daily_counts) <- Date
  daily_counts[is.na(daily_counts)] <- 0
  
  # pairs up columns in total to rows in "daily_counts"
  pairs <- cbind(c(1:nrow(daily_counts)),c(3:21))
  
  
  for(h in 1:nrow(daily_counts)){
    df.counts <- matches
    info <- sampleinfo
    col_num = pairs[h,][2] # this is the column number corresponding to the reaction type in the "total" data frame
    row_num = pairs[h,][1] # this is the row number corresponding to the reaction type in the "daily counts" results accounting
    
    # subsets whole dataframe to only include rows where a column entry is TRUE in the specified reaction column
    dat <- df.counts %>% 
      filter(df.counts[,col_num]!=0)
    
    num_exp <- nrow(dat)
    
    print(num_exp)
    
    for(i in 1:nrow(dat)){
      # selected profile IDs associated with these TRUE hits from "dat" dataset
      profile_ID_EFF <- dat$EFF_ID[i]
      profiles_EFF <- info[info$profile_ID %in% profile_ID_EFF,]
      profiles_EFF <- profiles_EFF[profiles_EFF$Place == "EFF",]
      
      profile_ID_INF <- dat$prof_ID_INF[i]
      profiles_INF <- info[info$profile_ID %in% profile_ID_INF,]
      profiles_INF <- profiles_INF[profiles_INF$Place == "INF",]
      
      # only counts as a reaction when the parent and TP are detected on the same date (merge INF and EFF profiles is an inner join and only keeps matching values)
      match_dates <- merge(profiles_EFF, profiles_INF, by = "Date")
      
      single_dates <- unique(match_dates$Date)
      
      for(date in colnames(daily_counts)){
        if(date %in% single_dates){
          daily_counts[row_num,date] <- daily_counts[row_num,date] + 1
          
        }
        
      }
      print(i)
    }
    
  }
  
  return(daily_counts)
}

##### Plotting function ####
daily_barplot_function <- function(counts, yaxis){
  # counts <- clust1_RT_filter[clust1_RT_filter$RT_match == TRUE,]
  # yaxis <- 5000

  counts$Reaction <- rownames(counts)
  m2 <- melt(counts)
  colnames(m2) <- c("Reaction","Date","Count")
  
  ggplot(data=m2, aes(x=Date, y=Count, fill = Reaction)) +
    geom_bar(stat="identity")
  
  m2$Reaction <- factor(m2$Reaction, levels=c("acetylation","hydroxylation","dihydroxylation","demethylation","deethylation","dehydrogenation","hydrogenation","dehydration","chlorine reduction",
                                                "deacetylation","glucuronidation","deglucuronidation","sulfonation","desulfonation","formylation","acylation","succinylation","fumarylation","malonylation"))

  
  m2$Date <- with(m2, ifelse(Date == "2020-11-18","Nov18",
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
  m2$Date <- factor(m2$Date, levels = c("Nov18","Nov19","Nov20","Nov21","Nov22",
                                                            "Nov23","Nov24","Nov25","Nov26","Nov27",
                                                            "Nov28","Nov29","Nov30","Dec1"))
  
  # figurename <- "reaction_counts_daily_230822.tiff"
  # 
  # tiff(figurename, units="in", width=10, height=8, res=300)
  gg <- ggplot(data=m2, aes(x=Date, y=Count, fill = Reaction))+
    geom_bar(stat="identity", size = 10)+
    # facet_wrap(~Place)+
    scale_fill_manual(values = c("#45D6FF","#FFC000","#C00000","#92D050","#FF83FA","#9575cd","#ff9797",
                                 "#ba632f","#2e7d32","#1c3ffd","#020873","#fa5b0f","#EE1289","#551A8B",
                                 "#99b3ff","#990099","#47d1d1","#e6e600","#c2c2a3"))+
    ylim(0,yaxis)+
    theme(text = element_text(size=18),
          # legend.position = "none",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.key.size = unit(0.7, 'cm'),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(color = "gray97"),
          panel.grid.minor = element_line(color = "gray97"),
          panel.border = element_rect(color = "grey", fill = NA))
  # dev.off()
  return(gg)
}

# function that calculates if each mz match meets RT criteria, and counts number of matches to use as filter later
RT_detect_criteria <- function(matches, workdat){

  # matches = clust1_ranged_matches
  # workdat = working_data
  
  # add empty RT matching row (to be filled with T/F)
  matches$RT_match <- rep(NA, nrow(matches))
  
  # try adding filter for detection ratio - but there are more than three detections per day?
  matches$detect_ratio <- rep(NA, nrow(matches))
  matches$matched_days <- rep(NA, nrow(matches))
  matches$INF_detects <- rep(NA, nrow(matches))
  matches$EFF_detects <- rep(NA, nrow(matches))
  
  for(i in 1:nrow(matches)){
    # i = 2
    INF_ID <- matches$prof_ID_INF[i]
    
    # take mean RT of all peaks detected with specified INF ID and change to minutes
    INF_RT <- mean(workdat$mean_RT[workdat$profile_ID == INF_ID])/60
    INF_detects <- workdat[workdat$profile_ID == INF_ID & workdat$Place == "INF",]
    INF_detects <- INF_detects %>% distinct(Date, .keep_all = TRUE) # only counts number of unique days detected (does not count triplicates)
    matches$INF_detects[i] <- nrow(INF_detects)

    EFF_ID <- matches$EFF_ID[i]
    EFF_RT <- mean(workdat$mean_RT[workdat$profile_ID == EFF_ID])/60
    EFF_detects <- workdat[workdat$profile_ID == EFF_ID & workdat$Place == "EFF",]
    EFF_detects <- EFF_detects %>% distinct(Date, .keep_all = TRUE)
    matches$EFF_detects[i] <- nrow(EFF_detects)
    
    # fill RT matching criteria for each INF/EFF pair
    if(EFF_RT < INF_RT + 2){
      matches$RT_match[i] <- TRUE
      
    }else{matches$RT_match[i] <- FALSE}
    
    matched_dates <- merge(EFF_detects, INF_detects, by = "Date")
    matches$matched_days[i] <- nrow(matched_dates)
    
    matched <- matches$matched_days[i]
    unmatch_INF <- matches$INF_detects[i]-matches$matched_days[i]
    unmatch_EFF <- matches$EFF_detects[i]-matches$matched_days[i]
    
    matches$detect_ratio[i] <- (2*matched)/(2*matched+(unmatch_INF)+(unmatch_EFF))
  }
  
  return(matches)
}

daily_counts_function_v2 <- function(matches, sampleinfo, bts){
  
  # matches <- test_set
  # sampleinfo <- base_data
  # bts <-btrules
  
  # here the variable "sampleinfo" is the "working_data" datframe with different dates that profiles are detected
  
  # Count the number of matches per day
  Date <-c("2020-11-18","2020-11-19","2020-11-20","2020-11-21",
           "2020-11-22","2020-11-23","2020-11-24","2020-11-25","2020-11-26","2020-11-27",
           "2020-11-28","2020-11-29","2020-11-30","2020-12-01")
  
  # dataframe to store the results of match counting per day
  daily_counts <- data.frame(matrix(nrow = 19, ncol = 14))
  rownames(daily_counts) <- bts$biotrans
  colnames(daily_counts) <- Date
  daily_counts[is.na(daily_counts)] <- 0
  
  # pairs up columns in total to rows in "daily_counts"
  pairs <- cbind(c(1:nrow(daily_counts)),c(3:21))
  
  
  for(h in 1:nrow(daily_counts)){
    # h = 1
    df.counts <- matches
    info <- sampleinfo
    col_num = pairs[h,][2] # this is the column number corresponding to the reaction type in the "total" data frame
    row_num = pairs[h,][1] # this is the row number corresponding to the reaction type in the "daily counts" results accounting
    
    
    # subsets whole dataframe to only include rows where a column entry is TRUE in the specified reaction column
    dat <- df.counts %>% 
      filter(df.counts[,col_num]!=0)
    
    if(nrow(dat) > 0){
      
      num_exp <- nrow(dat)
      
      print(num_exp)
      
      for(i in 1:nrow(dat)){
        # selected profile IDs associated with these TRUE hits from "dat" dataset
        profile_ID_EFF <- dat$EFF_ID[i]
        # profiles_EFF <- info[info$profile_ID %in% profile_ID_EFF,]
        # profiles_EFF <- profiles_EFF[profiles_EFF$Place == "EFF",]
        
        profile_ID_INF <- dat$prof_ID_INF[i]
        # profiles_INF <- info[info$profile_ID %in% profile_ID_INF,]
        # profiles_INF <- profiles_INF[profiles_INF$Place == "INF",]
        
        # only counts pair if INF and EFF profile ID are both present on the same day
        new_subset2 <- info[info$profile_ID == profile_ID_INF & info$Place == "INF" | info$profile_ID == profile_ID_EFF & info$Place == "EFF" ,]
        
        profiles_EFF <- new_subset2[new_subset2$Place == "EFF",]
        profiles_INF <- new_subset2[new_subset2$Place == "INF",]
        
        if(nrow(profiles_EFF) & nrow(profiles_INF) > 0){
          
          # only counts as a reaction when the parent and TP are detected on the same date (merge INF and EFF profiles is an inner join and only keeps matching values)
          match_dates <- merge(profiles_EFF, profiles_INF, by = "Date")
          
          single_dates <- unique(match_dates$Date)
          
          for(date in colnames(daily_counts)){
            if(date %in% single_dates){
              daily_counts[row_num,date] <- daily_counts[row_num,date] + 1
              
            }
          }
        }
        print(i)
      }
    }
  }
  
  return(daily_counts)
}

daily_barplot_function_left_top <- function(counts, yaxis){
  # counts <- clust1_RT_filter[clust1_RT_filter$RT_match == TRUE,]
  # yaxis <- 5000
  
  counts$Reaction <- rownames(counts)
  m2 <- melt(counts)
  colnames(m2) <- c("Reaction","Date","Count")
  
  ggplot(data=m2, aes(x=Date, y=Count, fill = Reaction)) +
    geom_bar(stat="identity")
  
  m2$Reaction <- factor(m2$Reaction, levels=c("acetylation","hydroxylation","dihydroxylation","demethylation","deethylation","dehydrogenation","hydrogenation","dehydration","chlorine reduction",
                                              "deacetylation","glucuronidation","deglucuronidation","sulfonation","desulfonation","formylation","acylation","succinylation","fumarylation","malonylation"))
  
  # figurename <- "reaction_counts_daily_230822.tiff"
  # 
  # tiff(figurename, units="in", width=10, height=8, res=300)
  gg <- ggplot(data=m2, aes(x=Date, y=Count, fill = Reaction))+
    geom_bar(stat="identity", size = 10)+
    # facet_wrap(~Place)+
    scale_fill_manual(values = c("#45D6FF","#FFC000","#C00000","#92D050","#FF83FA","#9575cd","#ff9797",
                                 "#ba632f","#2e7d32","#1c3ffd","#020873","#fa5b0f","#EE1289","#551A8B",
                                 "#99b3ff","#990099","#47d1d1","#e6e600","#c2c2a3"))+
    ylim(0,yaxis)+
    theme(text = element_text(size=18),
          # legend.text = element_text(size = 12),
          # legend.title = element_text(size = 14),
          # legend.key.size = unit(0.32, 'cm'),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(color = "gray97"),
          panel.grid.minor = element_line(color = "gray97"),
          panel.border = element_rect(color = "grey", fill = NA))
  # dev.off()
  return(gg)
}

daily_barplot_function_left <- function(counts, yaxis){
  # counts <- clust1_RT_filter[clust1_RT_filter$RT_match == TRUE,]
  # yaxis <- 5000
  
  counts$Reaction <- rownames(counts)
  m2 <- melt(counts)
  colnames(m2) <- c("Reaction","Date","Count")
  
  ggplot(data=m2, aes(x=Date, y=Count, fill = Reaction)) +
    geom_bar(stat="identity")
  
  m2$Reaction <- factor(m2$Reaction, levels=c("acetylation","hydroxylation","dihydroxylation","demethylation","deethylation","dehydrogenation","hydrogenation","dehydration","chlorine reduction",
                                              "deacetylation","glucuronidation","deglucuronidation","sulfonation","desulfonation","formylation","acylation","succinylation","fumarylation","malonylation"))
  
  # figurename <- "reaction_counts_daily_230822.tiff"
  # 
  # tiff(figurename, units="in", width=10, height=8, res=300)
  gg <- ggplot(data=m2, aes(x=Date, y=Count, fill = Reaction))+
    geom_bar(stat="identity", size = 10)+
    # facet_wrap(~Place)+
    scale_fill_manual(values = c("#45D6FF","#FFC000","#C00000","#92D050","#FF83FA","#9575cd","#ff9797",
                                 "#ba632f","#2e7d32","#1c3ffd","#020873","#fa5b0f","#EE1289","#551A8B",
                                 "#99b3ff","#990099","#47d1d1","#e6e600","#c2c2a3"))+
    ylim(0,yaxis)+
    theme(text = element_text(size=18),
          # legend.text = element_text(size = 12),
          # legend.title = element_text(size = 14),
          # legend.key.size = unit(0.32, 'cm'),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(color = "gray97"),
          panel.grid.minor = element_line(color = "gray97"),
          panel.border = element_rect(color = "grey", fill = NA))
  # dev.off()
  return(gg)
}

daily_barplot_function_right <- function(counts, yaxis){
  # counts <- clust1_RT_filter[clust1_RT_filter$RT_match == TRUE,]
  # yaxis <- 5000
  
  counts$Reaction <- rownames(counts)
  m2 <- melt(counts)
  colnames(m2) <- c("Reaction","Date","Count")
  
  ggplot(data=m2, aes(x=Date, y=Count, fill = Reaction)) +
    geom_bar(stat="identity")
  
  m2$Reaction <- factor(m2$Reaction, levels=c("acetylation","hydroxylation","dihydroxylation","demethylation","deethylation","dehydrogenation","hydrogenation","dehydration","chlorine reduction",
                                              "deacetylation","glucuronidation","deglucuronidation","sulfonation","desulfonation","formylation","acylation","succinylation","fumarylation","malonylation"))
  
  # figurename <- "reaction_counts_daily_230822.tiff"
  # 
  # tiff(figurename, units="in", width=10, height=8, res=300)
  gg <- ggplot(data=m2, aes(x=Date, y=Count, fill = Reaction))+
    geom_bar(stat="identity", size = 10)+
    # facet_wrap(~Place)+
    scale_fill_manual(values = c("#45D6FF","#FFC000","#C00000","#92D050","#FF83FA","#9575cd","#ff9797",
                                 "#ba632f","#2e7d32","#1c3ffd","#020873","#fa5b0f","#EE1289","#551A8B",
                                 "#99b3ff","#990099","#47d1d1","#e6e600","#c2c2a3"))+
    ylim(0,yaxis)+
    theme(text = element_text(size=18),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.key.size = unit(0.32, 'cm'),
          # legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(color = "gray97"),
          panel.grid.minor = element_line(color = "gray97"),
          panel.border = element_rect(color = "grey", fill = NA))
  # dev.off()
  return(gg)
}

daily_barplot_function_right_top <- function(counts, yaxis){
  # counts <- clust1_RT_filter[clust1_RT_filter$RT_match == TRUE,]
  # yaxis <- 5000
  
  counts$Reaction <- rownames(counts)
  m2 <- melt(counts)
  colnames(m2) <- c("Reaction","Date","Count")
  
  ggplot(data=m2, aes(x=Date, y=Count, fill = Reaction)) +
    geom_bar(stat="identity")
  
  m2$Reaction <- factor(m2$Reaction, levels=c("acetylation","hydroxylation","dihydroxylation","demethylation","deethylation","dehydrogenation","hydrogenation","dehydration","chlorine reduction",
                                              "deacetylation","glucuronidation","deglucuronidation","sulfonation","desulfonation","formylation","acylation","succinylation","fumarylation","malonylation"))
  
  # figurename <- "reaction_counts_daily_230822.tiff"
  # 
  # tiff(figurename, units="in", width=10, height=8, res=300)
  gg <- ggplot(data=m2, aes(x=Date, y=Count, fill = Reaction))+
    geom_bar(stat="identity", size = 10)+
    # facet_wrap(~Place)+
    scale_fill_manual(values = c("#45D6FF","#FFC000","#C00000","#92D050","#FF83FA","#9575cd","#ff9797",
                                 "#ba632f","#2e7d32","#1c3ffd","#020873","#fa5b0f","#EE1289","#551A8B",
                                 "#99b3ff","#990099","#47d1d1","#e6e600","#c2c2a3"))+
    ylim(0,yaxis)+
    theme(text = element_text(size=18),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.key.size = unit(0.32, 'cm'),
          # legend.position = "none",
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(color = "gray97"),
          panel.grid.minor = element_line(color = "gray97"),
          panel.border = element_rect(color = "grey", fill = NA))
  # dev.off()
  return(gg)
}


boxplot_function <- function(daily_counts_df,yaxis){
  # get boxplots
  daily_long <-melt(cbind(daily_counts_df, rownames(daily_counts_df)))
  colnames(daily_long) <-c("Reaction","Date","Count")
  
  ggplot(daily_long, aes(x = Reaction, y = Count))+
    geom_boxplot()
  
  
  daily_long$Reaction <- factor(daily_long$Reaction, levels=c("acetylation","hydroxylation","dihydroxylation","demethylation","deethylation","dehydrogenation","hydrogenation","dehydration","chlorine reduction",
                                                                          "deacetylation","glucuronidation","deglucuronidation","sulfonation","desulfonation","formylation","acylation","succinylation","fumarylation","malonylation"))
  
  
  gg <- ggplot(daily_long, aes(x=reorder(Reaction, -Count), y=Count, fill = Reaction)) + 
    geom_boxplot()+
    scale_fill_manual(values = c("#45D6FF","#FFC000","#C00000","#92D050","#FF83FA","#9575cd","#ff9797",
                                 "#ba632f","#2e7d32","#1c3ffd","#020873","#fa5b0f","#EE1289","#551A8B",
                                 "#99b3ff","#990099","#47d1d1","#e6e600","#c2c2a3"))+
    ylim(0,yaxis)+
    theme(text = element_text(size=18),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.key.size = unit(0.9, 'cm'),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 14),
          # axis.text.x = element_blank(),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(color = "gray97"),
          panel.grid.minor = element_line(color = "gray97"),
          panel.border = element_rect(color = "grey", fill = NA))
  
  return(gg)
}

save_function_Rdata <- function(objectname, folder){
  save(objectname, file = file.path("~/Project 3/All_objectives_P3WWMetabo/",folder,"/Robjects_toload/",objectname,".RData"))
  load(file.path("~/Project 3/All_objectives_P3WWMetabo/",folder,"/Robjects_toload/summ_peakarea.RData"))
}

