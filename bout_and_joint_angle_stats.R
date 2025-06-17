# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
#
# Behavior biomarker of SUDEP 

# Load Dependencies
library(tidyverse)
library(ggthemes)
library(ggsci)
library(rstatix)

# Get Bout data
successive <- read_csv("cleaned_successive_frames.csv")

# Split by fatality
fatal_filter <- successive$`Experimental Group`=="fatal"
f_successive <- successive[fatal_filter,]
nf_successive <- successive[!fatal_filter,]

# Get BGs to analyze bout
labels <- unique(f_successive$`B-SOiD Label`)

model_ks <- function(var) {
  model <- ks.test(x=f_successive[f_successive$`B-SOiD Label`==var,]$`Successive Frame Count`,
                   y=nf_successive[nf_successive$`B-SOiD Label`==var,]$`Successive Frame Count`)
  return(model)
}

z <- lapply(labels, model_ks)
p_values <- sapply(z, function(x) x$p)
statistics <- sapply(z, function(x) x$statistic)
p_adjusted <- p.adjust(p_values, method="bonferroni")
stats_df <- data.frame(bsoid=labels, statistics=statistics, p_values=p_values, p_adjusted=p_adjusted)

#write_csv(stats_df, "bout_stats.csv")


# Get joint angles for fatal and nonfatal for each seizure state
# Do KS test for distributions

gs_fatal_joints <- read_csv("final_GS_fatal_joint_angles.csv")
gs_nonfatal_joints <- read_csv("final_GS_nonfatal_joint_angles.csv")
gs_joints_stats <- ks.test(x=gs_fatal_joints$`Joint Angle`, y=gs_nonfatal_joints$`Joint Angle`)

ms_fatal_joints <- read_csv("final_MS_fatal_joint_angles.csv")
ms_nonfatal_joints <- read_csv("final_MS_nonfatal_joint_angles.csv")
ms_joints_stats <- ks.test(x=ms_fatal_joints$`Joint Angle`, y=ms_nonfatal_joints$`Joint Angle`)


p_fatal_joints <- read_csv("final_P_fatal_joint_angles.csv")
p_nonfatal_joints <- read_csv("final_P_nonfatal_joint_angles.csv")
p_joints_stats <- ks.test(x=p_fatal_joints$`Joint Angle`, y=p_nonfatal_joints$`Joint Angle`)

p_vals <- c(p_joints_stats$p.value, ms_joints_stats$p.value, gs_joints_stats$p.value)
stats <- c(p_joints_stats$statistic, ms_joints_stats$statistic, gs_joints_stats$statistic)
states <- c("P", "MS", "GS")
joints_df <- data.frame(time=states, statistic=stats, p.value=p_vals)

#write_csv(joints_df, "joints_stats.csv")
