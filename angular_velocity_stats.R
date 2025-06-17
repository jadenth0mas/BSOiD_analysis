# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
#
# Behavior biomarker of SUDEP 

# Load dependencies
library(tidyverse)
library(ggthemes)
library(ggsci)
library(rstatix)

# Load angular velocity data
fatal_angular <- read_csv("angular_velocity_fatal.csv")
nonfatal_angular <- read_csv("angular_velocity_nonfatal.csv")

# Split by extension and flexion within fatal and nonfatal
fatal_extension <- fatal_angular %>% filter(`Angular Velocity`>0) %>% mutate(group="fatal")
nonfatal_extension <- nonfatal_angular %>% filter(`Angular Velocity`>0) %>% mutate(group="nonfatal")

fatal_flexion <- fatal_angular %>% filter(`Angular Velocity`<0) %>% mutate(group="fatal")
nonfatal_flexion <- nonfatal_angular %>% filter(`Angular Velocity`<0) %>% mutate(group="nonfatal")

# Get full df of extension and flexion
extension <- bind_rows(fatal_extension, nonfatal_extension)
flexion <- bind_rows(fatal_flexion, nonfatal_flexion)

# Make fatal or nonfatal factor
extension$group = as.factor(extension$group)
flexion$group = as.factor(flexion$group)

# KS test by fatal or nonfatal for extension and flexion
extension_ks <- ks.test(x=extension[extension$group=="fatal",]$`Angular Velocity`,
                        y=extension[extension$group=='nonfatal',]$`Angular Velocity`)


flexion_ks <- ks.test(x=flexion[flexion$group=="fatal",]$`Angular Velocity`,
                              y=flexion[flexion$group=='nonfatal',]$`Angular Velocity`)

# DF of both results
s <- data.frame(var=c("Extension", "Flexion"), method=rep(extension_ks$method, 2),
           statistic=c(extension_ks$statistic, flexion_ks$statistic),
           p.val=c(extension_ks$p.value, flexion_ks$p.value))

#write_csv(s, "flexion_extension_ks.csv")