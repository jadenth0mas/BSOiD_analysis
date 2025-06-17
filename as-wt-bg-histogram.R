# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
#
# Fine Seizure Behavior Phenotyping of Angelman Syndrome Model Mice

# Load dependencies
library(tidyverse)
library(ggsci)
library(rstatix)

# Load data
df_fatal <- read_csv("mouse_data4.csv") %>% mutate(fatal=(Score==7))

# Split by seizure state
preictal <- df_fatal %>% filter(scorer<=3600) %>% mutate(time="preictal")
myoclonic_seizure <- df_fatal %>% filter(scorer>3600 & scorer<=(3600+((GST-MST)*30))) %>% mutate(time="myoclonic seizure")
gen_seizure <- df_fatal %>% filter(scorer>(3600+((GST-MST)*30))) %>% mutate(time="generalized seizure")

# Get percent usage of each bg for each individual
get_individual_df <- function(df) {
  p <- df %>% group_by(id, Genotype) %>% count(BSOiD_labels)
  p2 <- p %>% pivot_wider(names_from=BSOiD_labels, values_from=n, id_cols=c(id, Genotype))
  p2[is.na(p2)] <- 0
  
  nums_wide <- p2[,-c(1,2)]/rowSums(p2[,-c(1,2)])
  nums_wide <- mutate(nums_wide, id=p2$id, genotype=p2$Genotype)
  nums_wide <- relocate(nums_wide, c("id", "genotype"))
  
  return(nums_wide)
}

# Get BG usage for each time period
p <- get_individual_df(preictal) %>% mutate(time="Preictal") %>% dplyr::select(!c("53", "54", "57", "59"))
m <- get_individual_df(myoclonic_seizure) %>% mutate(time="Myoclonic Seizure") %>% dplyr::select(!c("53", "54", "57", "59"))
g <- get_individual_df(gen_seizure) %>% mutate(time="Generalized Seizure") %>% dplyr::select(!c("53", "54", "57", "59"))

total <- bind_rows(p, m, g)
total[is.na(total)] <- 0

# Get summary stat by genotype
get_summary_df <- function(df) {
  s_df <- data.frame()
  
  for (col in names(df)) {
    if (col != "time" && col != "genotype" && col!="id") {
      wt_mean <- mean(df[df$genotype == "WT", col], na.rm = TRUE)
      wt_std_error <- sd(df[df$genotype == "WT", col], na.rm = TRUE) / sqrt(sum(df$genotype == "WT", na.rm = T))
      
      as_mean <- mean(df[df$genotype =="AS", col], na.rm=T)
      as_std_error <- sd(df[df$genotype == "AS", col], na.rm = TRUE) / sqrt(sum(df$genotype == "AS", na.rm = T))
      
      delta <- abs(wt_mean - as_mean)/2
      as_minus_wt <- as_mean-wt_mean
      
      temp <- data.frame(
        
        column = col,
        group = c("Wild Type", "Angelman Syndrome"),
        mean = c(wt_mean, as_mean),
        std_error = c(wt_std_error, as_std_error),
        delta_total = c(delta),
        as_m_wt = c(as_minus_wt)
      )
      s_df <- rbind(s_df, temp)
    }
  }
  return(s_df)
}

# Get summary stats for each seizure state
p_summary <- get_summary_df(p)
ms_summary <- get_summary_df(m)
gs_summary <- get_summary_df(g)

# Make genotype an ordered factor
p_summary$group <- factor(p_summary$group, levels=c("Wild Type", "Angelman Syndrome"))
ms_summary$group <- factor(ms_summary$group, levels=c("Wild Type", "Angelman Syndrome"))
gs_summary$group <- factor(gs_summary$group, levels=c("Wild Type", "Angelman Syndrome"))

# Get fold Change between each
# ms_fold <- ms_summary %>% dplyr::select(c(column, mean, group)) %>%
#   pivot_wider(id_cols=column, names_from=group, values_from=mean) %>%
#   mutate(
#     as_wt_fold_percentage=(`Angelman Syndrome`-`Wild Type`)/`Wild Type`*100
#   )
# 
# gs_fold <- gs_summary %>% dplyr::select(c(column, mean, group)) %>%
#   pivot_wider(id_cols=column, names_from=group, values_from=mean) %>%
#   mutate(
#     as_wt_fold_percentage=(`Angelman Syndrome`-`Wild Type`)/`Wild Type`*100
#   )
# 
# write_csv(ms_fold, "as-wt-ms-fold.csv")
# write_csv(gs_fold, "as-wt-gs-fold.csv")


# Get numeric data for each seizure state
p_d <- p %>% dplyr::select(!time)
m_d <- m %>% dplyr::select(!time)
g_d <- g %>% dplyr::select(!time)

vars <- names(p_d)[-c(1, 2)]
nc <- paste0("bg_", vars)
new_cols <- c("id", "genotype", nc)
colnames(p_d) <- new_cols


m_vars <- names(m_d)[-c(1,2)]
m_nc <- paste0("bg_", m_vars)
m_new_cols <- c("id", "genotype", m_nc)
colnames(m_d) <- m_new_cols

g_vars <- names(g_d)[-c(1,2)]
g_nc <- paste0("bg_", g_vars)
g_new_cols <- c("id", "genotype", g_nc)
colnames(g_d) <- g_new_cols

num_p <- dplyr::select(p_d, !c("id", "genotype"))
num_m <- dplyr::select(m_d, !c("id", "genotype"))
num_g <- dplyr::select(g_d, !c("id", "genotype"))


adj_level <- 1 - (0.05 / 63)

# Wilcox test for P
model_wilcox_p <- function(var) {
  model <- wilcox_test(data=p_d, formula=as.formula(paste0(var,"~genotype")), p.adjust.method="bonferroni", conf.level = adj_level, detailed=T)
  return(model)
}

# Wilcox test for M
model_wilcox_m <- function(var) {
  model <- wilcox_test(data=m_d, formula=as.formula(paste0(var,"~genotype")), p.adjust.method="bonferroni", conf.level=adj_level, detailed=T)
  return(model)
}

# Wilcox test for G
model_wilcox_g <- function(var) {
  model <- wilcox_test(data=g_d, formula=as.formula(paste0(var,"~genotype")), p.adjust.method="bonferroni", conf.level = adj_level, detailed=T)
  return(model)
}

# model_kruskal_p <- function(var) {
#   model <- kruskal_test(data=p_d, formula=as.formula(paste0(var,"~genotype")))
#   return(model)
# }
# 
# model_kruskal_m <- function(var) {
#   model <- kruskal_test(data=m_d, formula=as.formula(paste0(var,"~genotype")))
#   return(model)
# }
# 
# model_kruskal_g <- function(var) {
#   model <- kruskal_test(data=g_d, formula=as.formula(paste0(var,"~genotype")))
#   return(model)
# }


# Mann whitney U
p_stats <- bind_rows(lapply(nc, model_wilcox_p))
m_stats <- bind_rows(lapply(m_nc, model_wilcox_m))
g_stats <- bind_rows(lapply(g_nc, model_wilcox_g))

# p_stats <- bind_rows(lapply(nc, model_kruskal_p))
# m_stats <- bind_rows(lapply(m_nc, model_kruskal_m))
# g_stats <- bind_rows(lapply(g_nc, model_kruskal_g))

# Get adjusted p-value
p_stats$p_adjusted <- p.adjust(p_stats$p, method="bonferroni")
m_stats$p_adjusted <- p.adjust(m_stats$p, method="bonferroni")
g_stats$p_adjusted <- p.adjust(g_stats$p, method="bonferroni")

# write_csv(p_stats, "kruskal_d4_preictal_bg_stats_63bgs.csv")
# write_csv(m_stats, "kruskal_d4_ms_bg_stats_63bgs.csv")
# write_csv(g_stats, "kruskal_d4_gs_bg_stats_63bgs.csv")

#write_csv(p_stats, "d4_preictal_bg_stats_63bgs.csv")
#write_csv(m_stats, "d4_ms_bg_stats_63bgs.csv")
#write_csv(g_stats, "d4_gs_bg_stats_63bgs.csv")


# Preictal histogram plot
p_summary %>% ggplot(aes(x=reorder(column, -as_m_wt), y=mean, color=group, ymin=mean-std_error, ymax=mean+std_error)) +
  geom_point(size=3) +
  geom_errorbar(width=.2) +
  labs(title="Mean and Standard Error Plot for Preictal WT and AS",
       x="Behavior Group", y="Average Percent Usage", color="Genotype") +
  theme_minimal() +
  scale_color_bmj() +
  theme(axis.text.x=element_text(angle=0, size=17, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=17, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=18),
        axis.ticks.length=unit(2, "mm")) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

# Figure 7D
ms_summary %>% ggplot(aes(x=reorder(column, -as_m_wt), y=mean, color=group, ymin=mean-std_error, ymax=mean+std_error)) +
  geom_point(size=3) +
  geom_errorbar(width=.2) +
  labs(title="MS",
       x="Behavior Group", y="Average Percent Usage", color="Group") +
  theme_minimal() +
  scale_color_bmj() +
  theme(axis.text.x=element_text(angle=0, size=16, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=16, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=18),
        axis.ticks.length=unit(2, "mm"),
        legend.text=element_text(angle=0, size=15), legend.position="none") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

# Figure 7F
gs_summary %>% ggplot(aes(x=reorder(column, -as_m_wt), y=mean, color=group, ymin=mean-std_error, ymax=mean+std_error)) +
  geom_point(size=3) +
  geom_errorbar(width=.2) +
  labs(title="GS",
       x="Behavior Group", y="Average Percent Usage",
       color="Group") +
  theme_minimal() +
  scale_color_bmj() +
  theme(axis.text.x=element_text(angle=0, size=16, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=16, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=18),
        axis.ticks.length=unit(2, "mm"),
        legend.text=element_text(angle=0, size=15), legend.position="none") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))


