# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
# 
# Behavior biomarker of SUDEP 

# Load dependencies
library(tidyverse)
library(ggsci)
library(ggrepel)

# Read data
df <- read_csv("mouse_data.csv")
df_fatal <- df %>% mutate(fatal=(Score==7))

# Partition data by seizure state
preictal <- df_fatal %>% filter(scorer<=3600) %>% mutate(time="preictal")
myoclonic_seizure <- df_fatal %>% filter(scorer>3600 & scorer<=(3600+((GST-MST)*30))) %>% mutate(time="myoclonic seizure")
gen_seizure <- df_fatal %>% filter(scorer>(3600+((GST-MST)*30))) %>% mutate(time="generalized seizure")

get_individual_df <- function(df) {
  p <- df %>% group_by(id) %>% count(BSOiD_labels)
  p2 <- p %>% pivot_wider(names_from=BSOiD_labels, values_from=n)
  p2[is.na(p2)] <- 0
  
  nums_wide <- p2[,-1]/rowSums(p2[,-1])
  nums_wide <- mutate(nums_wide, id=p2$id, strain=str_split_i(id, "-", 1))
  nums_wide <- relocate(nums_wide, c("id", "strain"))
  
  return(nums_wide)
}

p <- get_individual_df(preictal) %>% mutate(time="Preictal")
m <- get_individual_df(myoclonic_seizure) %>% mutate(time="Myoclonic Seizure")
g <- get_individual_df(gen_seizure) %>% mutate(time="Generalized Seizure")

total <- bind_rows(p, m, g) %>% left_join(unique(dplyr::select(df_fatal, c(id, fatal))))
# Remove BGs
total <- total %>% dplyr::select(!c("53", "54", "57", "59"))

# Get mean and sd by fatality
get_summary_df <- function(df) {
  # Initialize summary df
  s_df <- data.frame()
  
  for (col in names(df)) {
    # For BGs
    if (col != "time" && col != "fatal" && col!="id" && col!="strain") {  # Assuming the column specifying groups is named 'mortality'
      # Get mean and sd of BG col for fatal
      f_mean <- mean(df[df$fatal == T, col], na.rm = TRUE)
      f_std_error <- sd(df[df$fatal == T, col], na.rm = TRUE) / sqrt(sum(df$fatal == T, na.rm = T))
      # Get mean and sd for BG col for nonfatal
      nf_mean <- mean(df[df$fatal == F, col], na.rm=T)
      nf_std_error <- sd(df[df$fatal == F, col], na.rm = TRUE) / sqrt(sum(df$fatal == F, na.rm = T))
      
      # Get delta between fatal and nonfatal
      # Get actual difference between f and nf
      delta <- abs(f_mean - nf_mean)/2
      f_minus_nf <- f_mean-nf_mean

      temp <- data.frame(
        
        column = col,
        group = c("Fatal", "Non-Fatal"),
        mean = c(f_mean, nf_mean),
        std_error = c(f_std_error, nf_std_error),
        delta_total = c(delta),
        f_m_nf = c(f_minus_nf)
      )
      s_df <- rbind(s_df, temp)
    }
  }
  return(s_df)
}

# Get summary stats for P, MS, GS by fatality
p_summary <- total %>% filter(time=="Preictal") %>% get_summary_df()
m_summary <- total %>% filter(time=="Myoclonic Seizure") %>% get_summary_df()
g_summary <- total %>% filter(time=="Generalized Seizure") %>% get_summary_df()

gs_freq_data <- g_summary %>% dplyr::select(c(column, group, mean)) %>% pivot_wider(names_from=column, values_from=mean)
#write_csv(gs_freq_data, 'gs_fatal_frequency_for_ethnogram.csv')

# Kruskal test by fatality
# model_kruskal <- function(var, arg1) {
#   model <- kruskal_test(data=arg1, formula=as.formula(paste0(var, "~Fatal")))
#   return(model)
# }

# Wilcox test by fatality
conf_adj <- 1 - (0.05)/63

model_wilcox <- function(var, arg1) {
  model <- wilcox_test(data=arg1, detailed=T, conf.level=conf_adj, formula=as.formula(paste0(var,"~Fatal")))
  return(model)
}


# Get df bg usage with fatal column for P, MS, GS
p_fatal <- total %>% filter(time=="Preictal") %>% mutate(fatal=as.factor(fatal))
m_fatal <- total %>% filter(time=="Myoclonic Seizure") %>% mutate(fatal=as.factor(fatal))
g_fatal <- total %>% filter(time=="Generalized Seizure") %>% mutate(fatal=as.factor(fatal))

# Get cols in format for statistical tests
colnames(p_fatal) <- c("id", "strain", paste0("bg_", names(p_fatal[3:65])), "time", "Fatal")
colnames(m_fatal) <- c("id", "strain", paste0("bg_", names(m_fatal[3:65])), "time", "Fatal")
colnames(g_fatal) <- c("id", "strain", paste0("bg_", names(g_fatal[3:65])), "time", "Fatal")

# Get BGs names
wil_names <- names(p_fatal)[3:65]

# Mann WHitney U for P, MS, and GS
pf_wilcox <- lapply(wil_names, model_wilcox, arg1=p_fatal)
pf_pvals <- sapply(pf_wilcox, function(x) x$p)
pf_adj_pvals <- p.adjust(pf_pvals, method="bonferroni")
pf_wilcox <- bind_rows(pf_wilcox)
 
mf_wilcox <- lapply(wil_names, model_wilcox, arg1=m_fatal)
mf_pvals <- sapply(mf_wilcox, function(x) x$p)
mf_adj_pvals <- p.adjust(mf_pvals, method="bonferroni")
mf_wilcox <- bind_rows(mf_wilcox)
 
gf_wilcox <- lapply(wil_names, model_wilcox, arg1=g_fatal)
gf_pvals <- sapply(gf_wilcox, function(x) x$p)
gf_adj_pvals <- p.adjust(gf_pvals, method="bonferroni")
gf_wilcox <- bind_rows(gf_wilcox)


# pf_kruskal <- lapply(wil_names, model_kruskal, arg1=p_fatal)
# pf_pvals <- sapply(pf_kruskal, function(x) x$p)
# pf_adj_pvals <- p.adjust(pf_pvals, method="bonferroni")
#  
# mf_kruskal <- lapply(wil_names, model_kruskal, arg1=m_fatal)
# mf_pvals <- sapply(mf_kruskal, function(x) x$p)
# mf_adj_pvals <- p.adjust(mf_pvals, method="bonferroni")
#  
# gf_kruskal <- lapply(wil_names, model_kruskal, arg1=g_fatal)
# gf_pvals <- sapply(gf_kruskal, function(x) x$p)
# gf_adj_pvals <- p.adjust(gf_pvals, method="bonferroni")

# Get df of stat results


p_fatal_stat <- data.frame(behavior_group=wil_names, p_value=pf_pvals, p_adjusted=pf_adj_pvals, estimate = pf_wilcox$estimate, conf.low=pf_wilcox$conf.low, conf.high=pf_wilcox$conf.high)
m_fatal_stat <- data.frame(behavior_group=wil_names, p_value=mf_pvals, p_adjusted=mf_adj_pvals, estimate=mf_wilcox$estimate, conf.low=mf_wilcox$conf.low, conf.high=mf_wilcox$conf.high)
g_fatal_stat <- data.frame(behavior_group=wil_names, p_value=gf_pvals, p_adjusted=gf_adj_pvals, estimate = gf_wilcox$estimate, conf.low=gf_wilcox$conf.low, conf.high=gf_wilcox$conf.high)

# write_csv(p_fatal_stat, "kruskal_p_fatal_63bgs_stats.csv")
# write_csv(m_fatal_stat, "kruskal_m_fatal_63bgs_stats.csv")
# write_csv(g_fatal_stat, "kruskal_g_fatal_63bgs_stats.csv")

# Get labels for significance for P, MS, GS
p_fatal_stat$bg <- gsub("^bg_", "", p_fatal_stat$behavior_group)
p_fatal_stat <- dplyr::select(p_fatal_stat, c(p_adjusted, p_value, bg)) %>% mutate(significance = case_when(
  p_adjusted < 0.001 ~ "***",
  p_adjusted < 0.01 ~ "**",
  p_adjusted < 0.05 ~ "*", 
  TRUE ~ ""
))

m_fatal_stat$bg <- gsub("^bg_", "", m_fatal_stat$behavior_group)
m_fatal_stat <- dplyr::select(m_fatal_stat, c(p_adjusted, p_value, bg)) %>% mutate(significance = case_when(
  p_adjusted < 0.001 ~ "***",
  p_adjusted < 0.01 ~ "**",
  p_adjusted < 0.05 ~ "*", 
  TRUE ~ ""
))


g_fatal_stat$bg <- gsub("^bg_", "", g_fatal_stat$behavior_group)
g_fatal_stat <- dplyr::select(g_fatal_stat, c(p_adjusted, p_value, bg)) %>% mutate(significance = case_when(
  p_adjusted < 0.001 ~ "***",
  p_adjusted < 0.01 ~ "**",
  p_adjusted < 0.05 ~ "*", 
  TRUE ~ ""
))

# Get df in format to plot
p_summary <- p_summary %>% left_join(p_fatal_stat, by=join_by(column==bg))
m_summary <- m_summary %>% left_join(m_fatal_stat, by=join_by(column==bg))
g_summary <- g_summary %>% left_join(g_fatal_stat, by=join_by(column==bg))

ggplot(p_summary, aes(x=reorder(column, -f_m_nf), y=mean, color=group, ymin=mean-std_error, ymax=mean+std_error)) +
  geom_point(size=3) +
  geom_errorbar(width=0.2) +
  scale_color_nejm() +
  labs(x="Behavior Group", y="Percent Usage", title="Mean and Standard Error Plot for Preictal Fatal vs Non-Fatal") +
  theme_minimal() +
  geom_text(data=p_summary %>% filter(group=="Male"),
            aes(x=reorder(column, -f_m_nf), label=significance, y=.15), color="black", fontface="bold", size=6, angle=90, vjust=1, hjust="right") +
  theme(axis.text.x=element_text(angle=0, size=16, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=16, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=18),
        axis.ticks.length=unit(2, "mm"),
        legend.text=element_text(angle=0, size=15), legend.position="none") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

ggplot(m_summary, aes(x=reorder(column, -f_m_nf), y=mean, color=group, ymin=mean-std_error, ymax=mean+std_error)) +
  geom_point(size=3) +
  geom_errorbar(width=0.2) +
  scale_color_nejm() +
  labs(x="Behavior Group", y="Percent Usage", title="Mean and Standard Error Plot for MS Fatal vs Non-Fatal") +
  theme_minimal() +
  geom_text(data=m_summary %>% filter(group=="Male"),
            aes(x=reorder(column, -f_m_nf), label=significance, y=.15), color="black", fontface="bold", size=6, angle=90, vjust=1, hjust="right") +
  theme(axis.text.x=element_text(angle=0, size=26, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=16, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=18),
        axis.ticks.length=unit(2, "mm"),
        legend.text=element_text(angle=0, size=15), legend.position="none") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

# Figure 4A
ggplot(g_summary, aes(x=reorder(column, -f_m_nf), y=mean, color=group, ymin=mean-std_error, ymax=mean+std_error)) +
  geom_point(size=3) +
  geom_errorbar(width=0.2) +
  scale_color_nejm() +
  labs(x="Behavior Group", y="Percent Usage", color="Group", title="GS Fatal vs Non-Fatal") +
  theme_minimal() +
  geom_text(data=g_summary %>% filter(group=="Male"),
            aes(x=reorder(column, -f_m_nf), label=significance, y=.15), color="black", fontface="bold", size=6, angle=90, vjust=1, hjust="right") +
  theme(axis.text.x=element_text(angle=0, size=14, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=14, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.ticks.length=unit(2, "mm"),
        legend.text=element_text(angle=0, size=14)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
