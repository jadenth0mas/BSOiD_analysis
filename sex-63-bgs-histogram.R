# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
# 
# Behavior by Sex (Supplementary Figure 3)

# Load dependencies
library(tidyverse)
library(rstatix)

# Read in data
df <- read_csv("mouse_data.csv")
df_fatal <- df %>% mutate(fatal=(Score==7))

preictal <- df_fatal %>% filter(scorer<=3600) %>% mutate(time="preictal")
myoclonic_seizure <- df_fatal %>% filter(scorer>3600 & scorer<=(3600+((GST-MST)*30))) %>% mutate(time="myoclonic seizure")
gen_seizure <- df_fatal %>% filter(scorer>(3600+((GST-MST)*30))) %>% mutate(time="generalized seizure")

# Get percentage of behavior group usage for each mouse
get_individual_df <- function(df) {
  p <- df %>% group_by(id) %>% count(BSOiD_labels)
  p2 <- p %>% pivot_wider(names_from=BSOiD_labels, values_from=n)
  p2[is.na(p2)] <- 0
  
  nums_wide <- p2[,-1]/rowSums(p2[,-1])
  nums_wide <- mutate(nums_wide, id=p2$id, strain=str_split_i(id, "-", 1))
  nums_wide <- relocate(nums_wide, c("id", "strain"))
  
  return(nums_wide)
}

# Get usage of 
get_gender_df <- function(df) {
  p <- df %>% group_by(Sex) %>% count(BSOiD_labels)
  p2 <- p %>% pivot_wider(names_from=BSOiD_labels, values_from=n)
  p2[is.na(p2)] <- 0
  return(p2)
}

p <- get_individual_df(preictal) %>% mutate(time="Preictal")
m <- get_individual_df(myoclonic_seizure) %>% mutate(time="Myoclonic Seizure")
g <- get_individual_df(gen_seizure) %>% mutate(time="Generalized Seizure")

total <- bind_rows(p, m, g)
sex_df <- unique(dplyr::select(df_fatal, c(id, Sex)))
total <- left_join(total, sex_df)
total <- total %>% dplyr::select(!c("53", "54", "57", "59"))

p_sex <- total %>% filter(time=="Preictal")
m_sex <- total %>% filter(time=="Myoclonic Seizure")
g_sex <- total %>% filter(time=="Generalized Seizure")



get_summary_df <- function(df) {
  s_df <- data.frame()
  
  for (col in names(df)) {
    
    if (col != "time" &&  col != "Sex" && col!="strain" && col!="id") {  # Assuming the column specifying groups is named 'mortality'
      
      male_mean <- mean(df[df$Sex=="M", col], na.rm=T)
      male_std_error <- sd(df[df$Sex == "M", col], na.rm = TRUE) / sqrt(sum(df$Sex == "M", na.rm = T))
      
      female_mean <- mean(df[df$Sex=="F", col], na.rm=T)
      female_std_error <- sd(df[df$Sex == "F", col], na.rm = TRUE) / sqrt(sum(df$Sex == "F", na.rm = T))

      ab_m_f_delta <- male_mean-female_mean
      
      temp <- data.frame(
        
        column = col,
        group = c("Male", "Female"),
        mean = c(male_mean, female_mean),
        std_error = c(male_std_error, female_std_error),
        mf_delta <- c(ab_m_f_delta)
      )
      s_df <- rbind(s_df, temp)
    }
  }
  s_df <- s_df %>% rename(mf_delta=`mf_delta....c.ab_m_f_delta.`)
  return(s_df)
}

p_summary <- get_summary_df(p_sex)
m_summary <- get_summary_df(m_sex)
g_summary <- get_summary_df(g_sex)

#write_csv(p_summary, "p_sex_summary.csv")
#write_csv(m_summary, "ms_sex_summary.csv")
#write_csv(g_summary, "gs_sex_summary.csv")

# Wilcox test for BG~Sex
adj_level <- 1 - (0.05)/63
model_wilcox <- function(var, arg1) {
  model <- wilcox_test(data=arg1, detailed=T, conf.level=adj_level, formula=as.formula(paste0(var,"~Sex")))
  return(model)
}

# Kruskal teest for BG~Sex
# model_kruskal <- function(var, arg1) {
#   model <- kruskal_test(data=arg1, formula=as.formula(paste0(var, "~Sex")))
#   return(model)
# }

# Get columns names for behavior group
colnames(p_sex) <- c("id", "strain", paste0("bg_", names(p_sex[3:65])), "time", "Sex")
wil_names <- names(p_sex)[3:65]

# Get df for wilcox test by sex for preictal usage
sex_wilcox <- lapply(wil_names, model_wilcox, arg1=p_sex)
psex_df <- bind_rows(sex_wilcox)
p_values <- psex_df$p
p_adjusted <- p.adjust(p_values, method="bonferroni")
psex_df$p_adjusted <- p_adjusted

# Get df for wilcox test by sex for MS usage
colnames(m_sex) <- c("id", "strain", paste0("bg_", names(p_sex[3:65])), "time", "Sex")
wil_names <- names(m_sex)[3:65]

sex_wilcox <- lapply(wil_names, model_wilcox, arg1=m_sex)
msex_df <- bind_rows(sex_wilcox)
p_values <- msex_df$p
p_adjusted <- p.adjust(p_values, method="bonferroni")
msex_df$p_adjusted <- p_adjusted

# Get df for wilcox test by sex for GS usage
colnames(g_sex) <- c("id", "strain", paste0("bg_", names(p_sex[3:65])), "time", "Sex")
wil_names <- names(g_sex)[3:65]

sex_wilcox <- lapply(wil_names, model_wilcox, arg1=g_sex)
gsex_df <- bind_rows(sex_wilcox)
p_values <- gsex_df$p
p_adjusted <- p.adjust(p_values, method="bonferroni")
gsex_df$p_adjusted <- p_adjusted


# Kruskal

# sex_kruskal <- lapply(wil_names, model_kruskal, arg1=p_sex)
# psex_df <- bind_rows(sex_kruskal)
# p_values <- psex_df$p
# p_adjusted <- p.adjust(p_values, method="bonferroni")
# psex_df$p_adjusted <- p_adjusted
# 
# colnames(m_sex) <- c("id", "strain", paste0("bg_", names(p_sex[3:65])), "time", "Sex")
# wil_names <- names(m_sex)[3:65]
# 
# sex_kruskal <- lapply(wil_names, model_kruskal, arg1=m_sex)
# msex_df <- bind_rows(sex_kruskal)
# p_values <- msex_df$p
# p_adjusted <- p.adjust(p_values, method="bonferroni")
# msex_df$p_adjusted <- p_adjusted
# 
# colnames(g_sex) <- c("id", "strain", paste0("bg_", names(p_sex[3:65])), "time", "Sex")
# wil_names <- names(g_sex)[3:65]
# 
# sex_kruskal <- lapply(wil_names, model_kruskal, arg1=g_sex)
# gsex_df <- bind_rows(sex_kruskal)
# p_values <- gsex_df$p
# p_adjusted <- p.adjust(p_values, method="bonferroni")
# gsex_df$p_adjusted <- p_adjusted
# 
# write_csv(psex_df, "kruskal_p_63bg_sex_wilcox.csv")
# write_csv(msex_df, "kruskal_m_63bg_sex_wilcox.csv")
# write_csv(gsex_df, "kruskal_g_63bg_sex_wilcox.csv")

# Remove bg_ prefix for P BGs
psex_df$bg <- gsub("^bg_", "", psex_df$.y.)
# Add significance for P
p_wilcox_df <- dplyr::select(psex_df, c(p_adjusted, p, bg)) %>% mutate(significance = case_when(
  p_adjusted < 0.001 ~ "***",
  p_adjusted < 0.01 ~ "**",
  p_adjusted < 0.05 ~ "*", 
  TRUE ~ ""
))

# Remove bg_ prefix for MS BGs
msex_df$bg <- gsub("^bg_", "", msex_df$.y.)
# Add significance for P
m_wilcox_df <- dplyr::select(msex_df, c(p_adjusted, p, bg)) %>% mutate(significance = case_when(
  p_adjusted < 0.001 ~ "***",
  p_adjusted < 0.01 ~ "**",
  p_adjusted < 0.05 ~ "*", 
  TRUE ~ ""
))

# Remove bg_ prefix for GS BGs
gsex_df$bg <- gsub("^bg_", "", gsex_df$.y.)
# Add significance for G
g_wilcox_df <- dplyr::select(gsex_df, c(p_adjusted, p, bg)) %>% mutate(significance = case_when(
  p_adjusted < 0.001 ~ "***",
  p_adjusted < 0.01 ~ "**",
  p_adjusted < 0.05 ~ "*", 
  TRUE ~ ""
))

# Get df in format for plotting with labels
p_summary <- p_summary %>% left_join(p_wilcox_df, by=join_by(column==bg))
m_summary <- m_summary %>% left_join(m_wilcox_df, by=join_by(column==bg))
g_summary <- g_summary %>% left_join(g_wilcox_df, by=join_by(column==bg))

# Supplementary Figure 3A
ggplot(p_summary, aes(x=reorder(column, -mf_delta), y=mean, color=group, ymin=mean-std_error, ymax=mean+std_error)) +
  geom_point(size=3) +
  geom_errorbar(width=0.2) +
  labs(x="Behavior Group", y="Percent Usage", title="Preictal Mean and Standard Error Plot for Behavior Group by Sex", color="Sex") +
  theme_minimal() +
  geom_text(data=p_summary %>% filter(group=="Male"),
           aes(x=reorder(column, -mf_delta), label=significance, y=.15), color="black", fontface="bold", size=6, angle=90, vjust=1, hjust="right") +
  theme(axis.text.x=element_text(angle=90, size=12, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=12, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14))

# Supplementary Figure 3B
ggplot(m_summary, aes(x=reorder(column, -mf_delta), y=mean, color=group, ymin=mean-std_error, ymax=mean+std_error)) +
  geom_point(size=3) +
  geom_errorbar(width=0.2) +
  labs(x="Behavior Group", y="Percent Usage", title="Myoclonic Seizure Mean and Standard Error Plot for Behavior Group by Sex", color="Sex") +
  theme_minimal() +
  geom_text(data=m_summary %>% filter(group=="Male"),
            aes(x=reorder(column, -mf_delta), label=significance, y=.15), color="black", fontface="bold", size=6, angle=90, vjust=.8, hjust="right") +
  theme(axis.text.x=element_text(angle=90, size=12, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=12, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14))

# Supplementary Figure 3C
ggplot(g_summary, aes(x=reorder(column, -mf_delta), y=mean, color=group, ymin=mean-std_error, ymax=mean+std_error)) +
  geom_point(size=3) +
  geom_errorbar(width=0.2) +
  labs(x="Behavior Group", y="Percent Usage", title="Generalized Seizure Mean and Standard Error Plot for Behavior Group by Sex", color="Sex") +
  theme_minimal() +
  geom_text(data=g_summary %>% filter(group=="Male"),
            aes(x=reorder(column, -mf_delta), label=significance, y=.15), color="black", fontface="bold", size=6, angle=90, vjust=1, hjust="right") +
  theme(axis.text.x=element_text(angle=90, size=12, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=12, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14))
