# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
#
# Fine behavior discrimination of seizure states section

# Load Packages
library(tidyverse)
library(rstatix)
library(PMCMRplus)
library(MASS)
library(ggrepel)


# Read in dataset
df <- read_csv("mouse_data.csv")
# Set a SUDEP flag for if the racine score is 7
df_fatal <- df %>% mutate(fatal=(Score==7))

# Split df into time periods
preictal <- df_fatal %>% filter(scorer<=3600) %>% mutate(time="preictal")
myoclonic_seizure <- df_fatal %>% filter(scorer>3600 & scorer<=(3600+((GST-MST)*30))) %>% mutate(time="myoclonic seizure")
gen_seizure <- df_fatal %>% filter(scorer>(3600+((GST-MST)*30))) %>% mutate(time="generalized seizure")

# Get the df with observations as rows and behavior groups as columns, with the 
# value as a percentage of the given time

get_individual_df <- function(df) {
  # Group by mouse id and count each label occurance
  p <- df %>% group_by(id) %>% count(BSOiD_labels) %>% ungroup()
  p2 <- p %>% pivot_wider(names_from=BSOiD_labels, values_from=n)
  p2[is.na(p2)] <- 0
  
  
  # Divide the values by the row sums, so each cell is the percentage of the
  # time period that the mouse displays BG #
  nums_wide <- p2[,-1]/rowSums(p2[,-1])
  nums_wide <- mutate(nums_wide, id=p2$id, strain=str_split_i(id, "-", 1))
  nums_wide <- relocate(nums_wide, c("id", "strain"))
  
  # Return percentage data with mouse id and strain
  return(nums_wide)
}

# Get for each of preictal, MS, and GS
p <- get_individual_df(preictal) %>% mutate(time="Preictal")
m <- get_individual_df(myoclonic_seizure) %>% mutate(time="Myoclonic Seizure")
g <- get_individual_df(gen_seizure) %>% mutate(time="Generalized Seizure")

# Create a total dataframe for analysis with time column for filtering
total <- bind_rows(p, m, g)


# Remove BGs
total <- total %>% dplyr::select(!c("53", "54", "57", "59"))

# Initialize a dataframe for the summary
summary_stats <- data.frame()

# Select the numeric and time columns
data <- total %>% dplyr::select(!c(strain, id))


# For all columns
for (col in names(data)) {
  # If the column is a BG
  if (col != "time") {
    # Get mean and standard error for each behavior group for each time period
    
    preictal_mean <- mean(data[data$time == "Preictal", col], na.rm = TRUE)
    preictal_std_error <- sd(data[data$time == "Preictal", col], na.rm = TRUE) / sqrt(sum(data$time == "Preictal"))
    ms_mean <- mean(data[data$time == "Myoclonic Seizure", col], na.rm = TRUE)
    ms_std_error <- sd(data[data$time == "Myoclonic Seizure", col], na.rm = TRUE) / sqrt(sum(data$time == "Myoclonic Seizure"))
    gs_mean <- mean(data[data$time == "Generalized Seizure", col], na.rm = TRUE)
    gs_std_error <- sd(data[data$time == "Generalized Seizure", col], na.rm = TRUE) / sqrt(sum(data$time == "Generalized Seizure"))
    
    # Get differences between time periods
    p_minus_p <- (0)
    m_minus_p <- (ms_mean-preictal_mean)
    g_minus_p <- (gs_mean-preictal_mean)
    
    # Get absolute differences between time periods
    ab_P_MS_delta <- abs(preictal_mean - ms_mean)
    ab_MS_GS_delta <- abs(ms_mean - gs_mean)
    ab_P_GS_delta <- abs(preictal_mean-gs_mean)
    
    # Get average absolute difference
    dt <- (ab_P_GS_delta + ab_MS_GS_delta + ab_P_MS_delta)/3
    
    # Create a temporary dataframe with the information from BG
    temp <- data.frame(
      # BG identifier
      column = col,
      group = c("Preictal", "MS", "GS"),
      mean = c(preictal_mean, ms_mean, gs_mean),
      std_error = c(preictal_std_error, ms_std_error, gs_std_error),
      delta_pms = c(ab_P_MS_delta),
      delta_msgs = c(ab_MS_GS_delta),
      delta_total = c(dt),
      diff_from_p= c(p_minus_p, m_minus_p, g_minus_p)
    )
    # Add BG summary to summary stats for all BG
    summary_stats <- rbind(summary_stats, temp)
  }
}

# Convert time to factor
summary_stats$group <- factor(summary_stats$group, levels=c("Preictal", "MS", "GS"))

# Get fold change between 2 periods
get_fold <- function(x, y) {
  ifelse(x > y, (x - y) / y, 
         ifelse(y > x, (y - x) / x, 0))
}

# Get dataframe with fold change for each behavior group
d1_fold <- summary_stats %>% dplyr::select(c(column, mean, group)) %>%
  pivot_wider(id_cols=column, names_from=group, values_from=mean) %>%
  # Perform rowise computations
  rowwise() %>%
  mutate(
    p_m_fold = get_fold(Preictal, MS),
    p_g_fold = get_fold(Preictal, GS),
    m_g_fold = get_fold(MS, GS),
    # Get average fold
    avg_fold = mean(c(p_m_fold, p_g_fold, m_g_fold), na.rm = TRUE)
  ) %>%
  # Arrange in desc(avg_fold) when writing to file
  arrange(desc(avg_fold))

#write_csv(d1_fold, "d1-pmg-fold.csv")

# Get data in format for statistical tests
a_data <- total %>% dplyr::select(!c(strain))
# The behavior groups
vars <- names(a_data)[-c(1, 65)]
# Paste bg so no issue with numeric BGs
nc <- paste0("bg_", vars)
new_cols <- c("id", nc, "time")
colnames(a_data) <- new_cols

# Get numeric data
num_v <- a_data %>% dplyr::select(!c(id, time))

# Shapiro Test for Normality
st <- shapiro_test(num_v, vars=names(a_data)[-c(1, 65)])

# Friedman Test for Usage~time

# Compares all 3
model_friedman <- function(var) {
  model <- friedman_test(as.formula(paste0(var, "~time|id")), data=a_data)
  return(model)
}

# Pairwise Wilcox Rank Sum for PvsMS, PvsGS, MSvsGs
# With bonferroni Correction
adj_confidence <- 1 - (0.05)/(63)
model_wilcox <- function(var) {
  model <- wilcox_test(data=a_data, formula=as.formula(paste0(var,"~time")), conf.level=adj_confidence, p.adjust.method="bonferroni", paired=T, 
                       detailed=T)
  return(model)
}



# Friedman Test
m_names <- names(a_data)[-c(1, 65)]
z <- lapply(m_names, model_friedman)
p_values <- sapply(z, function(x) x$p)
p_adjusted <- p.adjust(p_values, method="bonferroni")
friedman_df <- data.frame(behavior_group=m_names, p_value=p_values, p_adjusted=p_adjusted)


# Wilcox Test
wil_names <- filter(friedman_df, p_adjusted<0.05)$behavior_group
wil <- lapply(wil_names, model_wilcox)
wil_df <- bind_rows(wil)

pairwise_wilcox <- wil_df %>% dplyr::select(`.y.`, group1, group2, p, estimate, conf.low, conf.high) %>%
  mutate(estimate=round(estimate*100, 2), conf.low=round(conf.low*100, 2), conf.high=round(conf.high*100, 2))


#write_csv(friedman_df, "64_bgs_behavior_significant_friedman.csv")
#write_csv(wil_df, "63bg_behavior_pairwise_wilcox.csv")


# Set up df for plotting significance
friedman_df$bg <- gsub("^bg_", "", friedman_df$behavior_group)
friedman_df <- dplyr::select(friedman_df, c(p_adjusted, p_value, bg)) %>% mutate(significance = case_when(
  p_adjusted < 0.001 ~ "***",
  p_adjusted < 0.01 ~ "**",
  p_adjusted < 0.05 ~ "*", 
  TRUE ~ ""
))

summary_stats <- summary_stats %>% left_join(friedman_df, by=join_by(column==bg))


# Figure 2A
summary_stats %>% ggplot(aes(x=reorder(column, -delta_total), y=mean, color=group, ymin=mean-std_error, ymax=mean+std_error)) +
  geom_point(size=4) +
  geom_errorbar(width=0.1) +
  labs(x="Behavior Group", y="Average Percent Usage",
       color="Group") +
  theme_minimal() +
  scale_color_manual(values=c("#A0AF84", "#F0746E", "#A52794")) +
  geom_text(data=summary_stats %>% filter(group=="Preictal"),
            aes(x=reorder(column, -delta_total), label=significance, y=.15), color="black", fontface="bold", size=7, angle=90, vjust=.8, hjust="right") +
  theme(axis.text.x=element_text(angle=0, size=22, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=22, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=18),
        axis.ticks.length=unit(2, "mm"),
        legend.text=element_text(angle=0, size=15)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

#ggsave("d1_4rm_bg_histogram.tiff", dpi=300)


wil_df <- wil_df %>% rename(group=.y.) %>% mutate(group=as.factor(str_remove(group, "bg_")))

p_v_m <- wil_df %>% filter(group1=="Myoclonic Seizure" & group2=="Preictal")
p_v_m <- p_v_m %>% dplyr::select(c(p.adj, group))

p1 <- ggplot(p_v_m, aes(y=-log10(p.adj), x=reorder(group, log10(p.adj)))) +
  geom_hline(yintercept=-log10(0.05), color="red", lty=2) +
  geom_point(aes(color=p.adj<0.05)) +
  theme_minimal() +
  scale_color_manual(values=c("Black", "Red")) +
  labs(x="Behavior Group", y=expression(-log[10](p)), title="Preictal v. Myoclonic Seizure") +
  theme(axis.text.x=element_text(angle=0, size=22, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=22, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=18),
        axis.ticks.length=unit(2, "mm"),
        legend.position="none") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ylim(c(0, 25))


p1

p_v_g <- wil_df %>% filter(group1=="Generalized Seizure" & group2=="Preictal")
p_v_g <- p_v_g %>% dplyr::select(c(p.adj, group))

p2 <- ggplot(p_v_g, aes(y=-log10(p.adj), x=reorder(group, log10(p.adj)))) +
  geom_hline(yintercept=-log10(0.05), color="red", lty=2) +
  geom_point(aes(color=p.adj<0.05)) +
  theme_minimal() +
  scale_color_manual(values=c("Black", "Red")) +
  labs(x="Behavior Group", y=expression(-log[10](p)), title="Preictal v. Generalized Seizure") +
  theme(axis.text.x=element_text(angle=0, size=22, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=22, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=18),
        axis.ticks.length=unit(2, "mm"),
        legend.position="none") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ylim(c(0, 25))

m_v_g <- wil_df %>% filter(group1=="Generalized Seizure" & group2=="Myoclonic Seizure")
m_v_g <- m_v_g %>% dplyr::select(c(p.adj, group))

p3 <- ggplot(m_v_g, aes(y=-log10(p.adj), x=reorder(group, log10(p.adj)))) +
  geom_hline(lty=2, yintercept=-log10(0.05), color="red") +
  geom_point(aes(color=p.adj<0.05)) +
  theme_minimal() +
  scale_color_manual(values=c("black", "red")) +
  labs(x="Behavior Group", y=expression(-log[10](p)), title="Myoclonic Seizure v. Generalized Seizure") +
  theme(axis.text.x=element_text(angle=0, size=22, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=22, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=18),
        axis.ticks.length=unit(2, "mm"),
        legend.position="none") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ylim(c(0, 25))

library(patchwork)

p3
p2
p1
