# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
# 
# Fine Seizure Behavior Progresses Over Time During Repeated Seizure Challenges



# Load dependencies
library(tidyverse)
library(ggh4x)
library(tidyverse)
library(rstatix)
library(ggthemes)
library(gridExtra)
library(lme4)
library(lmerTest)
library(MASS)
library(pROC)
# Load data
all_day_usage <- read_csv("d3_all_days_usage.csv")
# Remove BGs
all_day_usage <- all_day_usage %>% dplyr::select(!c("53", "54", "57", "59"))

# Get summary stats for each day
get_day_means <- function(data) {
  summary_stats <- data.frame()
  for (col in names(data)) {
    if (col != "day" & col != "id") {  # Assuming the column specifying groups is named 'mortality'
      
      d1_mean <- mean(data[data$day == 1, col][[1]], na.rm = TRUE)
      d1_std_error <- sd(data[data$day == 1, col][[1]], na.rm = TRUE) / sqrt(sum(data$day == 1, na.rm=T))
      
      d2_mean <- mean(data[data$day == 2, col][[1]], na.rm = TRUE)
      d2_std_error <- sd(data[data$day == 2, col][[1]], na.rm = TRUE) / sqrt(sum(data$day == 2, na.rm=T))
      
      d3_mean <- mean(data[data$day == 3, col][[1]], na.rm = TRUE)
      d3_std_error <- sd(data[data$day == 3, col][[1]], na.rm = TRUE) / sqrt(sum(data$day == 3, na.rm=T))
      
      d4_mean <- mean(data[data$day == 4, col][[1]], na.rm = TRUE)
      d4_std_error <- sd(data[data$day == 4, col][[1]], na.rm = TRUE) / sqrt(sum(data$day == 4, na.rm=T))
      
      d5_mean <- mean(data[data$day == 5, col][[1]], na.rm = TRUE)
      d5_std_error <- sd(data[data$day == 5, col][[1]], na.rm = TRUE) / sqrt(sum(data$day == 5, na.rm=T))
      
      d6_mean <- mean(data[data$day == 6, col][[1]], na.rm = TRUE)
      d6_std_error <- sd(data[data$day == 6, col][[1]], na.rm = TRUE) / sqrt(sum(data$day == 6, na.rm=T))
      
      d7_mean <- mean(data[data$day == 7, col][[1]], na.rm = TRUE)
      d7_std_error <- sd(data[data$day == 7, col][[1]], na.rm = TRUE) / sqrt(sum(data$day == 7, na.rm=T))
      
      d8_mean <- mean(data[data$day == 8, col][[1]], na.rm = TRUE)
      d8_std_error <- sd(data[data$day == 8, col][[1]], na.rm = TRUE) / sqrt(sum(data$day == 8, na.rm=T))
      
      # Get absolute delta between day 8 and day 1
      dt=abs(d8_mean-d1_mean)
      # Get delta between day 8 and day 1
      d8_d1=(d8_mean-d1_mean)
      
      temp <- data.frame(
        column = col,
        group = c("1", "2", "3", "4", "5", "6", "7", "8"),
        mean = c(d1_mean, d2_mean, d3_mean, d4_mean, d5_mean, d6_mean, d7_mean, d8_mean),
        std_error = c(d1_std_error, d2_std_error, d3_std_error, d4_std_error, d5_std_error, d6_std_error, d7_std_error, d8_std_error),
        delta_total = c(dt),
        d8_d1_dt=c(d8_d1)
      )
      summary_stats <- rbind(summary_stats, temp)
    }
  }
  return(summary_stats)
}

cc051_p_data <- all_day_usage %>% filter(time=="Preictal" & strain=="CC051") %>% dplyr::select(!c(strain, time))
cc051_m_data <- all_day_usage %>% filter(time=="Myoclonic Seizure" & strain=="CC051") %>% dplyr::select(!c(strain, time))
cc051_g_data <- all_day_usage %>% filter(time=="Generalized Seizure" & strain=="CC051") %>% dplyr::select(!c(strain, time))

cc051_p_summary <- get_day_means(cc051_p_data) %>% mutate(strain="CC051")
cc051_m_summary <- get_day_means(cc051_m_data) %>% mutate(strain="CC051")
cc051_g_summary <- get_day_means(cc051_g_data) %>% mutate(strain="CC051")

cc051_p_data <- cc051_p_data %>% relocate(day)
cc051_m_data <- cc051_m_data %>% relocate(day)
cc051_g_data <- cc051_g_data %>% relocate(day)

# C57
c57_p_data <- all_day_usage %>% filter(time=="Preictal" & strain=="C57") %>% dplyr::select(!c(strain, time))
c57_m_data <- all_day_usage %>% filter(time=="Myoclonic Seizure" & strain=="C57") %>% dplyr::select(!c(strain, time))
c57_g_data <- all_day_usage %>% filter(time=="Generalized Seizure" & strain=="C57") %>% dplyr::select(!c(strain, time))

c57_p_summary <- get_day_means(c57_p_data) %>% mutate(strain="C57")
c57_m_summary <- get_day_means(c57_m_data) %>% mutate(strain="C57")
c57_g_summary <- get_day_means(c57_g_data) %>% mutate(strain="C57")

c57_p_data <- c57_p_data %>% relocate(day)
c57_m_data <- c57_m_data %>% relocate(day)
c57_g_data <- c57_g_data %>% relocate(day)


# Logit function to do logit transformation on data
logit <- function(p) {
  return(log(p/(1-p)))
}

# Mixed effects model for variable~day within id
adj_level <- 1 - (0.05)/(63)

mixed_effects_model <- function(var, data) {
  model <- lmer(as.formula(paste0(var, "~day + (1|id)")), data=data)
  model_summary <- summary(model)$coefficients %>% data.frame() %>% slice(-1) %>% rename(p=`Pr...t..`) %>% mutate(bg=var)
  model_confint <- confint(model, "day", level=adj_level)
  model_summary$lcl <- model_confint[1]
  model_summary$ucl <- model_confint[2]
  return(model_summary)
}

lmer_stats <- function(data, summ) {
  set.seed(1)
  a_data <- data
  m_names <- names(data)[-c(1,2)]
  nc <- paste0('bg_', m_names)
  new_cols <- c("day", "id", nc)
  colnames(a_data) <- new_cols
  above_mean <- summ %>% filter(delta_total>mean(delta_total)) %>% dplyr::select(column) %>% unique() %>% unlist() %>% as.vector()
  above_mean <- paste0("bg_", above_mean)
  
  epsilon <- 0.001
  above_data <- a_data %>% dplyr::select(id, day, all_of(above_mean)) %>% filter(day!=28)
  lg <- above_data %>% mutate(
    across(above_mean, ~ifelse(.==0, epsilon, ifelse(.==1, 1-epsilon, .)))
  ) %>% mutate(across(above_mean, logit))
  
  stat_result <- lapply(above_mean, mixed_effects_model, data=lg)
  stat_result <- bind_rows(stat_result)
  rownames(stat_result) <- NULL
  stat_result$adjusted <- p.adjust(stat_result$p, method="bonferroni")
  
  return(stat_result)
}

bg_stat_result <- lmer_stats(cc051_p_data, cc051_p_summary)
m_bg_stat_result <- lmer_stats(cc051_m_data, cc051_m_summary)
g_bg_stat_result <- lmer_stats(cc051_g_data, cc051_g_summary)


c57_p_stat_result <- lmer_stats(c57_p_data, c57_p_summary)
c57_m_stat_result <- lmer_stats(c57_m_data, c57_m_summary)
c57_g_stat_result <- lmer_stats(c57_g_data, c57_g_summary)


#write_csv(m_bg_stat_result, "ms_d3_63bg_kindling_bg_mixed_effects_stats.csv")
#write_csv(g_bg_stat_result, "gs_d3_63bg_kindling_bg_mixed_effects_stats.csv")

#write_csv(c57_m_stat_result, "ms_c57_kindling_mixed_effects_stats.csv")
#write_csv(c57_g_stat_result, "gs_c57_kindling_mixed_effects_stats.csv")

# Figure 6B

cc051_m_plotting <- cc051_m_summary %>% filter(delta_total>mean(delta_total)) %>% 
  mutate(x_var=reorder(interaction(group, column), d8_d1_dt))

c51_m_x_levels <- levels(cc051_m_plotting$x_var)
columns <- sapply(strsplit(c51_m_x_levels, "\\."), function(x) x[2])

column_changes <- which(columns[-1] != columns[-length(columns)]) + 0.5

cc051_m_plotting %>% ggplot(aes(x=x_var, y=mean, group=column)) +
  geom_point() +
  geom_smooth(method="lm", se=F, color="red") +
  scale_x_discrete(guide=guide_axis_nested()) +
  geom_vline(xintercept=column_changes, linetype="dashed", color="lightgrey", alpha=.5) +
  theme_classic() +
  theme(panel.grid=element_blank(), 
        axis.line = element_line(), 
        axis.ticks=element_line(),
        axis.text.y=element_text(size=20, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text.x=element_blank(),
        ggh4x.axis.nesttext.x=element_text(size=20, angle=90),
        ) +
  labs(x="Behavior Group for 8 Day Kindling", y="Percent Usage", title="CC051 Myoclinic Seizure Kindling")


# Figure 6C
cc051_g_plotting <- cc051_g_summary %>% filter(delta_total>mean(delta_total)) %>% 
  mutate(x_var=reorder(interaction(group, column), d8_d1_dt))

x_levels <- levels(cc051_g_plotting$x_var)
columns <- sapply(strsplit(x_levels, "\\."), function(x) x[2])

column_changes <- which(columns[-1] != columns[-length(columns)]) + 0.5

cc051_g_plotting %>% ggplot(aes(x_var, y=mean, group=column)) +
  geom_point() +
  geom_smooth(method="lm", se=F, color="red") +
  scale_x_discrete(guide=guide_axis_nested()) +
  geom_vline(xintercept=column_changes, linetype="dashed", color='lightgrey', alpha=.5) +
  theme_classic() +
  labs(x="Behavior Group for 8 Day Kindling", y="Percent Usage", title="CC051 Generalized Seizure Kindling") +
  theme(panel.grid=element_blank(), 
        axis.line = element_line(), 
        axis.ticks=element_line(),
        axis.text.y=element_text(size=20, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text.x=element_blank(),
        ggh4x.axis.nesttext.x=element_text(size=20, angle=90))

# C57
c57_m_plotting <- c57_m_summary %>% filter(delta_total>mean(delta_total)) %>% 
  mutate(x_var=reorder(interaction(group, column), d8_d1_dt))

x_levels <- levels(c57_m_plotting$x_var)
columns <- sapply(strsplit(x_levels, "\\."), function(x) x[2])

column_changes <- which(columns[-1] != columns[-length(columns)]) + 0.5


c57_m_plotting %>% ggplot(aes(x=x_var, y=mean, group=column)) +
  geom_point() +
  geom_smooth(method="lm", se=F, color="red") +
  scale_x_discrete(guide=guide_axis_nested()) +
  geom_vline(xintercept=column_changes, linetype="dashed", color="lightgray", alpha=.5) +
  theme_classic() +
  theme(panel.grid=element_blank(), 
        axis.line = element_line(), 
        axis.ticks=element_line(),
        axis.text.y=element_text(size=20, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text.x=element_blank(),
        ggh4x.axis.nesttext.x=element_text(size=20, angle=90)) +
  labs(x="Behavior Group for 8 Day Kindling", y="Percent Usage", title="B6J Myoclinic Seizure Kindling")


c57_g_plotting <- c57_g_summary %>% filter(delta_total>mean(delta_total)) %>% 
  mutate(x_var=reorder(interaction(group, column), d8_d1_dt))

x_levels <- levels(c57_g_plotting$x_var)
columns <- sapply(strsplit(x_levels, "\\."), function(x) x[2])

column_changes <- which(columns[-1] != columns[-length(columns)]) + 0.5

c57_g_plotting %>% ggplot(aes(x=x_var, y=mean, group=column)) +
  geom_point() +
  geom_smooth(method="lm", se=F, color="red") +
  scale_x_discrete(guide=guide_axis_nested()) +
  geom_vline(xintercept=column_changes, linetype="dashed", color='lightgray', alpha=.5) +
  theme_classic() +
  labs(x="Behavior Group for 8 Day Kindling", y="Percent Usage", title="B6J Generalized Seizure Kindling") +
  theme(panel.grid=element_blank(), 
        axis.line = element_line(), 
        axis.ticks=element_line(),
        axis.text.y=element_text(size=20, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text.x=element_blank(),
        ggh4x.axis.nesttext.x=element_text(size=20, angle=90))


# LDA for kindling
cc051_full <- all_day_usage %>% filter(strain=="CC051" & day!=28)
cc051_m <- cc051_full %>% filter(time=="Myoclonic Seizure") %>% dplyr::select(!c(strain, time)) %>% relocate(c(id, day))


z <- cc051_m_plotting %>% dplyr::select(column) %>% drop() %>% unique() %>% unlist()

# MS
cc051_m_numeric <- cc051_m %>% dplyr::select(all_of(z))
days_m <- cc051_m$day


set.seed(1)
lda_model <- lda(x=cc051_m_numeric, grouping=days_m)
lda_preds <- predict(lda_model)
ld_scores <- data.frame(lda_preds$x, Day=days_m)

# Get centroids
centroids_ms <- ld_scores %>%
  group_by(Day) %>%
  summarise(across(starts_with("LD"), mean))

# Plot LDA results
ggplot(ld_scores, aes(x=LD1, y=LD2, color=Day)) +
  geom_point(size=3) +
  scale_color_viridis_c(direction=-1, option="C") +
  theme_clean() +
  labs(x="LD1 (45.69%)", y="LD2 (23.52%)", color="Day", title="CC051 MS Kindling LDA") +
  theme(panel.grid.major.y = element_blank())
  

ggplot(centroids_ms, aes(x=LD1, y=LD2, color=Day)) +
  geom_point(size=3) +
  scale_color_viridis_c(direction=-1, option="C") +
  theme_classic() +
  labs(color="Day", title="CC051 MS Kindling LDA Centroids") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial"))

# GS
cc051_g <- cc051_full %>% filter(time=="Myoclonic Seizure") %>% dplyr::select(!c(strain, time)) %>% relocate(c(id, day))

z2 <- cc051_g_plotting %>% dplyr::select(column) %>% drop() %>% unique() %>% unlist()

cc051_g_numeric <- cc051_g %>% dplyr::select(all_of(z2))
days_g <- cc051_g$day

lda_model_gs <- lda(x=cc051_g_numeric, grouping=days_g)
lda_preds_gs <- predict(lda_model_gs)
ld_scores_gs <- data.frame(lda_preds_gs$x, Day=days_g)

centroids_gs <- ld_scores_gs %>%
  group_by(Day) %>%
  summarise(across(starts_with("LD"), mean))

ggplot(ld_scores_gs, aes(x=LD1, y=LD2, color=Day)) +
  geom_point(size=3) + 
  scale_color_viridis_c(direction=-1, option="C") +
  theme_clean() +
  labs(x="LD1 (52.92%)", y="LD2 (15.69%)", color="Day", title="CC051 GS Kindling LDA") +
  theme(panel.grid.major.y = element_blank())


ggplot(centroids_gs, aes(x=LD1, y=LD2, color=Day)) +
  geom_point(size=3) +
  scale_color_viridis_c(direction=-1, option="C") +
  theme_classic() +
  labs(color="Day", title="CC051 GS Kindling LDA Centroids") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial"))


# Binary Early vs Late Days LDA
days_m_bin <- ifelse(days_m<=5, "Early", "Late")
days_g_bin <- ifelse(days_g<=5, "Early", "Late")

# MS
set.seed(1)
lda_ms_bin <- lda(x=cc051_m_numeric, grouping=days_m_bin, CV=T)
mean(lda_ms_bin$class == days_m_bin)


m_bin_conf <- table(Predicted = lda_ms_bin$class, Actual = days_m_bin)
m_bin_conf_df <- as.data.frame(m_bin_conf)
colnames(m_bin_conf_df) <- c("Predicted", "Actual", "Count")

ggplot(m_bin_conf_df, aes(x = Actual, y = Predicted, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "darkred") +
  labs(title = "MS Early/Late LDA Confusion Matrix", x = "Actual", y = "Predicted") +
  theme_clean() +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial"), legend.position="none")


roc(response=ifelse(days_m_bin=="Late", 1, 0), predictor=lda_ms_bin$posterior[,2], plot=T, xlab="False Positive Percentage",
    ylab="True Positive Percentage", auc=T, main="MS Early v Late ROC Curve")

# Get values for plotting
set.seed(1)
lda_ms_bin_full <- lda(x=cc051_m_numeric, grouping=days_m_bin)

lda_preds_ms_bin <- predict(lda_ms_bin_full)
ld_scores_ms_bin <- data.frame(lda_preds_ms_bin$x, Day=days_m_bin)

ggplot(ld_scores_ms_bin, aes(x=LD1, color=Day, y=0)) +
  geom_point() +
  labs(title="Binary Early vs. Late LDA for CC051 MS Kindling") +
  theme_classic()

ggplot(ld_scores_ms_bin, aes(x=LD1, fill=Day)) +
  geom_density(alpha=.8) +
  labs(title="Early vs. Late LDA for CC051 MS Kindling") +
  labs(y="Density") +
  theme_classic() +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial")) +
  xlim(c(-5, 8))

ms_bin_lda_wilcox <- wilcox_test(ld_scores_ms_bin, LD1~Day, p.adjust.method="bonferroni", detailed=T)
#write_csv(ms_bin_lda_wilcox, "MS LDA Early v Late LD1 Wilcox.csv")

# GS
set.seed(1)
lda_gs_bin <- lda(x=cc051_g_numeric, grouping=days_g_bin, CV=T)
mean(lda_gs_bin$class == days_g_bin)

g_bin_conf <- table(Predicted = lda_gs_bin$class, Actual = days_g_bin)
g_bin_conf_df <- as.data.frame(g_bin_conf)
colnames(g_bin_conf_df) <- c("Predicted", "Actual", "Count")

ggplot(g_bin_conf_df, aes(x = Actual, y = Predicted, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "darkred") +
  labs(title = "GS Early/Late LDA Confusion Matrix", x = "Actual", y = "Predicted") +
  theme_clean() +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial"), legend.position="none")

roc(response=ifelse(days_g_bin=="Late", 1, 0), predictor=lda_gs_bin$posterior[,2], plot=T, xlab="False Positive Percentage",
    ylab="True Positive Percentage", auc=T, main="GS Early v Late ROC Curve")


# Get values for plotting
set.seed(1)
lda_gs_bin_full <- lda(x=cc051_g_numeric, grouping=days_g_bin)
lda_preds_gs_bin <- predict(lda_gs_bin_full)
ld_scores_gs_bin <- data.frame(lda_preds_gs_bin$x, Day=days_g_bin)


ggplot(ld_scores_gs_bin, aes(x=LD1, color=Day, y=0)) +
  geom_point() +
  labs(title="Binary Early vs. Late LDA for CC051 GS Kindling") +
  theme_classic()

ggplot(ld_scores_gs_bin, aes(x=LD1, fill=Day)) +
  geom_density(alpha=.8) +
  labs(title="Early vs. Late LDA for CC051 GS Kindling") +
  theme_classic() +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial")) +
  xlim(c(-5, 5))

gs_bin_lda_wilcox <- wilcox_test(ld_scores_gs_bin, LD1~Day, p.adjust.method="bonferroni", detailed=T)

#write_csv(gs_bin_lda_wilcox, "GS LDA Early v Late LD1 Wilcox.csv")


#setwd('/Users/jadenthomas/Desktop/OSU/Gu Lab/Dataset 1 Analysis')

# P v G LDA
p_v_g <- read_csv("LD1_test_dataset_GS_vs_Preictal.csv")
p_v_g$`True Class` <- factor(p_v_g$`True Class`, levels=c("Preictal", "GS"))


ggplot(p_v_g, aes(x=LD1, fill=`True Class`)) +
  geom_density(alpha=0.8) +
  scale_fill_manual(values=c("#A0AF84", "#A52794"))+
  theme_classic() +
  labs(y="Density", fill="Seizure State", title="LD1 Density for test data") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial")) +
  xlim(c(-15, 12))

p_v_g <- p_v_g %>% rename(true_class=`True Class`)

pg_lda_wilcox <- wilcox_test(p_v_g, LD1~true_class, p.adjust.method="bonferroni")
#write_csv(pg_lda_wilcox, "LDA P v G LD1 Wilcox.csv")

table(Predicted = p_v_g$`Predicted Class`, Actual = p_v_g$true_class)

roc(response=ifelse(p_v_g$true_class=="GS", 1, 0), predictor=p_v_g$`P(GS)`, plot=T, xlab="False Positive Percentage",
    ylab="True Positive Percentage", auc=T, main="P v GS ROC Curve")

