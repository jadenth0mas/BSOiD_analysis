# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
#
# Behavior ID of Mouse Strains

# Load Packages
library(tidyverse)
library(rstatix)
library(PMCMRplus)
library(MASS)
library(ggrepel)
library(patchwork)


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

total <- total %>% dplyr::select(!c("53", "54", "57", "59"))

total %>% filter(time=="Myoclonic Seizure") %>% dplyr::select(strain) %>% drop() %>% table()

# Coat Color LDA
p_data <- total %>% filter(time=="Preictal")

p_cols <- preictal %>% dplyr::select(strain, Color) %>% unique()
p_cols <- p_cols %>% mutate(Color=ifelse(strain=="C57", "black", Color)) %>% drop_na()

color_fix_df <- mutate(p_cols, Color=case_when(str_detect(Color, "agouti") ~ "agouti", .default=as.character(Color)))
color_df <- color_fix_df %>% drop_na()

p_data <- p_data %>% left_join(color_df, by=join_by(strain==strain)) %>% dplyr::select(!c(strain, id, time))

p_data %>% group_by(Color) %>% summarize(across(`3`:`12`, \(x) mean(x, na.rm=T)))

# P
set.seed(1)
X <- p_data %>% dplyr::select(!c(Color, `12`))
color <- p_data$Color

preictal_strain_lda <- lda(x=X, grouping=color)
p_strain_predicts <- predict(preictal_strain_lda)
p_strain_scores <- data.frame(p_strain_predicts$x, Color=color)

centroids_p <- p_strain_scores %>%
  group_by(Color) %>%
  summarise(across(starts_with("LD"), mean))

ggplot(p_strain_scores, aes(x=LD1, y=LD2, color=as.factor(Color))) +
  geom_point(size=3) + 
  scale_color_manual(values=c("#b08968", "black", "gray"), labels=c("Agouti", "Black", "White")) +
  theme_classic() +
  labs(x="LD1 (61.65%)", y="LD2 (38.35%)", color="Coat Color", title="Preictal Coat Color LDA") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial"))


ggplot(centroids_p, aes(x=LD1, y=LD2, color=Color)) +
  geom_point(size=3) +
  scale_color_manual(values=c("#b08968", "black", "gray"), labels=c("Agouti", "Black", "White")) +
  theme_classic() +
  labs(color="Day", title="P Coat Color LDA Centroids") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial"))


# MS
m_data <- total %>% filter(time=="Myoclonic Seizure")

m_cols <- myoclonic_seizure %>% dplyr::select(strain, Color) %>% unique()
m_cols <- m_cols %>% mutate(Color=ifelse(strain=="C57", "black", Color)) %>% drop_na()

color_fix_df <- mutate(m_cols, Color=case_when(str_detect(Color, "agouti") ~ "agouti", .default=as.character(Color)))
color_df <- color_fix_df %>% drop_na()

m_data <- m_data %>% left_join(color_df, by=join_by(strain==strain)) %>% dplyr::select(!c(strain, id, time))

m_data %>% group_by(Color) %>% summarize(across(`3`:`12`, \(x) mean(x, na.rm=T)))


X <- m_data %>% dplyr::select(!c(Color))
color <- m_data$Color

set.seed(1)
m_strain_lda <- lda(x=X, grouping=color)

m_strain_predicts <- predict(m_strain_lda)
m_strain_scores <- data.frame(m_strain_predicts$x, Color=color)

centroids_m <- m_strain_scores %>%
  group_by(Color) %>%
  summarise(across(starts_with("LD"), mean))

ggplot(m_strain_scores, aes(x=LD1, y=LD2, color=as.factor(Color))) +
  geom_point(size=3) + 
  scale_color_manual(values=c("#b08968", "black", "gray"), labels=c("Agouti", "Black", "White")) +
  theme_classic() +
  labs(x="LD1 (58.02%)", y="LD2 (41.98%)", color="Coat Color", title="MS Coat Color LDA") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial"))


ggplot(centroids_m, aes(x=LD1, y=LD2, color=Color)) +
  geom_point(size=3) +
  scale_color_manual(values=c("#b08968", "black", "gray"), labels=c("Agouti", "Black", "White")) +
  theme_classic() +
  labs(color="Day", title="MS Coat Color LDA Centroids") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial"))


# GS
g_data <- total %>% filter(time=="Generalized Seizure")

g_cols <- gen_seizure %>% dplyr::select(strain, Color) %>% unique()
g_cols <- g_cols %>% mutate(Color=ifelse(strain=="C57", "black", Color)) %>% drop_na()

color_fix_df <- mutate(g_cols, Color=case_when(str_detect(Color, "agouti") ~ "agouti", .default=as.character(Color)))
color_df <- color_fix_df %>% drop_na()

g_data <- g_data %>% left_join(color_df, by=join_by(strain==strain)) %>% dplyr::select(!c(strain, id, time))

g_data %>% group_by(Color) %>% summarize(across(`3`:`12`, \(x) mean(x, na.rm=T)))


X <- g_data %>% dplyr::select(!c(Color))
color <- g_data$Color

set.seed(1)
g_strain_lda <- lda(x=X, grouping=color)

g_strain_predicts <- predict(g_strain_lda)
g_strain_scores <- data.frame(g_strain_predicts$x, Color=color)

centroids_g <- g_strain_scores %>%
  group_by(Color) %>%
  summarise(across(starts_with("LD"), mean))

ggplot(g_strain_scores, aes(x=LD1, y=LD2, color=as.factor(Color))) +
  geom_point(size=3) + 
  scale_color_manual(values=c("#b08968", "black", "gray"), labels=c("Agouti", "Black", "White")) +
  theme_classic() +
  labs(x="LD1 (66.09%)", y="LD2 (33.91%)", color="Coat Color", title="GS Coat Color LDA") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial"))

ggplot(centroids_g, aes(x=LD1, y=LD2, color=Color)) +
  geom_point(size=3) +
  scale_color_manual(values=c("#b08968", "black", "gray"), labels=c("Agouti", "Black", "White")) +
  theme_classic() +
  labs(color="Day", title="GS Coat Color LDA Centroids") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial"))


# All strains LDA

# P
p_all <- p %>% dplyr::select(!c(id, time, `12`))
p_all$strain <- ifelse(p_all$strain=="C57", "B6J", p_all$strain)

set.seed(1)
p_all_lda <- lda(strain~., data=p_all)
p_all_predicts <- predict(p_all_lda)
p_all_scores <- data.frame(p_all_predicts$x, Strain=p_all$strain)

centroids_p_all <- p_all_scores %>%
  group_by(Strain) %>%
  summarise(across(starts_with("LD"), mean))

p_all_1 <- ggplot(p_all_scores, aes(x=LD1, y=LD2, color=as.factor(Strain))) +
  geom_point(size=3) + 
  theme_classic() +
  labs(x="LD1 (20.50%)", y="LD2 (12.26%)", color="Strain", title="P Strains LDA") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial")) +
  lims(y=c(-6, 6), x=c(-22, 5))

p_all_2 <- ggplot(centroids_p_all, aes(x=LD1, y=LD2, color=Strain)) +
  geom_point(size=3) +
  theme_classic() +
  labs(color="Strain", title="P Strains LDA") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial")) +
  lims(y=c(-6, 6), x=c(-22, 5))

p_all_1 + p_all_2 + plot_layout(guides="collect")

# MS
ms_all <- m %>% dplyr::select(!c(id, time))
ms_all$strain <- ifelse(ms_all$strain=="C57", "B6J", ms_all$strain)

set.seed(1)
ms_all_lda <- lda(strain~., data=ms_all)
ms_all_predicts <- predict(ms_all_lda)
ms_all_scores <- data.frame(ms_all_predicts$x, Strain=ms_all$strain)

centroids_ms_all <- ms_all_scores %>%
  group_by(Strain) %>%
  summarise(across(starts_with("LD"), mean))

ms_all_1 <- ggplot(ms_all_scores, aes(x=LD1, y=LD2, color=as.factor(Strain))) +
  geom_point(size=3) + 
  theme_classic() +
  labs(x="LD1 (30.77%)", y="LD2 (8.69%)", color="Strain", title="MS Strains LDA") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial")) +
  lims(x=c(-5, 30), y=c(-12, 5))

ms_all_2 <- ggplot(centroids_ms_all, aes(x=LD1, y=LD2, color=Strain)) +
  geom_point(size=3) +
  theme_classic() +
  labs(color="Strain", title="MS Strains LDA") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial")) +
  lims(x=c(-5, 30), y=c(-12, 5))

ms_all_1 + ms_all_2 + plot_layout(guides="collect")


# GS
gs_all <- g %>% dplyr::select(!c(id, time))
gs_all$strain <- ifelse(gs_all$strain=="C57", "B6J", gs_all$strain)


set.seed(1)
gs_all_lda <- lda(strain~., data=gs_all)
gs_all_predicts <- predict(gs_all_lda)
gs_all_scores <- data.frame(gs_all_predicts$x, Strain=gs_all$strain)

centroids_gs_all <- gs_all_scores %>%
  group_by(Strain) %>%
  summarise(across(starts_with("LD"), mean))

gs_all_1 <- ggplot(gs_all_scores, aes(x=LD1, y=LD2, color=as.factor(Strain))) +
  geom_point(size=3) + 
  theme_classic() +
  labs(x="LD1 (12.04%)", y="LD2 (10.70%)", color="Strain", title="GS Strains LDA") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial")) +
  lims(x=c(-20, 5), y=c(-7, 8))

gs_all_2 <- ggplot(centroids_gs_all, aes(x=LD1, y=LD2, color=Strain)) +
  geom_point(size=3) +
  theme_classic() +
  labs(color="Strain", title="GS Strains LDA") +
  theme(plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=16),
        axis.text=element_text(size=14, family="Arial")) +
  lims(x=c(-20, 5), y=c(-7, 8))

gs_all_1 + gs_all_2 + plot_layout(guides="collect")
