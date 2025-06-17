# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
#
# Fine behavior discrimination of seizure states section

library(tidyverse)
library(rstatix)
library(lme4)
library(multcomp)
library(extrafont)
library(ggthemes)
library(gridExtra)

df <- read_csv("mouse_data.csv")
df_fatal <- df %>% mutate(fatal=(Score==7))

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


# Function to get entropy for given dataframe of behavior groups as columns and observations as rows
get_entropy <- function(df) {
  x <- df %>% dplyr::select(c(3:65))
  # Initialize entropy vector
  e <- rep(0, nrow(x))
  
  # For every observation
  for (i in 1:nrow(x)) {
    # Get row of behavior group usage percentage for observation i
    l <- x[i,]
    # Select behavior groups with nonzero usage
    q <- l[l!=0]
    # Initialize entropy value
    v <- 0
    # For every nonzero behavior group
    for (j in 1:length(q)) {
      # p * log2(p)
      temp <- q[j]*log2(q[j])
      # Sum to get entropy
      v = v + temp
    }
    # Entropy for observation i is v
    e[i] <- -v
  }
  
  # Return df with id strain and entropy for 1:nrow(x)
  entropyDF <- data.frame(id=df$id, strain=df$strain, entropy=e)
  return(entropyDF)
}

p <- get_individual_df(preictal) %>% mutate(time="Preictal") %>% dplyr::select(!c("53", "54", "57", "59"))
m <- get_individual_df(myoclonic_seizure) %>% mutate(time="Myoclonic Seizure") %>% dplyr::select(!c("53", "54", "57", "59"))
g <- get_individual_df(gen_seizure) %>% mutate(time="Generalized Seizure") %>% dplyr::select(!c("53", "54", "57", "59"))

total <- bind_rows(p, m, g)


# Get entropy for P, MS, and GS for every mouse
p_entropyDF <- get_entropy(p)
m_entropy <- get_entropy(m)
gs_entropy <- get_entropy(g)

# Add time column
p_entropyDF <- p_entropyDF %>% mutate(time="P")
m_entropy <- m_entropy %>% mutate(time="MS")
gs_entropy <- gs_entropy %>% mutate(time="GS")

# Combine into 1 df
pmg_e <- bind_rows(p_entropyDF, m_entropy, gs_entropy)

#write_csv(pmg_e, "d1_64bg_entropy.csv")

pmg_e$time <- factor(pmg_e$time, c("P", "MS", "GS"))

# Figure 2F
ggplot(pmg_e, aes(x=time, y=entropy, fill=time)) +
  geom_violin(linewidth=1.5, draw_quantiles = c(.5)) +
  labs(x="Seizure State", y="Entropy (bits)", title="Entropy by Seizure State") +
  geom_jitter(width=.1, alpha=0.3) +
  theme_classic() +
  scale_fill_manual(values=c("#A0AF84", "#F0746E", "#A52794")) +
  theme(legend.position="none", plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text=element_text(size=12, family="Arial"))

#ggsave("d1_4rm_entropy_violin.tiff", dpi=300)

# Get entropy for Figure 2E
entropy_obs <- df %>% filter(id=="CC012-2")
entropy_obs_p <- preictal %>% filter(id=="CC012-2")
entropy_obs_ms <- myoclonic_seizure %>% filter(id=="CC012-2")
entropy_obs_gs <- gen_seizure %>% filter(id=="CC012-2")

# Get entropy by time for CC012-2
#pmg_e %>% filter(id=="CC012-2") %>% print()

# Full plot of bg across time to get consistent color
full <- entropy_obs %>% filter(scorer>=3421 & scorer<=3960) %>% ggplot(aes(x=scorer, y=BSOiD_labels, color=as.factor(BSOiD_labels))) +
  geom_point() +
  geom_vline(xintercept=3600) +
  geom_vline(xintercept=3780) +
  labs(x="Frame", y="Behavior Group", title="Mouse CC012-2") +
  theme_minimal() +
  theme(legend.position="none", plot.title=element_text(family="Arial", size=20))

d <- ggplot_build(full)$data[[1]]
# Get vector of color by point
cols <- d$colour

# Get last 30 seconds of P for CC012-2
s <- entropy_obs_p %>% filter(scorer>=3421)

# Get colors to plot
p_cols <- cols[1:180]

# Change removed BGs to grey
p_cols[(s$BSOiD_labels=="53" | s$BSOiD_labels=="54" | s$BSOiD_labels=="57" | s$BSOiD_labels=="59")] <- "gray"

# Get full MS
s <- entropy_obs_ms

# Get colors to plot and change removed BGs to grey
ms_cols <- cols[181:360]
ms_cols[(s$BSOiD_labels=="53" | s$BSOiD_labels=="54" | s$BSOiD_labels=="57" | s$BSOiD_labels=="59")] <- "gray"

# Get first 30 seconds of GS
s <- entropy_obs_gs %>% filter(scorer<=3960)
# Get colors and change removed Bgs to grey
gs_cols <- cols[361:540]
gs_cols[(s$BSOiD_labels=="53" | s$BSOiD_labels=="54" | s$BSOiD_labels=="57" | s$BSOiD_labels=="59")] <- "gray"

# Get separate plots to join
p1 <- entropy_obs_p %>% filter(scorer>=3421) %>% ggplot(aes(x=scorer, y=BSOiD_labels, color=p_cols)) +
  geom_point() +
  scale_color_identity() +
  ylim(0, 66) +
  labs(y="Behavior Group", title="Mouse CC012-2 Preictal: 4.21 bits") +
  theme_minimal() +
  scale_x_continuous(breaks=c(3421, 3466, 3511, 3556, 3600), labels=c(0, 45, 90, 135, 180)) +
  theme(legend.position = "none", plot.title=element_text(family="Arial", size=15), 
        panel.border = element_rect(fill=NA, linewidth=1), axis.ticks = element_line(linewidth=1), 
        panel.grid=element_blank(),         
        axis.title.y = element_text(size=14),
        axis.text=element_text(size=12, family="Arial"),
        axis.title.x=element_blank())

p2 <- ggplot(entropy_obs_ms, aes(x=scorer, y=BSOiD_labels, color=ms_cols)) +
  geom_point() +
  scale_color_identity() +
  ylim(0, 66) +
  labs(y="Behavior Group", title="Mouse CC012-2 MS: 0.68 bits") +
  theme_minimal() +
  scale_x_continuous(breaks=c(3601, 3645, 3690, 3735, 3780), labels=c(0, 45, 90, 135, 180)) +
  theme(legend.position = "none", plot.title=element_text(family="Arial", size=15), 
        panel.border = element_rect(fill=NA, linewidth=1), axis.ticks = element_line(linewidth=1), 
        panel.grid=element_blank(), 
        axis.title.y = element_text(size=14),
        axis.text=element_text(size=12, family="Arial"),
        axis.title.x=element_blank())

p3 <- entropy_obs_gs %>% filter(scorer<=3960) %>% ggplot(aes(x=scorer, y=BSOiD_labels, color=gs_cols)) +
  geom_point() +
  scale_color_identity() +
  ylim(0, 66) +
  labs(y="Behavior Group", title="Mouse CC012-2 GS: 1.72 bits") +
  theme_minimal() +
  scale_x_continuous(breaks=c(3781, 3826, 3871, 3916, 3961), labels=c(0, 45, 90, 135, 180)) +
  theme(legend.position = "none", plot.title=element_text(family="Arial", size=15), 
        panel.border = element_rect(fill=NA, linewidth=1), axis.ticks = element_line(linewidth=1, ), 
        panel.grid=element_blank(), 
        axis.title.y = element_text(size=14),
        axis.text=element_text(size=12, family="Arial"),
        axis.title.x=element_blank())

# Figure 2E
l <- grid.arrange(p1, p2, p3, ncol=3)
l <- arrangeGrob(p1, p2, p3, ncol=3)
#ggsave("d1_4rm_entropy_example.tiff", l, dpi=300)

# Test for normality
shapiro_test(pmg_e, entropy)

# Friedman and post-hoc wilcox test
f <- friedman_test(pmg_e, entropy~time|id)
wt <- wilcox_test(pmg_e, entropy~time, paired=T, p.adjust.method = "bonferroni", detailed=T) %>% as.data.frame()

#write_csv(f, "d1_63bg_entropy_friedman.csv")
#write_csv(wt, "d1_63bg_entropy_wilcox.csv")
