# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
#
# Fine Seizure Behavior Progresses Over Time During Repeated Seizure Challenges

# Load Dependencies
library(tidyverse)
library(ggh4x)
library(ggthemes)
library(lme4)
library(lmerTest)
library(broom)

# Load Data
all_day_usage <- read_csv("d3_all_days_usage.csv")
# Remove BGs
all_day_usage <- all_day_usage %>% dplyr::select(!c("53", "54", "57", "59"))


# Get entropy for given df
get_entropy <- function(df) {
  x <- df %>% dplyr::select(!c(id, strain, time, day))
  e <- rep(0, nrow(x))
  
  
  for (i in 1:nrow(x)) {
    l <- x[i,]
    q <- l[l!=0]
    v <- 0
    for (j in 1:length(q)) {
      temp <- q[j]*log2(q[j])
      v = v + temp
    }
    e[i] <- -v
  }
  
  entropyDF <- data.frame(id=df$id, strain=df$strain, day=df$day, entropy=e)
  return(entropyDF)
}

# Function to get entropy for every day
get_days_entropy <- function(data) {
  d1_entropy <- data %>% filter(day==1) %>% get_entropy()
  d2_entropy <- data %>% filter(day==2) %>% get_entropy()
  d3_entropy <- data %>% filter(day==3) %>% get_entropy()
  d4_entropy <- data %>% filter(day==4) %>% get_entropy()
  d5_entropy <- data %>% filter(day==5) %>% get_entropy()
  d6_entropy <- data %>% filter(day==6) %>% get_entropy()
  d7_entropy <- data %>% filter(day==7) %>% get_entropy()
  d8_entropy <- data %>% filter(day==8) %>% get_entropy()
  
  e_df <- bind_rows(d1_entropy, d2_entropy, d3_entropy, d4_entropy, d5_entropy, d6_entropy,
                    d7_entropy, d8_entropy)
  return(e_df)
}

# Get dfs with entropy with id, strain, day, and entropy columns
total_kindle_entropy <- all_day_usage %>% filter(time=="Total") %>% get_days_entropy()
p_kindle_entropy <- all_day_usage %>% filter(time=="Preictal") %>% get_days_entropy()
m_kindle_entropy <- all_day_usage %>% filter(time=="Myoclonic Seizure") %>% get_days_entropy()
g_kindle_entropy <- all_day_usage %>% filter(time=="Generalized Seizure") %>% get_days_entropy()

# Get strain averages and se
p_split_entropy <- p_kindle_entropy %>% group_by(strain, day) %>% mutate(m=mean(entropy), se=sd(entropy)/sqrt(n()))
m_split_entropy <- m_kindle_entropy %>% group_by(strain, day) %>% mutate(m=mean(entropy), se=sd(entropy)/sqrt(n()))
g_split_entropy <- g_kindle_entropy %>% group_by(strain, day) %>% mutate(m=mean(entropy), se=sd(entropy)/sqrt(n()))

# Split by strain C57
c57_p_entropy <- p_split_entropy %>% filter(strain=="C57") %>% mutate(time="Preictal")
c57_m_entropy <- m_split_entropy %>% filter(strain=="C57") %>% mutate(time="Myoclonic Seizure")
c57_g_entropy <- g_split_entropy %>% filter(strain=="C57") %>% mutate(time="Generalized Seizure")

c57_combo <- bind_rows(c57_p_entropy, c57_m_entropy, c57_g_entropy) %>% mutate(time=factor(time, levels=c("Preictal", "Myoclonic Seizure", "Generalized Seizure")))


# Split by strain CC051
cc051_p_entropy <- p_split_entropy %>% filter(strain=="CC051") %>% mutate(time="Preictal")
cc051_m_entropy <- m_split_entropy %>% filter(strain=="CC051") %>% mutate(time="Myoclonic Seizure")
cc051_g_entropy <- g_split_entropy %>% filter(strain=="CC051") %>% mutate(time="Generalized Seizure")

cc051_combo <- bind_rows(cc051_p_entropy, cc051_m_entropy, cc051_g_entropy) %>% mutate(time=factor(time, levels=c("Preictal", "Myoclonic Seizure", "Generalized Seizure")))

set.seed(2024)


# Statistical test for entropy by day in id for C57

aov(entropy~as.factor(day), data=c57_p_entropy)

p_c57_etest <- lmer(entropy~as.factor(day) + (1|id), data=c57_p_entropy)
p_c57_ci <- confint(p_c57_etest, 3:10) %>% as.data.frame()
m_c57_etest <- lmer(entropy~as.factor(day) + (1|id), data=c57_m_entropy)
m_c57_ci <- confint(m_c57_etest, 3:10) %>% as.data.frame()
g_c57_etest <- lmer(entropy~as.factor(day) + (1|id), data=c57_g_entropy)
g_c57_ci <- confint(g_c57_etest, 3:10) %>% as.data.frame()

c57_p_anova <- anova(p_c57_etest)
c57_m_anova <- anova(m_c57_etest)
c57_g_anova <- anova(g_c57_estest)

#write_csv(c57_p_anova, "c57_kindling_entropy_p_anova.csv")
#write_csv(c57_m_anova, "c57_kindling_entropy_m_anova.csv")
#write_csv(c57_g_anova, "c57_kindling_entropy_g_anova.csv")

p_summ_c57 <- summary(p_c57_etest)$coefficients %>% data.frame() %>% rename(p=`Pr...t..`) %>% mutate(lcl = p_c57_ci[,1], ucl = p_c57_ci[,2])
m_summ_c57 <- summary(m_c57_etest)$coefficients %>% data.frame() %>% rename(p=`Pr...t..`) %>% mutate(lcl = m_c57_ci[,1], ucl = m_c57_ci[,2])
g_summ_c57 <- summary(g_c57_estest)$coefficients %>% data.frame() %>% rename(p=`Pr...t..`) %>% mutate(lcl = g_c57_ci[,1], ucl = g_c57_ci[,2])

#write_csv(p_summ_c57, "c57_kindling_entropy_p_lmersummary.csv")
#write_csv(m_summ_c57, "c57_kindling_entropy_m_lmersummary.csv")
#write_csv(g_summ_c57, "c57_kindling_entropy_g_lmersummary.csv")

# Statistical test for entropy by day in id
p_entropy_test <- lmer(entropy~as.factor(day) + (1|id), data=cc051_p_entropy)
p_entropy_ci <- confint(p_entropy_test, 3:10) %>% as.data.frame()

m_entropy_test <- lmer(entropy~as.factor(day) + (1|id), data=cc051_m_entropy)
m_entropy_ci <- confint(m_entropy_test, 3:10) %>% as.data.frame()

g_entropy_test <- lmer(entropy~as.factor(day) + (1|id), data=cc051_g_entropy)
g_entropy_ci <- confint(g_entropy_test, 3:10) %>% as.data.frame()


anova(p_entropy_test)
anova(m_entropy_test)
anova(g_entropy_test)

p_summ <- summary(p_entropy_test)$coefficients %>% data.frame() %>% rename(p=`Pr...t..`) %>% mutate(lcl = p_entropy_ci[,1], ucl = p_entropy_ci[,2])
m_summ <- summary(m_entropy_test)$coefficients %>% data.frame() %>% rename(p=`Pr...t..`) %>% mutate(lcl = m_entropy_ci[,1], ucl = m_entropy_ci[,2])
g_summ <- summary(g_entropy_test)$coefficients %>% data.frame() %>% rename(p=`Pr...t..`) %>% mutate(lcl = g_entropy_ci[,1], ucl = g_entropy_ci[,2])

#write_csv(p_summ, "d4_kindling_entropy_p_lmersummary.csv")
#rite_csv(m_summ, "d4_kindling_entropy_m_lmersummary.csv")
#write_csv(g_summ, "d4_kindling_entropy_g_lmersummary.csv")

# Join MS and GS
mg_summary <- bind_rows(m_summ, g_summ)
mg_summary$p_adjusted <- p.adjust(mg_summary$p)

c57mg_summary <- bind_rows(m_summ_c57, g_summ_c57) %>% mutate(p_adjusted=p.adjust(p))

#write_csv(mg_summary, "d3_63bgs_kindling_entropy.csv")

# Figure 6D
p1 <- ggplot(cc051_combo, aes(x=as.factor(day), y=m, fill=time, ymin=m-se, ymax=m+se, color=time, group=time)) +
  geom_point() +
  geom_line() +
  #facet_wrap(vars(time)) +
  scale_color_manual(values=c("#A0AF84", "#F0746E", "#A52794")) +
  scale_fill_manual(values=c("#A0AF84", "#F0746E", "#A52794")) +
  geom_ribbon(alpha=.1) +
  theme_clean() +
  labs(x="Day", y="Entropy", color="Time", fill="Time", title="CC051 Kindling") +
  theme(axis.text.x=element_text(size=12, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=12, family="Arial"),
        axis.title = element_text(size=14),
        strip.text=element_text(family="Arial", size=20),
        panel.grid.major.y=element_blank()) +
  ylim(c(0, 4.5))


# C57 Entropy
p2 <- ggplot(c57_combo, aes(x=as.factor(day), y=m, fill=time, ymin=m-se, ymax=m+se, color=time, group=time)) +
  geom_point() +
  geom_line() +
  #facet_wrap(vars(time)) +
  scale_color_manual(values=c("#A0AF84", "#F0746E", "#A52794")) +
  scale_fill_manual(values=c("#A0AF84", "#F0746E", "#A52794")) +
  geom_ribbon(alpha=.1) +
  theme_clean() +
  labs(x="Day", y="Entropy", color="Time", fill="Time", title="B6J Kindling") +
  theme(axis.text.x=element_text(size=12, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=12, family="Arial"),
        axis.title = element_text(size=14),
        strip.text=element_text(family="Arial", size=20),
        panel.grid.major.y=element_blank()) +
  ylim(c(0, 4.5))

library(patchwork)
p1 + p2 + plot_layout(guides="collect")

#ggsave("cc051_kindling_entropy.tiff", p1, device="tiff", height=563, width=795, units="px", dpi=72)
