# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
# 
# Behavior biomarker of SUDEP 

library(tidyverse)
library(gridExtra)
library(rstatix)

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

get_entropy <- function(df) {
  x <- df %>% dplyr::select(c(3:65))
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
  
  entropyDF <- data.frame(id=df$id, strain=df$strain, fatal=df$fatal, entropy=e)
  return(entropyDF)
}


p <- get_individual_df(preictal) %>% mutate(time="Preictal") %>% left_join(unique(dplyr::select(df_fatal, c(id, fatal)))) %>% dplyr::select(!c("53", "54", "57", "59"))
m <- get_individual_df(myoclonic_seizure) %>% mutate(time="Myoclonic Seizure") %>% left_join(unique(dplyr::select(df_fatal, c(id, fatal)))) %>% dplyr::select(!c("53", "54", "57", "59"))
g <- get_individual_df(gen_seizure) %>% mutate(time="Generalized Seizure") %>% left_join(unique(dplyr::select(df_fatal, c(id, fatal)))) %>% dplyr::select(!c("53", "54", "57", "59"))

total <- bind_rows(p, m, g)



p_entropy <- get_entropy(p) %>% mutate(fatal=factor(fatal, levels=c(TRUE, FALSE), labels=c("Fatal", "Non-Fatal")))
m_entropy <- get_entropy(m) %>% mutate(fatal=factor(fatal, levels=c(TRUE, FALSE), labels=c("Fatal", "Non-Fatal")))
g_entropy <- get_entropy(g) %>% mutate(fatal=factor(fatal, levels=c(TRUE, FALSE), labels=c("Fatal", "Non-Fatal")))

fatal_gu <- g_entropy %>% left_join(select(g, !c(time, fatal, strain)), by=join_by(id==id))

p_entropy <- mutate(p_entropy, time="Preictal")
m_entropy <- mutate(m_entropy, time="Myoclonic Seizure")
g_entropy <- mutate(g_entropy, time="Generalized Seizure")

pmg_e <- bind_rows(p_entropy, m_entropy, g_entropy)

#write_csv(pmg_e, "d1_fatal_63bg_entropy.csv")

pmg_e$time <- factor(pmg_e$time, c("P", "MS", "GS"))

p_wt <- wilcox_test(p_entropy, entropy~fatal, p.adjust.method="bonferroni", detailed=T)
m_wt <- wilcox_test(m_entropy, entropy~fatal, p.adjust.method="bonferroni", detailed=T)
g_wt <- wilcox_test(g_entropy, entropy~fatal, p.adjust.method="bonferroni", detailed=T)

# p_kt <- kruskal_test(p_entropy, entropy~fatal)
# m_kt <- kruskal_test(m_entropy, entropy~fatal)
# g_kt <- kruskal_test(g_entropy, entropy~fatal)


# write_csv(p_kt, "kruskal_d1_fatal_p_entropy_stats.csv")
# write_csv(m_kt, "kruskal_d1_fatal_m_entropy_stats.csv")
# write_csv(g_kt, "kruskal_d1_fatal_g_entropy_stats.csv")

#write_csv(p_wt, "d1_fatal_p_entropy_stats.csv")
#write_csv(m_wt, "d1_fatal_m_entropy_stats.csv")
#write_csv(g_wt, "d1_fatal_g_entropy_stats.csv")

ggplot(p_entropy, aes(x=factor(fatal, level=c("Non-Fatal", "Fatal")), y=entropy, fill=fatal)) +
  geom_violin(linewidth=1.5, draw_quantiles=c(.5)) +
  geom_jitter(width=.1, alpha=0.3) +
  scale_fill_nejm() +
  theme_classic() +
  labs(title="Preictal Entropy for Fatal and Non-Fatal", x="Fatality", y="Entropy") +
  theme(legend.position="none", plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text=element_text(size=12, family="Arial"))

ggplot(m_entropy, aes(x=factor(fatal, level=c("Non-Fatal", "Fatal")), y=entropy, fill=fatal)) +
  geom_violin(linewidth=1.5, draw_quantiles=c(0.5)) +
  geom_jitter(width=.1, alpha=0.3) +
  scale_fill_nejm() +
  theme_classic() +
  labs(title="Myoclonic Seizure Entropy for Fatal and Non-Fatal", x="Fatality", y="Entropy") +
  theme(legend.position="none", plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text=element_text(size=12, family="Arial"))

# Figure 4B
ggplot(g_entropy, aes(x=factor(fatal, level=c("Non-Fatal", "Fatal")), y=entropy, fill=fatal)) +
  geom_violin(linewidth=1.5, draw_quantiles=c(0.5)) +
  geom_jitter(width=.1, alpha=0.3) +
  scale_fill_nejm() +
  theme_classic() +
  labs(title="Generalized Seizure Entropy for Fatal and Non-Fatal", x="Fatality", y="Entropy") +
  theme(legend.position="none", plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text=element_text(size=12, family="Arial"))

