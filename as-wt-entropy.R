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
library(ggpubr)

# Load Data
df_fatal <- read_csv("mouse_data4.csv") %>% mutate(fatal=(Score==7))

# Partition by 
preictal <- df_fatal %>% filter(scorer<=3600) %>% mutate(time="preictal")
myoclonic_seizure <- df_fatal %>% filter(scorer>3600 & scorer<=(3600+((GST-MST)*30))) %>% mutate(time="myoclonic seizure")
gen_seizure <- df_fatal %>% filter(scorer>(3600+((GST-MST)*30))) %>% mutate(time="generalized seizure")

# Get BG usage for each individual in df
get_individual_df <- function(df) {
  p <- df %>% group_by(id, Genotype) %>% count(BSOiD_labels)
  p2 <- p %>% pivot_wider(names_from=BSOiD_labels, values_from=n, id_cols=c(id, Genotype))
  p2[is.na(p2)] <- 0
  
  nums_wide <- p2[,-c(1,2)]/rowSums(p2[,-c(1,2)])
  nums_wide <- mutate(nums_wide, id=p2$id, genotype=p2$Genotype)
  nums_wide <- relocate(nums_wide, c("id", "genotype"))
  
  return(nums_wide)
}

# Get BG usage for each individual by seizure state
p <- get_individual_df(preictal) %>% mutate(time="Preictal") %>% dplyr::select(!c("53", "54", "57", "59"))
m <- get_individual_df(myoclonic_seizure) %>% mutate(time="Myoclonic Seizure") %>% dplyr::select(!c("53", "54", "57", "59"))
g <- get_individual_df(gen_seizure) %>% mutate(time="Generalized Seizure") %>% dplyr::select(!c("53", "54", "57", "59"))

t <- get_individual_df(df_fatal) %>% mutate(time="Full")

total <- bind_rows(p, m, g)
total[is.na(total)] <- 0

# Get entropy for a given df
get_entropy <- function(df) {
  x <- df %>% dplyr::select(!c(id, genotype, time))
  e <- rep(0, nrow(x))
  
  
  for (i in 1:nrow(x)) {
    l <- x[i,]
    # Select nonzero behavior group usage
    q <- l[l!=0]
    v <- 0
    for (j in 1:length(q)) {
      temp <- q[j]*log2(q[j])
      v = v + temp
    }
    e[i] <- -v
  }
  
  entropyDF <- data.frame(id=df$id, genotype=df$genotype, entropy=e)
  return(entropyDF)
}

# Get entropy for each seizure state
p_entropy <- get_entropy(p)
m_entropy <- get_entropy(m)
g_entropy <- get_entropy(g)

p_entropy$genotype <- factor(p_entropy$genotype, levels=c("WT", "AS"))
m_entropy$genotype <- factor(m_entropy$genotype, levels=c("WT", "AS"))
g_entropy$genotype <- factor(g_entropy$genotype, levels=c("WT", "AS"))



# Mann Whitney U
p_wt <- wilcox_test(p_entropy, entropy~genotype, p.adjust.method="bonferroni", detailed=T)
m_wt <- wilcox_test(m_entropy, entropy~genotype, p.adjust.method="bonferroni", detailed=T)
g_wt <- wilcox_test(g_entropy, entropy~genotype, p.adjust.method="bonferroni", detailed=T)



# Kruskal Wallis
# p_kw <- kruskal_test(p_entropy, entropy~genotype)
# m_kw <- kruskal_test(m_entropy, entropy~genotype)
# g_kw <- kruskal_test(g_entropy, entropy~genotype)
# 
# write_csv(p_kw, "kruskal_d4_63bgs_aswt_p_entropy_stats.csv")
# write_csv(m_kw, "kruskal_d4_63bgs_aswt_m_entropy_stats.csv")
# write_csv(g_kw, "kruskal_d4_63bgs_aswt_g_entropy_stats.csv")

#write_csv(p_wt, "d4_63bgs_aswt_p_entropy_stats.csv")
#write_csv(m_wt, "d4_63bgs_aswt_m_entropy_stats.csv")
#write_csv(g_wt, "d4_63bgs_aswt_g_entropy_stats.csv")

# Preictal entropy
ggplot(p_entropy, aes(x=genotype, y=entropy, fill=genotype)) +
  geom_violin(linewidth=1.5, draw_quantiles=c(.5)) +
  geom_jitter(width=.1, alpha=.3) +
  scale_fill_bmj() +
  theme_classic() +
  labs(title="Preictal Entropy for AS and WT", x="Genotype", y="Entropy") +
  theme(legend.position="none", plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text=element_text(size=12, family="Arial")) +
  scale_y_continuous(limits=c(0, 5), breaks=seq(0, 5, by=1))


# Figure 7E
ggplot(m_entropy, aes(x=genotype, y=entropy, fill=genotype)) +
  geom_violin(linewidth=1.5, draw_quantiles=c(0.5)) +
  geom_jitter(width=.1, alpha=0.3) +
  scale_fill_bmj() +
  theme_classic() +
  labs(title="Myoclonic Seizure Entropy for AS and WT", x="Genotype", y="Entropy") +
  theme(legend.position="none", plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text=element_text(size=12, family="Arial")) +
  scale_y_continuous(limits=c(0, 5), breaks=seq(0, 5, by=1))

# Figure 7G
ggplot(g_entropy, aes(x=genotype, y=entropy, fill=genotype)) +
  geom_violin(linewidth=1.5, draw_quantiles=c(0.5)) +
  geom_jitter(width=.1, alpha=0.3) +
  scale_fill_bmj() +
  theme_classic() +
  labs(title="Generalized Seizure Entropy for AS and WT", x="Genotype", y="Entropy") +
  theme(legend.position="none", plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text=element_text(size=12, family="Arial")) +
  scale_y_continuous(limits=c(0, 5), breaks=seq(0, 5, by=1))
