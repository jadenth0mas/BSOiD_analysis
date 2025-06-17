# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
# 
# Fine behavior discrimination of seizure states section

# Load Packages
library(tidyverse)
library(rstatix)
library(lme4)
library(multcomp)
library(extrafont)
library(ggthemes)
library(syndRomics)

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

p <- get_individual_df(preictal) %>% mutate(time="Preictal")
m <- get_individual_df(myoclonic_seizure) %>% mutate(time="Myoclonic Seizure")
g <- get_individual_df(gen_seizure) %>% mutate(time="Generalized Seizure")

total <- bind_rows(p, m, g)

# df with the 4 removed bgs still in it
total_4in <- total

# Remove the 4 bgs
total <- total %>% dplyr::select(!c("53", "54", "57", "59"))
pca_data <- dplyr::select(total, 3:65) %>% as.matrix()

# PCA on the numeric data
time_pca <- prcomp(pca_data)
# Get summary fo PCA
#summary(time_pca)


# Set seed for permutation based significance
set.seed(1)
# Get standardized loadings using syndromics package
original_loadings <- stand_loadings(time_pca, pca_data)

# Set seed again for premutation based significance
set.seed(1)
# Do permutation based significance test on PCA data for Variance Accounted FOr
per <- permut_pc_test(time_pca, pca_data, P=10000, ndim=63, statistic="VAF")
per_results <- per$results

#write_csv(per_results, "d1_63_pcs_permutation_results.csv")
#per_results <- read_csv("d1_63_pcs_permutation_results.csv")
rownames(per_results) <- colnames(time_pca$rotation)

# Variance accounted for plot
VAF_plot(time_pca, pca_data, ndim=1:5, resample_ci=per_results)


set.seed(1)
# Permuation test for standardized loadings for PC1 and PC2
s_per<-permut_pc_test(time_pca, pca_data, P=1000, ndim=2, statistic = 's.loadings')

plot(s_per, plot_resample = T) +
   theme(axis.text.x=element_text(angle=90))

#write_csv(s_per$results, "63_permutation_laoding.csv")
s_per$results$Variables <- factor(s_per$results$Variables)

ggplot(s_per$results, aes(y=Variables)) +
  geom_col(aes(x=original, fill=original)) +
  geom_errorbar(aes(x=mean, xmin=ci_low, xmax=ci_high)) +
  facet_wrap(vars(component)) +
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  theme_minimal() +
  labs(fill="Original PC Loading")

splot<-syndromic_plot(time_pca, pca_data, ndim = 2,cutoff = 0.45, text_size = 5)

# Supplementary Figure 2A and 2B
splot$PC1
splot$PC2

loadings <- time_pca$rotation
loadings_long <- loadings %>% data.frame() %>% rownames_to_column(var="Group") %>% pivot_longer(cols=2:63)

loadings_long$positive <- loadings_long$value>0
pc_loadings <- loadings_long %>% filter(name=="PC1" | name=="PC2" | name=="PC3") %>% ggplot(aes(x=value, y=as.factor(as.numeric(Group)))) +
  geom_col(aes(fill=positive)) +
  labs(x="Loading", y="Behavior Group") +
  theme_clean() + 
  scale_fill_colorblind() +
  theme(legend.position="none") +
  facet_wrap(vars(name))

pc_loadings

#write_csv(loadings_long, "d1_64bg_loadings_long.csv")
#write_csv(pcs, "d1_64_pcpredicts.csv")

pcs <- predict(time_pca, dplyr::select(total, (3:65))) %>% as.data.frame()
pcs$time <- total$time
pcs$strain <- total$strain
pcs$id <- total$id

pcs$time <- factor(pcs$time, levels=c("Preictal", "Myoclonic Seizure", "Generalized Seizure"))


# Figure 2B
pca_scatter <- ggplot(pcs, aes(x=PC1, y=PC2, color=time)) +
  geom_point() +
  labs(x="PC1: 13.01%", y="PC2: 9.097%", title="PCA for All Mice", color="Seizure State", fill="Seizure State") +
  theme_minimal() +
  stat_ellipse(aes(fill=time), linetype=2, alpha=0.15, geom="polygon", level=0.95) +
  theme(panel.border = element_rect(fill=NA, linewidth=1), axis.ticks = element_line(linewidth=1, ), panel.grid=element_line(linewidth=1),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text=element_text(size=12, family="Arial")) +
  scale_color_manual(values=c("#A0AF84", "#F0746E", "#A52794")) +
  scale_fill_manual(values=c("#A0AF84", "#F0746E", "#A52794"))

pca_scatter
#ggsave("d1_4rm_pca_scatter.tiff",  dpi=300)

# Figure 2C
ggplot(pcs, aes(x=time, y=PC1)) +
  geom_violin(aes(fill=time), linewidth=1.5, draw_quantiles = c(.5)) +
  labs(x="Seizure State", y="Score on PC1", title="PC1 Scores Across Seizure States") +
  geom_jitter(width=.1, alpha=0.3) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=, binwidth=1/60, alpha=.8) +
  theme_classic() +
  scale_fill_manual(values=c("#A0AF84", "#F0746E", "#A52794")) +
  theme(legend.position="none", plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text=element_text(size=12, family="Arial"))

#ggsave("d1_4rm_pc1_violin.tiff", dpi=300)


# Figure 2D
ggplot(pcs, aes(x=time, y=PC2)) +
  geom_violin(aes(fill=time), linewidth=1.5, draw_quantiles=c(0.5)) +
  labs(x="Seizure State", y="Score on PC2", title="PC2 Scores Across Seizure States") +
  geom_jitter(width=.1, alpha=0.3) +
  theme_classic() +
  scale_fill_manual(values=c("#A0AF84", "#F0746E", "#A52794")) +
  theme(legend.position="none", plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        axis.text=element_text(size=12, family="Arial"))

#ggsave("d1_4rm_pc2_violin.tiff", dpi=300)

# Test normality
pcs %>% group_by(time) %>% shapiro_test(PC1)



# Statistical test for PC1
f <- friedman_test(pcs, PC1~time|id)
z <- wilcox_test(pcs, PC1~time, paired=T, p.adjust.method = "bonferroni") %>% as.data.frame()

#write_csv(f, "d1_63bg_friedman_test_pc1.csv")
#write_csv(z, "d1_63bg_wilcox_test_pc1.csv")

# Statistical Test for PC2
f2 <- friedman_test(pcs, PC2~time|id)
z2 <- wilcox_test(pcs, PC2~time, paired=T, p.adjust.method = "bonferroni") %>% as.data.frame()

#write_csv(f2, "d1_63bgs_friedman_pc2.csv")
#write_csv(z2, "d1_63bgs_wilcox_pc2.csv")
