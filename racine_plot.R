# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
# 
# Behavior ID of Mouse Strains

# Load dependencies
library(tidyverse)
library(ggsci)
library(ggthemes)
library(rstatix)

# Load data
df <- read_csv("mouse_data.csv")
df_fatal <- df %>% mutate(fatal=(Score==7))


# Partition based on seizure state
preictal <- df_fatal %>% filter(scorer<=3600) %>% mutate(time="preictal")
myoclonic_seizure <- df_fatal %>% filter(scorer>3600 & scorer<=(3600+((GST-MST)*30))) %>% mutate(time="myoclonic seizure")
gen_seizure <- df_fatal %>% filter(scorer>(3600+((GST-MST)*30))) %>% mutate(time="generalized seizure")

# Get strain with median racine
summed_df <- df %>% group_by(strain) %>% summarize(median_racine=median(Score), sd_racine=sd(Score)) %>% ungroup()

# Join df with Racine score
joined_df <- left_join(df, summed_df)
racine_scores <- df_fatal %>% dplyr::select(strain, id, Score) %>% unique()

#racine_kruskal <- kruskal_test(Score~strain, data=racine_scores)


# Mann U Whitney Test
racine_pairwise_mannu <- pairwise.wilcox.test(racine_scores$Score, racine_scores$strain, p.adjust.method = "bonferroni")

#write_csv(racine_kruskal, "racine-kruskal.csv")
#write_csv(data.frame(racine_pairwise_mannu$p.value), "racine-pairwise-mann-u-whitney.csv")

joined_df <- joined_df %>% select(strain, id, Score, median_racine, sd_racine) %>% unique()
joined_df <- joined_df %>% mutate(strain=ifelse(strain=="C57", "B6J", strain))


# Figure 3D
ggplot(joined_df, aes(x=reorder(strain, median_racine, decreasing=T), y=Score, group=strain)) +
  geom_jitter(width=0.2, height=0.2, alpha=.7) +
  geom_crossbar(color="black", aes(ymin=median_racine, ymax=median_racine, y=median_racine), width=.75) +
  labs(x="Strain", y="Racine Score") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size=12, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=12, family="Arial"),
        axis.title = element_text(size=14), legend.position="none") +
  scale_y_continuous(limits=c(1, 8), breaks=seq(1, 7, by=1))


#write_csv(joined_df, "d1_racine_data.csv")
