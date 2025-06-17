# 
# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
# Behavior ID of Mouse Strains

# Load dependencies
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(dendextend)
library(pvclust)

# Get full dataframe
df <- read_csv("mouse_data.csv")


# Separate into time frames
preictal <- df %>% filter(scorer<=3600)
myoclonic_seizure <- df %>% filter(scorer>3600 & scorer<=(3600+((GST-MST)*30)))
gen_seizure <- df %>% filter(scorer>(3600+((GST-MST)*30)))


# Define Functions

# get_strain_df
# Inputs
#   df - A Data Frame with columns strain, and BSOiD_labels of individual mice with behavioral group classification
#       by frame
# Outputs
#   nums_long - A long Data Frame with the percentage of behavioral group by Mouse Strain for the BSOiD labels

get_strain_df <- function(df) {
  p <- df %>% group_by(strain) %>% count(BSOiD_labels)
  p2 <- p %>% pivot_wider(id_cols=strain, names_from=BSOiD_labels, values_from=n)
  p2[is.na(p2)] <- 0
  
  nums_wide <- p2[,-1]/rowSums(p2[,-1])
  nums_wide$strain <- p2$strain
  nums_wide <- relocate(nums_wide, "strain")
  return(nums_wide)
}

# get_heatmap
# Inputs
#   long_df - A wide
#   title - The title for the plot, default is NULL
# Outputs
#   A pheatmap object as defined in pheatmap package

get_heatmap <- function(nums_wide, annotation_df, color_list, title=NULL) {
  # Get numeric values
  nums_wide_numeric <- nums_wide[,-1]
  # Scale numeric values by behavior group
  nums_wide_numeric <- scale(nums_wide_numeric)
  # Set rownames to strain
  rownames(nums_wide_numeric) <- nums_wide[,1]
  mat <- as.matrix(nums_wide_numeric)
  
  return(pheatmap(t(mat), scale="none", annotation_col=annotation_df, annotation_colors=color_list, xlab="Mouse Strain", ylab="Behavior Group", main=title,
                  fontsize=8))
}

# get dfs for seizure states
p <- get_strain_df(preictal) %>% dplyr::select(!c("53", "54", "57", "59"))
m <- get_strain_df(myoclonic_seizure) %>% dplyr::select(!c("53", "54", "57", "59"))
g <- get_strain_df(gen_seizure) %>% dplyr::select(!c("53", "54", "57", "59"))


strs <- p$strain
strs[1] <- "B6J"

p$strain <- strs
m$strain <- strs
g$strain <- strs

# Get color and strain combos
p_cols <- preictal %>% dplyr::select(strain, Color) %>% unique()
p_cols <- p_cols %>% mutate(Color=ifelse(strain=="C57", "black", Color)) %>% drop_na()
p_cols <- p_cols %>% mutate(strain=ifelse(strain=="C57", "B6J", strain))

# Set Color for Tiles
c_list <- list(Color=c(agouti="#b08968", black="black", white="white"))
color_fix_df <- mutate(p_cols, Color=case_when(str_detect(Color, "agouti") ~ "agouti", .default=as.character(Color)))
color_df <- color_fix_df %>% drop_na()
color_df <- column_to_rownames(color_df, var="strain")



# Figure 3A
h1 <- get_heatmap(p, color_df, c_list, "Pre-Ictal")
#ggsave("d1_63_p_dendrogramheatmap.tiff",  dpi=300)

# Figure 3B
h2 <- get_heatmap(m, color_df, c_list, "Myoclonic Seizure")
#ggsave("d1_63_m_dendrogramheatmap.tiff",  dpi=300)

# Figure 3C
h3 <- get_heatmap(g, color_df, c_list, "Generalized Seizure")
#ggsave("d1_63_g_dendrogramheatmap.tiff",  dpi=300)

# Get scaled values for P, M, and GS to get dendrograms with p-values
p_nwn <- p[,-1]
p_nwn <- scale(p_nwn)
rownames(p_nwn) <- strs
p_mat <- as.matrix(p_nwn)

m_nwn <- m[,-1]
m_nwn <- scale(m_nwn)
rownames(m_nwn) <- strs
m_mat <- as.matrix(m_nwn)

g_nwn <- g[,-1]
g_nwn <- scale(g_nwn)
rownames(g_nwn) <- strs
g_mat <- as.matrix(g_nwn)

# Do bootsrap hierarchical clustering with complete linkage and euclidean distance
p_clust <- pvclust(t(p_mat), method.hclust="complete", method.dist="euclidean", iseed=2024)
m_clust <- pvclust(t(m_mat), method.hclust="complete", method.dist="euclidean", iseed=2024)
g_clust <- pvclust(t(g_mat), method.hclust="complete", method.dist="euclidean", iseed=2024)

# Get dfs with p-values and standard errors
p_edge <- p_clust$edges
m_edge <- m_clust$edges
g_edge <- g_clust$edges

#write_csv(p_edge, "d1_preictal_clustering_values.csv")
#write_csv(m_edge, "d1_ms_clustering_values.csv")
#write_csv(g_edge, "d1_gs_clustering_values.csv")

# Supplementary Figure 4A
plot(p_clust)
pvrect(p_clust)

# Supplementary Figure 4B
plot(m_clust)
pvrect(m_clust)

# Supplementary Figure 4C
plot(g_clust)
pvrect(g_clust)

#write_csv(rownames_to_column(as.data.frame(t(p_mat))), "preictal_standized_usage_stain.csv")
#write_csv(rownames_to_column(as.data.frame(t(m_mat))), "ms_standized_usage_stain.csv")
#write_csv(rownames_to_column(as.data.frame(t(g_mat))), "gs_standized_usage_stain.csv")



nums_wide_numeric <- p[,-1]
# Scale numeric values by behavior group
nums_wide_numeric <- scale(nums_wide_numeric) %>% as.data.frame()
nums_wide_numeric$strain <- p$strain
p_sum <- nums_wide_numeric %>% mutate(color=color_df$Color)
p_sum %>% group_by(color) %>% summarize(bg_66=mean(`66`))
