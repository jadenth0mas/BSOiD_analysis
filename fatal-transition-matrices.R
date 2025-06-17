# R Code for Manuscript as Part of Gu Lab
# 
# Behavior decoding delineates seizure microfeatures and associated sudden death risks in mouse models of epilepsy
# 
# Behavior biomarker of SUDEP 

library(tidyverse)
library(ggsci)
library(reshape2)
library(ggthemes)
library(viridis)
library(circlize)
library(ggrepel)

# Load data
df <- read_csv("mouse_data.csv")
df_fatal <- df %>% mutate(fatal=(Score==7))

preictal <- df_fatal %>% filter(scorer<=3600) %>% mutate(time="preictal")
myoclonic_seizure <- df_fatal %>% filter(scorer>3600 & scorer<=(3600+((GST-MST)*30))) %>% mutate(time="myoclonic seizure")
gen_seizure <- df_fatal %>% filter(scorer>(3600+((GST-MST)*30))) %>% mutate(time="generalized seizure")

# Set removed BGs to NA
preictal[preictal$BSOiD_labels=="53" | preictal$BSOiD_labels=="54" | preictal$BSOiD_labels== "57" | preictal$BSOiD_labels== "59",]$BSOiD_labels <- NA
myoclonic_seizure[myoclonic_seizure$BSOiD_labels=="53" | myoclonic_seizure$BSOiD_labels=="54" | myoclonic_seizure$BSOiD_labels== "57" | myoclonic_seizure$BSOiD_labels== "59",]$BSOiD_labels <- NA
gen_seizure[gen_seizure$BSOiD_labels=="53" | gen_seizure$BSOiD_labels=="54" | gen_seizure$BSOiD_labels== "57" | gen_seizure$BSOiD_labels== "59",]$BSOiD_labels <- NA

# Get the next BG within each id
gen_seizure$nextGroup <- ifelse(lead(gen_seizure$id)==gen_seizure$id, lead(gen_seizure$BSOiD_labels), NA)

# Filter by fatality
g_fatal <- filter(gen_seizure, fatal==T)
g_nonfatal <- filter(gen_seizure, fatal==F)

# Get fatal transition matrix
f_tm_pre <- as.matrix(table(g_fatal$BSOiD_labels, g_fatal$nextGroup))
f_tm <- f_tm_pre

# Remove diagonal to look at other transitions
diag(f_tm) <- NA

# Get sum in each row 
row_sums <- rowSums(f_tm, na.rm = TRUE)

# Get probability matrix, dividing each row by the row sum
prob_matrix <- sweep(f_tm, 1, row_sums, FUN = "/")

# Melt matrix into 2 columns
melted_fatal_prob <- melt(prob_matrix, na.rm=T)

# Figure 4C
# Fatal
ggplot(melt(prob_matrix, na.rm=T), aes(x=as.factor(Var1), y=as.factor(Var2), fill=value)) +
  geom_tile() +
  labs(fill="Probability", title="Fatal Transition Probability Matrix", x="Behavior Group", y="Behavior Group") +
  scale_fill_viridis(option="F") +
  #geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
  theme_clean() +
  theme(axis.text.x=element_text(angle=90, size=20, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=20, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        panel.grid.major.y = element_blank()) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2), expand = expansion(add = 1)) +
  scale_y_discrete(guide = guide_axis(n.dodge = 2), expand = expansion(add = 1))


# Get nonfatal transition matrix
nf_tm_pre <- as.matrix(table(g_nonfatal$BSOiD_labels, g_nonfatal$nextGroup))
nf_tm <- nf_tm_pre

diag(nf_tm) <- NA


# Get rowsums and probability
row_sums_nf <- rowSums(nf_tm, na.rm = TRUE)
prob_matrix_nf <- sweep(nf_tm, 1, row_sums_nf, FUN = "/")

# Melt matrix
melted_nonfatal_prob <- melt(prob_matrix_nf, na.rm=T)

# Figure 4C
# Nonfatal
ggplot(melted_nonfatal_prob, aes(x=as.factor(Var1), y=as.factor(Var2), fill=value)) +
  geom_tile() +
  labs(fill="Probability", title="Non-Fatal Transition Probability Matrix", x="Behavior Group", y="Behavior Group") +
  scale_fill_viridis(option="F", limits=c(0, 1)) +
  #geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
  theme_clean() +
  theme(axis.text.x=element_text(angle=90, size=20, vjust=.5, family="Arial"), panel.grid=element_blank(), axis.line = element_line(), axis.ticks=element_line(),
        axis.text.y=element_text(size=20, family="Arial"),
        plot.title=element_text(family="Arial", size=20),
        axis.title = element_text(size=14),
        panel.grid.major.y = element_blank()) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2), expand = expansion(add = 1)) +
  scale_y_discrete(guide = guide_axis(n.dodge = 2), expand = expansion(add = 1))


# Transitions per mouse for F and NF

melted_tm <- melt(f_tm, na.rm=T)
melted_tm$value <- melted_tm$value/length(unique(g_fatal$id))

melted_nftm <- melt(nf_tm, na.rm=T)
melted_nftm$value <- melted_nftm$value/length(unique(g_nonfatal$id))

#write_csv(melted_nftm, "non-fatal-transition-probabilities-fixed.csv")
#write_csv(melted_tm, "fatal-transition-probabilities-fixed.csv")


# Get total number of transitions for each
total_fatal <- sum(f_tm, na.rm=T)
total_nonfatal <- sum(nf_tm, na.rm=T)

# Scale fatal transition matrix since nonfatal more obs
scaled_ftm <- (f_tm / total_fatal) * total_nonfatal

# Get integer scaled transitions
integer_ftm <- floor(scaled_ftm)
integer_nftm <- floor(nf_tm)

# Get BGs available for fatal in both row and column
fatal_rows <- rownames(integer_ftm)
fatal_cols <- colnames(integer_ftm)

# Subset nonfatal to only include the BG rows and cols in fatal
non_fatal_subset <- integer_nftm[fatal_rows, fatal_cols, drop = FALSE]

# Get integer total of scaled and subsetted data
total_fatal_transitions <- sum(integer_ftm, na.rm=T)
total_non_fatal_transitions <- sum(non_fatal_subset, na.rm=T)

# Label fatal and nonfatal counts to use
fatal_counts <- integer_ftm
non_fatal_counts <- non_fatal_subset

# Initialize chi-square result matrix
p_values <- matrix(NA, nrow = nrow(fatal_counts), ncol = ncol(fatal_counts),
                   dimnames = list(fatal_rows, fatal_cols))

conf_ints <- matrix(NA, nrow = nrow(fatal_counts), ncol = ncol(fatal_counts),
                   dimnames = list(fatal_rows, fatal_cols))

# Minimum threshold for transitions
min_transitions <- 5

set.seed(2024)

adj_level <- 1 - (0.05)/(382)

# Perform chi-square tests on matching cells
for (i in 1:nrow(fatal_counts)) {
  for (j in 1:ncol(fatal_counts)) {
    # If both fatal and nonfatal above threshold and not NA
    if (!is.na(fatal_counts[i, j]) && !is.na(non_fatal_counts[i,j])) {
      if (fatal_counts[i, j] >= min_transitions && non_fatal_counts[i, j] >= min_transitions) {
        # Create a 2x2 contingency table
        observed <- matrix(
          c(fatal_counts[i, j], total_fatal_transitions - fatal_counts[i, j],
            non_fatal_counts[i, j], total_non_fatal_transitions - non_fatal_counts[i, j]),
          nrow = 2, byrow = TRUE
        )
      
        # Chi-Square
        test <- chisq.test(observed, correct = TRUE)
        ci <- prop.test(x=observed[,1], n=observed[,2], conf.level=adj_level, alternative="two.sided", correct=T)
        # Put p-value in initialized matrix
        p_values[i, j] <- test$p.value
        conf_ints[i, j] <- paste(ci$conf.int[1], ci$conf.int[2], sep=", ")
      } else {
      # Else have p-value NA for below threshold counts
      p_values[i, j] <- NA
      conf_ints[i, j] <- NA
      }
    } else {
      # Else have p-value NA if count is NA
      p_values[i, j] <- NA
      conf_ints[i, j] <- NA
    }
    
  }
}

# Get results
chi_square_results <- p_values %>% as.data.frame() %>% rownames_to_column() %>% pivot_longer(cols=2:54) %>% filter(!is.na(value)) %>% mutate(adjusted_pvalue=p.adjust(value))
prop_test_results <- conf_ints %>% as.data.frame() %>% rownames_to_column() %>% pivot_longer(cols=2:54) %>% filter(!is.na(value)) %>% rename(ci=value)

chi_results <- chi_square_results %>% left_join(prop_test_results, by=join_by(rowname==rowname, name==name))

#write_csv(chi_square_results, "fatal_vs_nonfatal_chi_square.csv")

# Compute volcano plot
fatal_probs <- melt(fatal_counts)
non_fatal_probs <- melt(non_fatal_counts)

non_fatal_probs <- non_fatal_probs %>% rename("nf-usage"=value) %>% unite("transition", Var1, Var2, sep="->")
fatal_probs <- fatal_probs %>% rename("f-usage"=value) %>%  unite("transition", Var1, Var2, sep = "->")

total_probs <- left_join(fatal_probs, non_fatal_probs, by=join_by("transition"))

# Show direction of transition
chi_square_results <- chi_square_results %>% unite("transition", rowname, name, sep="->")


joined_p <- chi_square_results %>% left_join(total_probs)
joined_p <- joined_p %>% mutate(fold_change=`f-usage`/`nf-usage`, log2_fold_change = log2(fold_change),
                                significant = adjusted_pvalue < 0.05 & abs(log2_fold_change) > 1,
                                extra_significant = adjusted_pvalue < 0.01 & abs(abs(log2_fold_change)>2.5),
                                label = ifelse(extra_significant, transition, NA),
                                regulation=ifelse(significant, ifelse(log2_fold_change>0, "Fatal Up", "Fatal Down"), "Not Significant"))


# Figure 4D
ggplot(joined_p, aes(x=log2_fold_change, y=-log10(adjusted_pvalue), color=regulation)) +
  geom_point(size=4) +
  scale_color_manual(values=c("#0072B5FF", "#BC3C29FF", "gray")) +
  geom_point(data=filter(joined_p, extra_significant), aes(color=regulation)) +
  labs(x="Log2 Fold Change", y="-Log10 Adjusted P-Value") +
  #geom_text_repel(aes(label = label), color="black", size=12, max.overlaps = 10,      # Controls the max number of overlapping labels shown
  #               box.padding = 0.3,      # Padding around each label
  #                point.padding = 0.4,    # Distance between point and label
  #                nudge_y = 0.2, na.rm = TRUE) +
  theme_minimal() +
  theme(legend.position = "none", plot.title=element_text(family="Arial", size=15), 
        panel.border = element_rect(fill=NA, linewidth=1), axis.ticks = element_line(linewidth=1), 
        panel.grid=element_blank(), 
        axis.title.y = element_text(size=30),
        axis.text=element_text(size=12, family="Arial"),
        axis.title.x=element_text(size=30))

#write_csv(joined_p, "directional_fatal_stats_and_adjusted_counts.csv")
