#libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(RColorBrewer)
library(ggsignif)
library(ggpubr)
library(tidyverse)

#data
metabolites_df <- read.csv("/Users/renata.m/Desktop/Thesis/Rashidi_Data_2022/normalized_samples_matchsampleday.csv")  # Adjust the filename and path as needed

metabolites_df <- metabolites_df %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# List of metabolites of interest
metabolite_columns <- c('X1628', 'X1629', 'X1648', 'X136','X342')



# Renaming columns
metabolites_df <- metabolites_df %>%
  rename('glycochenodeoxycholate' = 'X1628', 'taurochenodeoxycholate' = 'X1629', 'taurocholate' = 'X1648', 
         'cholate' = 'X136', 'glycocholate' = 'X342')

#columns with new names
metabolite_columns_new <- c('glycochenodeoxycholate', 'taurochenodeoxycholate', 'taurocholate',
                            'cholate', 'glycocholate')

#long format
metabolites_long <- metabolites_df %>%
  pivot_longer(cols = all_of(metabolite_columns_new), names_to = "Metabolite", values_to = "Level")

#color palette
palette <- brewer.pal(3, "Paired")  # "Paired" palette with 3 colors

#boxplot
p <- ggplot(metabolites_long, aes(x = CDiff_infection, y = Level, fill = CDiff_infection)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  facet_wrap(~ Metabolite, scales = "free_y") +
  labs(x = " ", y = "Normalized Metabolite Levels") +
  theme_minimal() +
  scale_fill_manual(values = palette, 
                    name = "Group",
                    labels = c("False" = "Non-CDI", "True" = "CDI")) +
  theme(plot.title = element_text(size = 12), 
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 12))

# Add significance p-value
p <- p + stat_compare_means(method = "t.test", label = "p.signif", 
                            comparisons = list(c("True", "False")),
                            size = 4)

# show plot
print(p)

#############
#spearman correlation alpha diversity with metabolites 
shannon_entropy <- read_tsv('/Users/renata.m/QIIME2/qiime2-final_analysis-silva/diversity_metrics/exported-shannon-vector/alpha-diversity.tsv')

#join df
combined_df <- merge(metabolites_df, shannon_entropy, by = "ID")

#Spearman correlation
cor_test_results <- sapply(metabolite_columns_new, function(met) {
  test <- cor.test(combined_df$shannon_entropy, combined_df[[met]], method = "spearman", use = "complete.obs")
  c(correlation = test$estimate, p.value = test$p.value)
})
#df
correlation_df <- as.data.frame(t(cor_test_results))
rownames(correlation_df) <- metabolite_columns_new

correlation_df_sorted <- correlation_df[order(correlation_df$correlation, decreasing = TRUE), ]

correlation_df_sorted$significance <- cut(correlation_df_sorted$p.value,
                                          breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                          labels = c("***", "**", "*", "ns"))


#plot correlations
ggplot(correlation_df_sorted, aes(x = correlation_df_sorted$correlation, y = reorder(rownames(correlation_df_sorted), correlation_df_sorted$correlation))) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  geom_text(aes(label = significance), hjust = -0.5, color = "black", size = 5) +
  labs(x = "Spearman Correlation", y = "Metabolites", 
       title = " ") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) 
