library(tidyverse)
library(ggplot2)
library(ggtext)
library(RColorBrewer)

#feature table
Tab <- read_tsv("/Users/renata.m/QIIME2/qiime2-final_analysis-silva_abx/feature-table.tsv",
                col_types = cols(otu_id = col_character(),
                                 .default = col_number()))
Tab


#taxonomy table 
dat <- Tab %>%
  pivot_longer(-otu_id, names_to = "ID", values_to = "count")
dat


#taxonomy and feature table merged
Tax <- read_tsv("/Users/renata.m/QIIME2/qiime2-final_analysis-silva_abx/Taxonomy/taxonomy/taxonomy.tsv", 
                col_types = cols(.default = "character")) %>%
  select("Feature ID", "Taxon") %>%
  rename_all(tolower) %>%
  mutate(taxon = str_replace_all(taxon, "\\(\\d\\)", ""),
         taxon = str_replace(taxon, ":$", "")) %>% 
  separate(taxon, 
           into= c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";")

colnames(Tax)[colnames(Tax) == 'feature id'] <- 'otu_id'
Tax


#metadata
metadata <- read_tsv("/Users/renata.m/QIIME2/qiime2-final_analysis-silva_abx/fondue-output/metadata_abx_before_after_5d.tsv")
metadata


#otu count
otu_rel_abund <- inner_join(metadata, dat, by="ID") %>%
  inner_join(., Tax, by="otu_id") %>%
  group_by(ID) %>%
  mutate(rel_abund = count/sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "otu_id"),
               names_to = "level", 
               values_to = "taxon")
otu_rel_abund



#relative abundnace of taxon
taxon_rel_abund <- otu_rel_abund %>%
  filter(level=="genus") %>%
  group_by(CDiff_infection, ID,SampleDayRelativeToD1Chemo, taxon) %>%
  summarise(rel_abund=sum(rel_abund),.groups = "drop") %>%
  group_by(CDiff_infection,SampleDayRelativeToD1Chemo, taxon) %>%
  summarise(mean_rel_abund= 100*mean(rel_abund), .groups = "drop") %>% 
  mutate(taxon=str_replace(taxon,
                           "^(\\S*)$", "*\\1*"))
taxon_rel_abund

taxon_pool <- taxon_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) <4, 
            mean = mean(mean_rel_abund), 
            .groupd="drop")
print(taxon_pool)

# List of taxa to plot
taxa_to_plot <- c(
  ' g__Actinomyces',
  ' g__Lactococcus',
  ' g__Subdoligranulum', 
  ' g__Barnesiella',
  ' g__Coprococcus',
  ' g__Dorea',
  ' g__Faecalibacterium',
  ' g__Dialister',
  ' g__Roseburia',
  ' g__Hungatella',
  ' g__Prevotella',
  ' g__Enterococcus'
)

# make specific genera next to each other
taxon_rel_abund$taxon <- factor(taxon_rel_abund$taxon, levels = c(
  ' g__Actinomyces',
  ' g__Lactococcus',
  ' g__Barnesiella',
  ' g__Dorea',
  ' g__Subdoligranulum', 
  ' g__Coprococcus',
  ' g__Faecalibacterium',
  ' g__Enterococcus',
  ' g__Dialister',
  ' g__Roseburia',
  ' g__Hungatella',
  ' g__Prevotella'
))

# Filter
taxa_to_plot <- taxon_rel_abund %>%
  filter(taxon %in% taxa_to_plot)

#plot each taxon
plot <- ggplot(taxa_to_plot, aes(x = SampleDayRelativeToD1Chemo, y = mean_rel_abund, color = CDiff_infection)) + 
  geom_line() + 
  facet_wrap(~taxon, scales = "free_y",
             ncol = 4) + 
  theme_classic() + 
  scale_color_manual(values = palette,
                     name = "Group",
                     labels = c("Non-CDI", "CDI")) +
  labs(x = "Time (Days Since Experiment Start)", y = "Relative Abundance (%)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10))

# Display the plot
print(plot)
