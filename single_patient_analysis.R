#libraries
library(tidyverse)
library(ggplot2)
library(ggtext)
library(RColorBrewer)

#data
feature_table <- read_tsv("/Users/renata.m/QIIME2/qiime2-final_analysis-silva_abx/feature-table.tsv",
                          col_types = cols(otu_id = col_character(),
                                           .default = col_number()))

# Drop the sample (column) duplicate 
feature_table <- feature_table %>% select(-SRR7048601)

taxonomy <- read_tsv("/Users/renata.m/QIIME2/qiime2-final_analysis-silva_abx/Taxonomy/taxonomy/taxonomy.tsv", col_types = cols(.default = "character")) %>%
  select("Feature ID", "Taxon") %>%
  rename_all(tolower) %>%
  mutate(taxon = str_replace_all(taxon, "\\(\\d\\)", ""),
         taxon = str_replace(taxon, ":$", "")) %>%
  separate(taxon, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";")

# Drop the sample (column) duplicate 
colnames(taxonomy)[colnames(taxonomy) == 'feature id'] <- 'otu_id'

metadata <- read_tsv("/Users/renata.m/QIIME2/qiime2-final_analysis-silva_abx/fondue-output/shannon_sra_metadata.tsv",
                     col_types = cols(Patient_ID = col_character())) 

metadata <- metadata[metadata$ID != 'SRR7048601', ]

antibiotics <- read_csv("/Users/renata.m/Desktop/Thesis/Rashidi_Data_2022/AntibioticsTable.csv", skip = 3,
                        col_types = cols(Patient_ID = col_character(), AntibacterialABx = col_character(), 
                                         Route = col_character(),.default = col_number()))
# Filter data for patient 7D010
metadata_patient <- metadata %>% filter(Patient_ID == "7D010")
antibiotics_patient <- antibiotics %>% filter(Patient_ID == "7D010")


# Merge feature table with metadata by ID
feature_metadata <- feature_table %>%
  pivot_longer(-otu_id, names_to = "ID", values_to = "count") %>%
  inner_join(metadata_patient %>% select(ID, Patient_ID), by = "ID")


##Relative Abundance ##

# Calculate relative abundance for patient 7D010
otu_rel_abund <- inner_join(feature_metadata,metadata_patient, by = "ID") %>%
  inner_join(.,taxonomy, by="otu_id") %>%
  group_by(ID) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "otu_id"),
               names_to = "level", values_to = "taxon")

taxon_rel_abund <- otu_rel_abund %>%
  filter(level=="genus") %>%
  group_by(SampleDayRelativeToD1Chemo, ID, taxon) %>%
  summarise(rel_abund=sum(rel_abund),.groups = "drop") %>%
  group_by(SampleDayRelativeToD1Chemo, taxon) %>%
  summarise(mean_rel_abund= 100*mean(rel_abund), .groups = "drop")%>%
  ungroup()


taxon_rel_abund

#pool low abundance data
taxon_pool <- taxon_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) <6, mean = mean(mean_rel_abund),
            .groupd="drop")

taxon_pool

taxon_rel_abund <- inner_join(taxon_rel_abund, taxon_pool, by="taxon")%>% 
  mutate(taxon = if_else(pool, "Other", taxon))%>%
  group_by(SampleDayRelativeToD1Chemo, taxon)%>%
  summarise(mean_rel_abund = sum(mean_rel_abund),
            .groups = "drop") %>%
  mutate(
    taxon = if_else(taxon == "Other", "Other", str_replace(taxon, "(\\S*)$", "*\\1*")))
 

#antibiotic treatment
antibiotics_timeline <- ggplot(antibiotics_patient, aes(x = StartDayRelativeToD1Chemo, y = AntibacterialABx, color = Route)) +
  geom_segment(aes(xend = StopDayRelativeToD1Chemo, yend = AntibacterialABx), size = 5) +
  scale_color_manual(values = c("IV" = "skyblue", "PO" = "orange")) +
  labs(x = NULL, y = "Antibiotics") +
  theme_minimal()+
  scale_x_continuous(limits = c(-5, 30))+
  theme(legend.position = "right")

antibiotics_timeline


#shannon entropy
#data frame for vertical lines (for clincial events)
vline_data <- data.frame(
  Event = c("First C. diff Day", "Diarrhea Start Day", "Diarrhea End Day", "First NF Day"),
  xintercept = c(metadata_patient$FirstCDiffDayRelativeToD1Chemo[1],
                 metadata_patient$DiarrheaStartDayRelativeToD1Chemo[1]+ 0.1,
                 metadata_patient$DiarrheaEndDayRelativeToD1Chemo[1],
                 metadata_patient$FirstNFDayRelativeToD1Chemo[1])
)

# Corrected Shannon Entropy plot with a legend for the vertical lines
shannon_entropy_plot <- ggplot(metadata_patient, aes(x = SampleDayRelativeToD1Chemo, y = shannon_entropy)) +
  geom_line() +
  geom_point() +
  geom_vline(data = vline_data, aes(xintercept = xintercept, color = Event), 
             linetype = "dashed", show.legend = TRUE) +
  scale_color_manual(values = c("First C. diff Day" = "red", 
                                "Diarrhea Start Day" = "green", 
                                "Diarrhea End Day" = "purple", 
                                "First NF Day" = "blue"),
                     labels = c("red" = "First C. diff Day", 
                                "green" = "Diarrhea Start Day", 
                                "purple" = "Diarrhea End Day", 
                                "blue" = "First NF Day"),
                     name = "Event") +
  labs(x = NULL, y = "Shannon Entropy") +
  theme_minimal() +
  scale_x_continuous(limits = c(-5,30)) +
  theme(legend.position = "right")

shannon_entropy_plot

#
nb.cols <- length(unique(taxon_rel_abund$taxon))
mycolors <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols)

stacked_bar_plot <- ggplot(taxon_rel_abund, aes(x = SampleDayRelativeToD1Chemo, y = mean_rel_abund, fill = taxon)) +
  geom_col() +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(x = "Days Relative to D1 Chemo", y = "Relative Abundance (%)", fill = "Taxon") +
  theme_minimal() +
  theme(legend.text = element_markdown(), 
        legend.key.size = unit(10, "pt"), 
        legend.position = "right")+
  scale_x_continuous(limits = c(-5,30))

stacked_bar_plot


#combine all plots
#install.packages("patchwork")
library(patchwork)
combined_plot <- antibiotics_timeline / shannon_entropy_plot / stacked_bar_plot + plot_layout(ncol = 1)
combined_plot

ggsave("combined_plot.tiff", combined_plot, width = 20, height = 10)




