#Metacoder
library(microViz)
library(metacoder)
library(cowplot)
library(ggplot2)


#BIPOLARIS
setwd = ("C:/Users/2022lj002/Documents/R/C1/propre/16S/namco bipo/metacoder") 
path <- "C:/Users/2022lj002/Documents/R/C1/propre/16S/namco bipo/metacoder"
list.files(path)
hmp_otus <- read.csv2("C:/Users/2022lj002/Documents/R/C1/propre/16S/namco bipo/metacoder/metacoder bipo.csv", sep=";", header=TRUE)
hmp_otus
hmp_samples <- read.csv2("C:/Users/2022lj002/Documents/R/C1/propre/16S/namco bipo/metacoder/metadata_all_16S_C1.csv", sep=";", header=TRUE)
hmp_samples



tm_obj_bipo = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
                             class_key = c(tax_rank = "taxon_rank", tax_name = "taxon_name"),
                             class_regex = "^(.+)__(.+)$")


# Convert counts to proportions
tm_obj_bipo$data$otu_table <- calc_obs_props(tm_obj_bipo, data = "tax_data", cols = hmp_samples$SampleID)
#tm_obj_bipo$data$otu_table <- tm_obj_bipo$data$tax_data

tm_obj_bipo$data$otu_table


# Get per-taxon counts
tm_obj_bipo$data$tax_table <- calc_taxon_abund(tm_obj_bipo, data = "otu_table", cols = hmp_samples$SampleID)
#tm_obj_bipo$data$tax_table <- tm_obj_bipo$data$otu_table
tm_obj_bipo

hmp_samples

unique(hmp_samples$Condition)

comparisons <- list(c("Disease","Healthy"))

# Calculate difference between groups
tm_obj_bipo$data$diff_table <- compare_groups(tm_obj_bipo, data = "tax_table",
                                              cols = hmp_samples$SampleID,
                                              groups = hmp_samples$Condition,
                                              combinations = comparisons)
print(tm_obj_bipo$data$diff_table)

## Make comparison plots

plot_comp <- function(comp_pair) {
  set.seed(1) # so plot layout and sizes always look the same
  tm_obj_bipo %>%
    filter_obs(data = "diff_table", treatment_1 %in% comp_pair & treatment_2 %in% comp_pair) %>%
    heat_tree(node_size = n_obs,
              node_size_range = c(0.01, 0.05),
              node_color = log2_median_ratio,
              node_label = taxon_names,
              node_color_range = diverging_palette(),
              node_color_trans = "linear",
              node_color_interval = c(-3, 3),
              edge_color_interval = c(-3, 3),
              node_size_axis_label = "Number of ASV",
              node_color_axis_label = "Log2 ratio median proportions",
              title = paste0(comp_pair[1], ' vs. ', comp_pair[2]),
              make_node_legend = TRUE, 
              make_edge_legend = TRUE)
}


comp_plots <- lapply(comparisons, plot_comp)
comp_plots


set.seed(1) # so plot layout and sizes always look the same
key_plot <- tm_obj_bipo %>%
  heat_tree(node_size = n_obs,
            node_size_range = c(0.01, 0.05),
            node_label = taxon_names,
            node_color = "grey",
            node_color_range = diverging_palette(),
            node_color_trans = "linear",
            node_color_interval = c(-5, 5),
            edge_color_interval = c(-5, 5),
            node_size_axis_label = "Number of OTUs",
            node_color_axis_label = "Log2 ratio median proportions")

combined_plot <- plot_grid(plot_grid(plotlist = comp_plots), key_plot, ncol = 1)
combined_plot




#pour enlever les labels NA au niveau de l'espèce
tm_obj_bipo %>% 
  filter_taxa(taxon_ranks == "s", supertaxa = TRUE) %>% # subset to the species rank
  filter_taxa(taxon_ranks != "s" | !grepl("NA", taxon_names)) %>% 
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,
            node_size_range = c(0.01, 0.05),
            node_color = log2_median_ratio,
            node_color_axis_label = "Log2 ratio median proportions",
            node_color_interval = c(-3, 3),
            edge_color_interval = c(-3, 3),
            node_color_range = (c("#2EC4B6",'#E5DCC5','#FC6471')),
            layout = "davidson-harel", initial_layout = "reingold-tilford")

























# Load required libraries
library(microViz)
library(metacoder)
library(cowplot)
library(ggplot2)

# Set working directory and define path
setwd("path/to/directory")
path <- "path/to/directory"
list.files(path)

# Import necessary data files
hmp_otus <- read.csv2("path/to/directory/metacoder_data.csv", sep=";", header=TRUE)
hmp_samples <- read.csv2("path/to/directory/metadata.csv", sep=";", header=TRUE)

# Parse taxonomic data
tm_obj <- parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
                         class_key = c(tax_rank = "taxon_rank", tax_name = "taxon_name"),
                         class_regex = "^(.+)__(.+)$")

# Convert counts to proportions
tm_obj$data$otu_table <- calc_obs_props(tm_obj, data = "tax_data", cols = hmp_samples$SampleID)

# Get per-taxon counts
tm_obj$data$tax_table <- calc_taxon_abund(tm_obj, data = "otu_table", cols = hmp_samples$SampleID)

# Define condition comparisons
comparisons <- list(c("Disease","Healthy"))

# Calculate difference between groups
tm_obj$data$diff_table <- compare_groups(tm_obj, data = "tax_table",
                                         cols = hmp_samples$SampleID,
                                         groups = hmp_samples$Condition,
                                         combinations = comparisons)

# Function to generate comparison plots
plot_comp <- function(comp_pair) {
  set.seed(1)
  tm_obj %>%
    filter_obs(data = "diff_table", treatment_1 %in% comp_pair & treatment_2 %in% comp_pair) %>%
    heat_tree(node_size = n_obs,
              node_size_range = c(0.01, 0.05),
              node_color = log2_median_ratio,
              node_label = taxon_names,
              node_color_range = diverging_palette(),
              node_color_trans = "linear",
              node_color_interval = c(-3, 3),
              edge_color_interval = c(-3, 3),
              node_size_axis_label = "Number of ASVs",
              node_color_axis_label = "Log2 ratio median proportions",
              title = paste0(comp_pair[1], ' vs. ', comp_pair[2]),
              make_node_legend = TRUE, 
              make_edge_legend = TRUE)
}

# Generate plots
comp_plots <- lapply(comparisons, plot_comp)

# Generate key reference plot
set.seed(1)
key_plot <- tm_obj %>%
  heat_tree(node_size = n_obs,
            node_size_range = c(0.01, 0.05),
            node_label = taxon_names,
            node_color = "grey",
            node_color_range = diverging_palette(),
            node_color_trans = "linear",
            node_color_interval = c(-5, 5),
            edge_color_interval = c(-5, 5),
            node_size_axis_label = "Number of OTUs",
            node_color_axis_label = "Log2 ratio median proportions")

# Combine comparison plots with reference plot
combined_plot <- plot_grid(plot_grid(plotlist = comp_plots), key_plot, ncol = 1)
print(combined_plot)

# Remove NA labels at the species level
tm_obj %>% 
  filter_taxa(taxon_ranks == "s", supertaxa = TRUE) %>%
  filter_taxa(taxon_ranks != "s" | !grepl("NA", taxon_names)) %>% 
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,
            node_size_range = c(0.01, 0.05),
            node_color = log2_median_ratio,
            node_color_axis_label = "Log2 ratio median proportions",
            node_color_interval = c(-3, 3),
            edge_color_interval = c(-3, 3),
            node_color_range = c("#2EC4B6",'#E5DCC5','#FC6471'),
            layout = "davidson-harel", initial_layout = "reingold-tilford")

