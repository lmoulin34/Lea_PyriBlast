# Load required libraries
library(ggplot2)
library(ggpubr)

# Set working directory and define path
setwd("path/to/directory")
path <- "path/to/directory"
list.files(path)

# Import dataset
data <- read.csv("path/to/directory/dataset.csv", sep = ";", header = TRUE)

# Display dataset summary
summary(data)
str(data)

# Boxplot for Bipolaris Disease Severity
ggplot(bipolaris, aes(x = Condition, y = Severity, fill = Condition)) +
  scale_fill_manual(values = c("myc" = "#788BFF", "ctrl" = "#44CC00")) +
  geom_boxplot(show.legend = TRUE) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(x = "Condition", y = "Disease Severity") +
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "white") +
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     hide.ns = FALSE, size = 5) + 
  theme_minimal()

# Boxplot for Pyricularia Disease Severity
ggplot(pyricularia, aes(x = Condition, y = Index, fill = Condition)) +
  scale_fill_manual(values = c("myc" = "#788BFF", "ctrl" = "#44CC00")) +
  geom_boxplot(show.legend = TRUE) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(x = "Condition", y = "Disease Severity") +   
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "white") +
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     hide.ns = FALSE, size = 5) + 
  theme_minimal()
