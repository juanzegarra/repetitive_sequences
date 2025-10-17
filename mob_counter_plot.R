#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(viridis)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

df <- read.table(input_file, header = TRUE, sep = "\t")

df_sum <- df %>%
  group_by(Family) %>%
  summarise(Value = sum(Copies))

p <- ggplot(df_sum, aes(
    x = reorder(Family, Value),
    y = Value,
    fill = Family
)) +
  geom_col(position = "identity", width = 0.64) +
  labs(
    title = element_blank(),
    x = "Famílias de Transposons",
    y = "N° de linhagens"
  ) +
  theme(
    panel.background = element_rect(fill = "#F5F5F5"),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
    axis.ticks = element_blank()
  ) +
  scale_fill_viridis_d(option = "plasma")

ggsave(output_file, p, width = 8, height = 6, dpi = 300)

