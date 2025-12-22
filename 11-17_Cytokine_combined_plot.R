library(immunarch)
library(ggplot2)
library(tidyverse)
library(stringr)
library(dplyr)
library(readxl)
library(readr)
library(tidyr)
library(ggpubr)

#read in
meta <- read_xlsx("PatientID.xlsx") 
imm_data <- repLoad("./Combined_data//")
CK_data <- read_csv("CK_data(Sheet1).csv") #cytokine data



CK1 <- CK_data[CK_data$Response != "NA", ] %>%
  mutate(
    IL13 = as.numeric(IL13)
  ) %>%
  filter(!is.na(Response_group))   # <--- removes NA rows


p1 <- ggplot(CK1, aes(x = Response_group, y = IL2, fill = Response_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(
    aes(color = as.factor(PatientID)),
    width = 0.15,
    size = 2
  ) +
  
  labs(
    title = "IL2 by Response Group",
    x = "Response Group",
    y = "IL2",
    color = "Patient ID"
    ) +
  
  theme_classic2(base_size = 22) +
  theme(legend.position = "right") + 
  stat_compare_means(
    
    method = "wilcox.test",
    label = "p.signif",
    size = 14,
    bracket.size = 1.5,
    comparisons = list(c("Responder", "Non-responder")),
    hide.ns = TRUE,
    vjust = .5   # push the stars slightly downward
  )

p2 <- ggplot(CK1, aes(x = Response_group, y = IL6, fill = Response_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(
    aes(color = as.factor(PatientID)),
    width = 0.15,
    size = 2
  ) +
  
  labs(
    title = "IL6 by Response Group",
    x = "Response Group",
    y = "IL6",
    color = "Patient ID"
  ) +
  
  theme_classic2(base_size = 22) +
  theme(legend.position = "right")

p3 <- ggplot(CK1, aes(x = Response_group, y = IL10, fill = Response_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(
    aes(color = as.factor(PatientID)),
    width = 0.15,
    size = 2
  ) +
  
  labs(
    title = "IL10 by Response Group",
    x = "Response Group",
    y = "IL10",
    color = "Patient ID"
  ) +
  
  theme_classic2(base_size = 22) +
  theme(legend.position = "right")


p4 <- ggplot(CK1, aes(x = Response_group, y = IL13, fill = Response_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(
    aes(color = as.factor(PatientID)),
    width = 0.15,
    size = 2
  ) +
  
  labs(
    title = "IL13 by Response Group",
    x = "Response Group",
    y = "IL13",
    color = "Patient ID", 
  ) +
  
  theme_classic2(base_size = 22) +
  theme(legend.position = "right")

library(patchwork)

p1 <- p1 + guides(fill = "none")
p2 <- p2 + guides(fill = "none")
p3 <- p3 + guides(fill = "none")
p4 <- p4 + guides(fill = "none")

combined_plot <- (p1 + p2) / (p3 + p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", guides(fill = "none"), plot.title = element_text(face = "bold"))

combined_plot

ggsave(
  "CK_combined_plot.png", 
  plot = combined_plot, 
  width = 12, height = 12, 
  dpi = 800,
  bg = "white" # 600 DPI is publication quality
)
