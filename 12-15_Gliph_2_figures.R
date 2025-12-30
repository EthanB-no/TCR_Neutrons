library(immunarch)
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)

#imm_data <- repLoad("./Combined_data/")
#imm_meta <- read_xlsx("PatientID.xlsx") 
gliph_2_clusters <- read_csv("GLIPH_output_new.csv")

timepoint_levels <- c("Pre-TX", "Post Rad", "4 Week", "16 Week")

# Step 1: Compute patient-normalized proportions
cluster_patient <- gliph_2_clusters %>%
  mutate(
    Patient = sub(":.*", "", Sample),
    Timepoint = str_trim(sub(".*:", "", Sample)),
    Timepoint = factor(Timepoint, levels = timepoint_levels)
  ) %>%
  filter(pattern != "single") %>%
  filter(!is.na(Timepoint)) %>%   # drops malformed labels safely
  group_by(Patient, Timepoint, pattern) %>%
  summarise(total_freq = sum(Freq, na.rm = TRUE), .groups = "drop") %>%
  group_by(Patient, Timepoint) %>%
  mutate(prop = total_freq / sum(total_freq)) %>%
  ungroup()


cluster_timepoint <- cluster_patient %>%
  group_by(Timepoint, pattern) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    n_patients = n_distinct(Patient),
    .groups = "drop"
  )

top15_per_timepoint <- cluster_timepoint %>%
  group_by(Timepoint) %>%
  slice_max(order_by = mean_prop, n = 8, with_ties = FALSE) %>%
  ungroup()

p1 <- ggplot(top15_per_timepoint,
       aes(x = Timepoint, y = mean_prop, fill = pattern)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Timepoint",
    y = "Mean Cluster Abundance",
    fill = "GLIPH-2 motif",
    title = "Top 8 GLIPH-2 clusters per timepoint"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
p1 
ggsave(
  filename =("Gliph-2_clusters_over_time.png"),
  plot = p1,
  dpi = 400,
  width = 6,
  height = 6,
  units = "in"
)