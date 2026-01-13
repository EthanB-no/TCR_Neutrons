library(immunarch)
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(ggalluvial)

#imm_data <- repLoad("./Combined_data/")
#imm_meta <- read_xlsx("PatientID.xlsx") 
gliph_2_clusters <- read_csv("GLIPH_output_new.csv")

timepoint_levels <- c("Pre-TX", "Post Rad", "4 Week", "16 Week")

# Step 1: Compute patient-normalized proportions
# Step 1: Compute patient-normalized proportions
cluster_freq <- gliph_2_clusters %>%
  mutate(
    Patient = sub(":.*", "", Sample),
    Timepoint = str_trim(sub(".*:", "", Sample)),
    Timepoint = factor(Timepoint, levels = timepoint_levels)
  ) %>%
  filter(pattern != "single", !is.na(Timepoint)) %>%
  group_by(Patient, Timepoint, pattern) %>%
  summarise(total_freq = sum(Freq, na.rm = TRUE), .groups = "drop")


top_10_timepoint <- cluster_freq %>%
  group_by(Patient, Timepoint) %>%
  mutate(rank = rank(-total_freq, ties.method = "min")) %>% # largest first
  filter(rank <= 10)

library(ggalluvial)
plot_patient_sankey <- function(patient_id) {
  ggplot(
    top_10_timepoint %>% filter(Patient == patient_id),
    aes(
      x = Timepoint,
      stratum = pattern,
      alluvium = pattern,
      y = total_freq,
      fill = pattern       # <-- this maps colors to pattern
    )
  ) +
    geom_flow(alpha = 0.65, color = "grey40", width = 0.15) +
    geom_stratum(width = 0.35, color = "black") +
    labs(
      title = paste("Top GLIPH-2 clusters over time: Patient ", patient_id),
      x = "Timepoint",
      y = "Total frequency",
      fill = "Pattern"     # <-- this sets the legend title
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.2),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}



patients <- unique(top_10_timepoint$Patient)

for (p in patients) {
  g <- plot_patient_sankey(p)
  
  # Save each plot as PNG
  ggsave(
    filename = paste0("Sankey_top10_", p, ".png"),
    plot = g,
    width = 6,
    height = 6,
    dpi = 400
  )
  
  # Optional: print in R for quick inspection
  print(g)
}
