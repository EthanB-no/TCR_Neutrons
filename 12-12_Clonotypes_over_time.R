library(immunarch)
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(lme4)


imm_data <- repLoad("./Combined_data/")
imm_meta <- read_xlsx("PatientID.xlsx") 


# Diversity (Gini-Simpson)
ginisimp <- repDiversity(imm_data$data, .method = "gini.simp")

#add info from metadata 
clonality_annotated <- ginisimp %>%
  left_join(
    imm_meta %>% select(Sample, PatientID, Timepoint, Response), 
    by = "Sample"
  ) %>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("Pre-TX", "Post Rad", "4 Week", "16 Week")))



# 1. Extract the raw clonotype data
clonotype_data <- imm_data$data %>%
  bind_rows(.id = "Sample")

clonotype_data <- clonotype_data %>%
  mutate(Clones = as.numeric(Clones))

# 2. Join with metadata
clonotype_annotated <- clonotype_data %>%
  left_join(
    imm_meta %>% select(Sample, PatientID, Timepoint),
    by = "Sample"
  ) %>%
  # Filter for the relevant time points
  filter(Timepoint %in% c("Pre-TX", "4 Week"))

patient_ids <- unique(clonotype_annotated$PatientID)


patient_data <- clonotype_annotated %>% filter(PatientID == pid) %>%
  # Select key columns
  select(PatientID, Timepoint, CDR3.aa, Clones)

clonotype_for_pivot <- patient_data %>%
  group_by(Timepoint, CDR3.aa) %>%
  summarise(Clones = sum(Clones), .groups = "drop")

  

clonotype_tracked_patient <- clonotype_for_pivot %>%
  
  pivot_wider(
    names_from = Timepoint,
    values_from = Clones,
    names_prefix = "Clones_"
  ) %>%
  
  # Handle missing clonotypes (replace NA with 0)
  mutate(
    `Clones_Pre-TX` = replace_na(`Clones_Pre-TX`, 0),
    `Clones_4 Week` = replace_na(`Clones_4 Week`, 0)
  )

if (!all(c("Clones_Pre-TX", "Clones_4 Week") %in% colnames(clonotype_tracked_patient))) {
  message("Skipping patient ", pid, ": missing timepoint")
  next
}

if (sum(clonotype_tracked_patient$`Clones_Pre-TX`) == 0 ||
    sum(clonotype_tracked_patient$`Clones_4 Week`) == 0) {
  message("Skipping patient ", pid, ": zero clones")
  next
}
# Calculate Total Clones for the patient at each time point
N_pretx <- sum(clonotype_tracked_patient$`Clones_Pre-TX`)
N_4week <- sum(clonotype_tracked_patient$`Clones_4 Week`)

clonotype_analysis_patient <- clonotype_tracked_patient %>%
  mutate(
    # Calculate proportions
    Prop_PreTX = `Clones_Pre-TX` / N_pretx,
    Prop_4Week = `Clones_4 Week` / N_4week,
    
    # Define a minimum proportion (pseudo-count) for FC calculation stability
    min_prop = min(Prop_PreTX[Prop_PreTX > 0]) / 10,
    
    # Calculate Fold Change (FC)
    FC = (Prop_4Week + min_prop) / (Prop_PreTX + min_prop),
    
    # Classify the clonotype's fate
    Status = case_when(
      # Significantly Expanded/Contracted (must be present in both)
      FC > 1.5 & `Clones_Pre-TX` >= 1 & `Clones_4 Week` >= 1 ~ "Expanded (FC > 1.5)",
      FC < 0.5 & `Clones_Pre-TX` >= 1 & `Clones_4 Week` >= 1 ~ "Contracted (FC < 0.5)",
      # New/Lost clonotypes (only in one time point)
      `Clones_Pre-TX` == 0 & `Clones_4 Week` > 0 ~ "New (Post-TX)",
      `Clones_Pre-TX` > 0 & `Clones_4 Week` == 0 ~ "Lost (Post-TX)",
      .default = "Unchanged"
    )
  ) %>%
  
  # Apply Log Transformation for Plotting
  mutate(
    Log_Clones_PreTX = log10(`Clones_Pre-TX` + 1),
    Log_Clones_4Week = log10(`Clones_4 Week` + 1)
  )

# 5. Create the Colored Scatter Plot
max_val <- max(clonotype_analysis_patient$Log_Clones_PreTX, clonotype_analysis_patient$Log_Clones_4Week)

scatter_plot_colored <- clonotype_analysis_patient %>%
  ggplot(aes(x = Log_Clones_PreTX, y = Log_Clones_4Week, color = Status)) +
  
  geom_point(alpha = 0.6, size = 2.5) +
  
  # Identity line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  
  # Labels
  labs(
    title = paste("Clonotype Tracking (Pre-TX vs 4 Week) Patient", pid),
    x = expression(log[10]~"(Pre-TX Clones + 1)"),
    y = expression(log[10]~"(4 Week Clones + 1)"),
    color = "Clonotype Fate (FC Threshold)"
  ) +
  
  # Custom Colors
  scale_color_manual(
    values = c(
      "Expanded (FC > 1.5)" = "orange",
      "Contracted (FC < 0.5)" = "royalblue",
      "New (Post-TX)" = "red",
      "Lost (Post-TX)" = "black",
      "Unchanged" = "gray80"
    )
  ) +
  
  # Enforce equal axis limits
  coord_equal(xlim = c(0, max_val), ylim = c(0, max_val)) +
  
  theme_classic2(base_size = 18 ) +
  theme(plot.title = element_text(hjust = .15, face = "bold"))

print(scatter_plot_colored) 



ggsave(filename = "patient_1_clonotype_tracking.png", device = png, plot = scatter_plot_colored, dpi = 400, scale = 1, unit = "in", width = 8,
     height = 8)