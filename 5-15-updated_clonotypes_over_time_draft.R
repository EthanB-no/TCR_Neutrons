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
  filter(Timepoint %in% c("Pre-TX", "4 Week", "Post Rad"))

patient_ids <- unique(clonotype_annotated$PatientID)

#data storage
all_clonotype_status <- list()

for (pid in patient_ids){
  patient_data <- clonotype_annotated %>% filter(PatientID == pid) %>%
    # Select key columns
    select(PatientID, Timepoint, CDR3.aa, Clones)
  
  clonotype_for_pivot <- patient_data %>%
    group_by(Timepoint, CDR3.aa) %>%
    summarise(Clones = sum(Clones), .groups = "drop")
  
  message("Processing patient: ", pid)
  
  
  
  clonotype_tracked_patient <- clonotype_for_pivot %>%
    
    pivot_wider(
      names_from = Timepoint,
      values_from = Clones,
      names_prefix = "Clones_"
    ) 
  
  if (!all(c("Clones_Pre-TX", "Clones_4 Week", "Clones_Post Rad") %in% colnames(clonotype_tracked_patient))) {
    message("Skipping patient ", pid, ": missing timepoints")
    next
  }
  
  # Handle missing clonotypes (replace NA with 0)
  clonotype_tracked_patient <- clonotype_tracked_patient %>% mutate(
    `Clones_Pre-TX` = replace_na(`Clones_Pre-TX`, 0),
    `Clones_4 Week` = replace_na(`Clones_4 Week`, 0), 
    `Clones_Post Rad` = replace_na(`Clones_Post Rad`, 0)
  )
  
  
  if (sum(clonotype_tracked_patient$`Clones_Pre-TX`) == 0 ||
      sum(clonotype_tracked_patient$`Clones_4 Week`) == 0 ||
      sum(clonotype_tracked_patient$`Clones_Post Rad`) == 0){
    message("Skipping patient ", pid, ": zero clones")
    next
  }
  # Calculate Total Clones for the patient at each time point
  N_pretx <- sum(clonotype_tracked_patient$`Clones_Pre-TX`)
  N_4week <- sum(clonotype_tracked_patient$`Clones_4 Week`)
  N_PostRad <- sum(clonotype_tracked_patient$`Clones_Post Rad`)
  
  clonotype_analysis_patient <- clonotype_tracked_patient %>%
    mutate(
      
      # Calculate Fold Change (FC) IN RAW COUNTS
      FC_4week = log2((`Clones_4 Week` + 1) / (`Clones_Pre-TX` + 1)),
      FC_postrad = log2((`Clones_Post Rad` + 1) / (`Clones_Pre-TX` + 1)),
      
      
      # Classify the clonotype's fate for post RAD
      Status_postrad = case_when(
        # Significantly Expanded/Contracted (must be present in both)
        FC_postrad > log2(1.5) & `Clones_Pre-TX` >= 1 & `Clones_Post Rad` >= 1 ~ "Expanded (FC > 1.5)",
        FC_postrad < log2(0.5) & `Clones_Pre-TX` >= 1 & `Clones_Post Rad` >= 1 ~ "Contracted (FC < 0.5)",
        # New/Lost clonotypes (only in one time point)
        `Clones_Pre-TX` == 0 & `Clones_Post Rad` > 0 ~ "New (Post-TX)",
        `Clones_Pre-TX` > 0 & `Clones_Post Rad` == 0 ~ "Lost (Post-TX)",
        .default = "Unchanged"
      ),
      # Classify the clonotype's fate for 4 week 
      Status_4week = case_when(
        # Significantly Expanded/Contracted (must be present in both)
        FC_4week > log2(1.5) & `Clones_Pre-TX` >= 1 & `Clones_4 Week` >= 1 ~ "Expanded (FC > 1.5)",
        FC_4week < log2(0.5) & `Clones_Pre-TX` >= 1 & `Clones_4 Week` >= 1 ~ "Contracted (FC < 0.5)",
        # New/Lost clonotypes (only in one time point)
        `Clones_Pre-TX` == 0 & `Clones_4 Week` > 0 ~ "New (Post-TX)",
        `Clones_Pre-TX` > 0 & `Clones_4 Week` == 0 ~ "Lost (Post-TX)",
        .default = "Unchanged"
      )
      
    ) %>%
    
    # Apply Log Transformation for Plotting
    mutate(
      Log_Clones_PreTX = log10(`Clones_Pre-TX` + 1),
      Log_Clones_4Week = log10(`Clones_4 Week` + 1),
      Log_Clones_Postrad = log10(`Clones_Post Rad` + 1)
    )
  
  
  
  
  # 5. Create the Colored Scatter Plot
  max_val <- max(clonotype_analysis_patient$Log_Clones_PreTX, clonotype_analysis_patient$Log_Clones_4Week)
  
  scatter_plot_colored_4week <- clonotype_analysis_patient %>%
    ggplot(aes(x = Log_Clones_PreTX, y = Log_Clones_4Week, color = Status_4week)) +
    
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
    theme(plot.title = element_text(hjust = .15, face = "bold", size = 18)) +
    theme(legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12))
  
  print(scatter_plot_colored_4week) 
  
  
  
  ggsave(
    filename = paste0("patient_", pid, "_clonotype_tracking_4week.png"),
    plot = scatter_plot_colored_4week,
    dpi = 400,
    width = 8,
    height = 8,
    units = "in")
    
    ###############PLOT for POSTRAD timepoint########################################
    scatter_plot_colored_POSTRAD <- clonotype_analysis_patient %>%
      ggplot(aes(x = Log_Clones_PreTX, y = Log_Clones_Postrad, color = Status_postrad)) +
      
      geom_point(alpha = 0.6, size = 2.5) +
      
      # Identity line
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
      
      # Labels
      labs(
        title = paste("Clonotype Tracking (Pre-TX vs Postrad) Patient", pid),
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
      theme(plot.title = element_text(hjust = .15, face = "bold", size = 18)) +
      theme(legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 12))
    
    print(scatter_plot_colored_POSTRAD) 
    
    
    
    ggsave(
      filename = paste0("patient_", pid, "_clonotype_trackingPOSTRAD.png"),
      plot = scatter_plot_colored_POSTRAD,
      dpi = 400,
      width = 8,
      height = 8,
      units = "in"
  )
  
  
  # Add patient column
  clonotype_analysis_patient$PatientID <- pid
  
  # Store in the list
  all_clonotype_status[[pid]] <- clonotype_analysis_patient %>%
    select(PatientID, CDR3.aa, Status_4week, Status_postrad)
  
}

cohort_status <- bind_rows(all_clonotype_status) %>%
  pivot_longer(
    cols = starts_with("Status_"),
    names_to = "Timepoint",
    values_to = "Status"
  ) %>%
  mutate(
    Timepoint = recode(
      Timepoint,
      "Status_4week"   = "4 Week",
      "Status_postrad" = "Post Rad"
    ),
    FateGroup = case_when(
      Status %in% c("Expanded (FC > 1.5)", "New (Post-TX)") ~ "Expanded/New",
      Status %in% c("Contracted (FC < 0.5)", "Lost (Post-TX)") ~ "Contracted/Lost",
      TRUE ~ "Unchanged"
    ),
    Subgroup = case_when(
      Status == "Expanded (FC > 1.5)" ~ "Expanded",
      Status == "New (Post-TX)" ~ "New",
      Status == "Contracted (FC < 0.5)" ~ "Contracted",
      Status == "Lost (Post-TX)" ~ "Lost",
      TRUE ~ NA_character_
    )
  )

status_summary <- cohort_status %>%
  filter(FateGroup != "Unchanged") %>%   # optional
  group_by(Timepoint, FateGroup, Subgroup) %>%
  summarise(Count = n(), .groups = "drop")


p1 <- ggplot(
  status_summary,
  aes(x = FateGroup, y = Count, fill = Subgroup)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~ Timepoint) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  scale_fill_manual(values = c(
    "New" = "red",
    "Expanded" = "orange",
    "Lost" = "black",
    "Contracted" = "royalblue"
  )) +
  labs(
    title = "Cohort-level Clonotype Fate Across Timepoints",
    x = "Clonotype Fate Group",
    y = "# of Clonotypes (log10)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )

p1


ggsave(
  filename = "cohort_level_clonotype_fate_all_timepoints.png",
  plot = p1,
  dpi = 400,
  width = 9,
  height = 5,
  units = "in"
)

