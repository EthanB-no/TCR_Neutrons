library(tidyverse)
library(immunarch)
library(readxl)
library(ggpubr)

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

library(tidyverse)
library(ggpubr)


# Function to process a comparison between two timepoints
process_clonotype_comparison <- function(time1, time2) {
  
  # Use the already joined clonality_annotated (or join again if using raw clonotype data)
  # Bind the raw clonotype tables into one tibble
  clonotype_data <- bind_rows(imm_data$data, .id = "Sample") %>%
    mutate(Clones = as.numeric(Clones)) %>%
    left_join(imm_meta %>% select(Sample, PatientID, Timepoint), by = "Sample")
  
  patient_ids <- unique(clonotype_data$PatientID)
  all_status <- list()
  
  for (pid in patient_ids) {
    patient_data <- clonotype_data %>% filter(PatientID == pid)
    
    if (nrow(patient_data) == 0) next
    
    # Pivot wider
    clonotype_wide <- patient_data %>%
      group_by(Timepoint, CDR3.aa) %>%
      summarise(Clones = sum(Clones), .groups = "drop") %>%
      pivot_wider(names_from = Timepoint, values_from = Clones, names_prefix = "Clones_") %>%
      mutate(across(starts_with("Clones_"), ~replace_na(.x, 0)))
    
    if (!all(c(paste0("Clones_", time1), paste0("Clones_", time2)) %in% colnames(clonotype_wide))) next
    
    # Fold change and status
    clonotype_wide <- clonotype_wide %>%
      mutate(
        FC = log2((!!sym(paste0("Clones_", time2)) + 1) / (!!sym(paste0("Clones_", time1)) + 1)),
        Status = case_when(
          FC > log2(1.5) & !!sym(paste0("Clones_", time1)) >= 1 & !!sym(paste0("Clones_", time2)) >= 1 ~ "Expanded",
          FC < log2(0.5) & !!sym(paste0("Clones_", time1)) >= 1 & !!sym(paste0("Clones_", time2)) >= 1 ~ "Contracted",
          !!sym(paste0("Clones_", time1)) == 0 & !!sym(paste0("Clones_", time2)) > 0 ~ "New",
          !!sym(paste0("Clones_", time1)) > 0 & !!sym(paste0("Clones_", time2)) == 0 ~ "Lost",
          TRUE ~ "Unchanged"
        ),
        Comparison = paste(time1, "vs", time2)
      ) %>%
      select(PatientID, CDR3.aa, Status, Comparison)
    
    all_status[[pid]] <- clonotype_wide
  }
  
  bind_rows(all_status)
}

# Combine both comparisons
cohort_status <- bind_rows(
  process_clonotype_comparison("Pre-TX", "Post Rad"),
  process_clonotype_comparison("Pre-TX", "4 Week")
)

# Summarize for plotting
status_summary <- cohort_status %>%
  filter(Status != "Unchanged") %>%
  group_by(Comparison, Status) %>%
  summarise(Count = n(), .groups = "drop")

status_summary$Status <- factor(status_summary$Status, levels = c("New", "Expanded", "Lost", "Contracted"))

# Grouped bar chart: 4 bars (2 comparisons × 2 statuses)
p_combined <- ggplot(status_summary, aes(x = Comparison, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c(
    "New" = "red",
    "Expanded" = "orange",
    "Lost" = "black",
    "Contracted" = "royalblue"
  )) +
  labs(
    title = "Cohort-level Clonotype Fate",
    x = "Comparison",
    y = "# of Clonotypes",
    fill = "Fate"
  ) +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))

p_combined

ggsave("cohort_level_clonotype_fate_4bars.png", plot = p_combined, width = 8, height = 6, dpi = 400)
