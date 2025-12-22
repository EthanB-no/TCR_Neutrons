library(immunarch)
library(dplyr)
library(immunarch)
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(lme4)

imm_data<- repLoad("./Combined_data/")
imm_meta <- read_csv("./combined_metadata_3.csv") 

keep_samples <- sapply(imm_data$data, function(x) sum(x$Clones) >= 5000)
filtered_data <- imm_data$data[keep_samples]

# Downsample to a fixed number of clones (e.g., 5000)
imm_data_sub <- repSample(filtered_data, .method = "downsample", .n = 5000)

# Check totals
sapply(imm_data_sub, function(x) sum(x$Clones))


# Diversity (Gini-Simpson)
ginisimp <- repDiversity(imm_data_sub, .method = "gini.simp")

clonality_annotated <- ginisimp %>%
  left_join(
    imm_meta %>% select(Sample, PatientID, Timepoint, Response, Dataset, UC_status, Arm, ),
    by = "Sample"
  ) %>% 
  mutate(Timepoint = factor(Timepoint,
                            levels = c("Pre-TX", "Post Rad", "4 Week", "16 Week"))) 

clonality_annotated <- subset(clonality_annotated, Timepoint != "Pre-TX"  | is.na(Timepoint))

clonality_annotated <-subset(clonality_annotated, UC_status != "BN") %>%
  mutate(Dataset = case_when(
    Dataset %in% c("FH") ~ "Hutch",
    Dataset %in% c("Z") ~ "Zhang",
    Dataset %in% c("Atezo") ~ "Atezo"))

clonality_annotated <-subset(clonality_annotated, Dataset != "Zhang") 

clonality_annotated <- subset(clonality_annotated, Arm == "concurrent" | Arm == "Hutch")
# Convergence
# Flatten + clean
combined <- bind_rows(purrr::imap(imm_data_sub, function(df, sname) {
  df %>%
    mutate(
      Sample = sname,
      PatientID = imm_meta$PatientID[imm_meta$Sample == sname],
      Timepoint = imm_meta$Timepoint[imm_meta$Sample == sname],
    )
})) %>% filter(!is.na(CDR3.aa), !is.na(CDR3.nt))

# Convergence summary
convergence_summary <- combined %>%
  group_by(PatientID, Timepoint, CDR3.aa) %>%
  summarise(num_unique_nt = n_distinct(CDR3.nt), .groups = "drop") %>%
  filter(num_unique_nt > 1) %>%
  count(PatientID, Timepoint, name = "num_convergent_clonotypes") %>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("Pre-TX", "Post Rad", "4 Week", "16 Week")))

convergence_summary <- subset(convergence_summary, Timepoint != "Pre-TX" | is.na(Timepoint))

# Merge with response info
convergence_with_response <- convergence_summary %>%
  left_join(imm_meta %>% select(PatientID, Dataset, UC_status, Arm), by = "PatientID") %>%
  mutate(Dataset = na_if(Dataset, "NA")) %>%  
  mutate(Dataset = case_when(
    Dataset %in% c("FH") ~ "Hutch",
    Dataset %in% c("Z") ~ "Zhang",
    Dataset %in% c("Atezo") ~ "Atezo"
    
  )) %>%
  distinct(PatientID, num_convergent_clonotypes, .keep_all = TRUE)

convergence_with_response <-subset( convergence_with_response, UC_status != "BN") 

convergence_with_response <-subset( convergence_with_response, Dataset != "Zhang") 
#convergence_with_response <-subset( convergence_with_response, Arm == "concurrent") 
# Helper function for plots
make_boxplot <- function(df, xvar, yvar, title, ylab) {
  ggplot(df, aes_string(x = xvar, y = yvar, fill = xvar)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, size = .8) +
    geom_jitter(aes(color = as.factor(PatientID)), show.legend = FALSE, width = 0.15, size = 4) +
    labs(title = title, x = "Cohort", y = ylab, color = "Cohort") +
    theme_classic2(base_size = 22) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      size = 10,
      bracket.size = 1.5,
      comparisons = list(c("Hutch", "Atezo")),
      hide.ns = FALSE,
      vjust = 0.1   # push the stars slightly downward
    )
}


# Figures
p1 <- make_boxplot(clonality_annotated,
                   "Dataset", "Value",
                   "Gini Index by Cohort",
                   "Gini-Simpson Value")

p2 <- make_boxplot(convergence_with_response,
                   "Dataset", "num_convergent_clonotypes",
                   "Convergent Clonotypes by Cohort",
                   "Number of Convergent Clonotypes")

# Combine with patchwork
p1 + p2

combined_plot <- p1 + p2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right", plot.title = element_text(face = "bold", size = 26)) &
  theme(axis.text = element_text(size = 26))
        

combined_plot

ggsave(
  "combined_COHORT_plot.png", 
  plot = combined_plot, 
  width = 15, height = 12, 
  dpi = 800  # 600 DPI is publication quality
)
###############################
# ARM PLOTS
###########################################
# Arm-level Gini & Convergence comparisons
###########################################

# Make sure Arm is a factor with correct order
clonality_annotated$Arm <- factor(clonality_annotated$Arm,
                                  levels = c("neoadjuvant", "concurrent", "Hutch"))
convergence_with_response$Arm <- factor(convergence_with_response$Arm,
                                        levels = c("neoadjuvant", "concurrent", "Hutch"))

# Function for Arm boxplots
make_arm_boxplot <- function(df, yvar, title, ylab) {
  ggplot(df, aes(x = Arm, y = .data[[yvar]], fill = Arm)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, size = 0.9) +
    geom_jitter(aes(color = as.factor(PatientID)), width = 0.15, size = 3.6) +
    labs(
      title = title,
      x = "Trial Arm",
      y = ylab,
      color = "Patient ID"
    ) +
    theme_classic2(base_size = 22) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      size = 8,
      hide.ns = FALSE,
      bracket.size = 1.3,
      comparisons = list(
        c("neoadjuvant", "concurrent"),
        c("neoadjuvant", "Hutch"),
        c("concurrent", "Hutch")
      )
    )
}

### 1. Gini-Simpson by Arm
p_gini_arm <- make_arm_boxplot(
  clonality_annotated,
  "Value",
  "Gini-Simpson Diversity by Trial Arm",
  "Gini-Simpson Index"
)

### 2. TCR Convergence by Arm
p_conv_arm <- make_arm_boxplot(
  convergence_with_response,
  "num_convergent_clonotypes",
  "TCR Convergence by Trial Arm",
  "Number of Convergent Clonotypes"
)

### 3. Combined panel
arm_combined_plot <- p_gini_arm + p_conv_arm +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Save
ggsave(
  "arm_comparison_plot.png",
  plot = arm_combined_plot,
  width = 19,
  height = 12,
  dpi = 800
)

arm_combined_plot
