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

#keep_samples <- sapply(imm_data$data, function(x) sum(x$Clones) >= 5000)
#filtered_data <- imm_data$data[keep_samples]

# Downsample to a fixed number of clones (e.g., 5000)
#imm_data_sub <- repSample(filtered_data, .method = "downsample", .n = 5000)

# Check totals
#sapply(imm_data_sub, function(x) sum(x$Clones))


# Diversity (Gini-Simpson)
ginisimp <- repDiversity(imm_data$data, .method = "gini.simp")

clonality_annotated <- ginisimp %>%
  left_join(
    imm_meta %>% select(Sample, PatientID, Timepoint, Response, Dataset, UC_status),
    by = "Sample"
  ) %>% 
  mutate(Timepoint = factor(Timepoint,
                            levels = c("Pre-TX", "Post Rad", "4 Week", "16 Week"))) 


clonality_annotated <- subset(clonality_annotated, Dataset != "Pre-TX"  | is.na(Timepoint))
clonality_annotated <-subset(clonality_annotated, UC_status != "BN") %>%
  mutate(Dataset = case_when(
    Dataset %in% c("FH") ~ "Hutch",
    Dataset %in% c("Z") ~ "Zhang",
    Dataset %in% c("Atezo") ~ "Atezo"))

# Convergence
# Flatten + clean
combined <- bind_rows(purrr::imap(imm_data$data, function(df, sname) {
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
  left_join(imm_meta %>% select(PatientID, Dataset, UC_status), by = "PatientID") %>%
  mutate(Dataset = na_if(Dataset, "NA")) %>%  
  mutate(Dataset = case_when(
    Dataset %in% c("FH") ~ "Hutch",
    Dataset %in% c("Z") ~ "Zhang",
    Dataset %in% c("Atezo") ~ "Atezo"
    
  )) %>%
  distinct(PatientID, num_convergent_clonotypes, .keep_all = TRUE)

convergence_with_response <-subset( convergence_with_response, UC_status != "BN") 

# Helper function for plots
make_boxplot <- function(df, xvar, yvar, title, ylab) {
  ggplot(df, aes_string(x = xvar, y = yvar, fill = xvar)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, size = .8) +
    geom_jitter(aes(color = as.factor(PatientID)), width = 0.15, size = 4) +
    labs(title = title, x = "Cohort", y = ylab, color = "Patient ID") +
    theme_classic2(base_size = 22) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      size = 10,
      bracket.size = 1.5,
      comparisons = list(c("Zhang", "Hutch", "Atezo")),
      hide.ns = FALSE,
      vjust = 0.1   # push the stars slightly downward
    )
}


# Figures
p1 <- make_boxplot(clonality_annotated,
                   "Dataset", "Value",
                   "Gini-Simpson Index by Cohort raw",
                   "Gini-Simpson Value")

p2 <- make_boxplot(convergence_with_response,
                   "Dataset", "num_convergent_clonotypes",
                   "Convergent Clonotypes by Cohort raw",
                   "Number of Convergent Clonotypes")

# Combine with patchwork
p1 + p2

combined_plot <- p1 + p2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right") 
combined_plot

ggsave(
  "combined_plot_no_downample.png", 
  plot = combined_plot, 
  width = 17, height = 12, 
  dpi = 800  # 600 DPI is publication quality
)
###############################
# Longitudinal line plots