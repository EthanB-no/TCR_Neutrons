library(immunarch)
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)


imm_data <- repLoad("./Combined_data/")
imm_meta <- read_xlsx("PatientID.xlsx") 


imm_data_filt <- imm_data
imm_data_filt$data <- imm_data$data[!grepl("_TCRB$", names(imm_data$data))]

# Convergence
# Flatten + clean
combined <- bind_rows(purrr::imap(imm_data_filt$data, function(df, sname) {
  df %>%
    mutate(
      Sample = sname,
      PatientID = imm_meta$PatientID[imm_meta$Sample == sname],
      Timepoint = imm_meta$Timepoint[imm_meta$Sample == sname]
    )
})) %>% filter(!is.na(CDR3.aa), !is.na(CDR3.nt))

convergent_clonotypes <- combined %>%
  group_by(PatientID, Timepoint, CDR3.aa) %>%
  filter(n_distinct(CDR3.nt) > 1) %>%   # convergence definition
  summarise(
    nt_sequences = list(unique(CDR3.nt)),
    num_unique_nt = n_distinct(CDR3.nt),
    total_clones = sum(Clones),
    .groups = "drop"
  )


cluster_burden <- convergent_clonotypes %>%
  group_by(PatientID, Timepoint, type) %>%
  summarise(
    n_convergent = n(),
    .groups = "drop"
  )

cohort_conv <- cluster_burden %>%
  group_by(Timepoint) %>%
  summarise(
    median_convergent = median(n_convergent),
    mean_convergent   = mean(n_convergent),
    n_samples         = n(),
    .groups = "drop"
  )

cluster_burden_no_pre <- cluster_burden %>%
  filter(Timepoint != "Pre-TX")

cluster_burden_no_pre <- cluster_burden_no_pre %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c("Post Rad", "4 Week", "16 Week")
    )
  )

ggplot(cluster_burden_no_pre,
       aes(x = Timepoint, y = n_convergent)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  theme_classic() +
  labs(
    y = "Convergent clonotypes per timepointx",
    x = "Timepoint",
    title = "Cohort-level convergent TCRs across treatment"
  )
