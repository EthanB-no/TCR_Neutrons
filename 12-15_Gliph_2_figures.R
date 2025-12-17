library(immunarch)
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)



imm_data <- repLoad("./Combined_data/")
imm_meta <- read_xlsx("PatientID.xlsx") 
gliph_2_clusters <- read_csv("GLIPH_output_new.csv")

gliph_2_clusters <- select(gliph_2_clusters, Fisher_score, Sample, TcRb, type)

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

convergent_clusters <- inner_join(convergent_clonotypes, gliph_2_clusters, by = c("CDR3.aa" = "TcRb"))


