library(immunarch)
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(lme4)


imm_data <- repLoad("./Combined_data/")
imm_meta <- read_xlsx("PatientID.xlsx") 

imm_data_filt <- imm_data
imm_data_filt$data <- imm_data$data[!grepl("_TCRB$", names(imm_data$data))]

# Diversity (Gini-Simpson)
ginisimp <- repDiversity(imm_data_filt$data, .method = "gini.simp")

#add info from metadata 
clonality_annotated <- ginisimp %>%
  left_join(
    imm_meta %>% select(Sample, PatientID, Timepoint, Response), 
    by = "Sample"
  ) %>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("Pre-TX", "Post Rad", "4 Week", "16 Week")))

#group by response
gini_response <- clonality_annotated %>%
  mutate(Response = na_if(Response, "NA")) %>%
  mutate(Response_Group = case_when(
    Response %in% c("CR", "PR") ~ "Responder",
    Response %in% c("PD", "SD") ~ "Non-responder"
  )) %>%
  filter(!is.na(Response_Group)) 

# Convergence
 