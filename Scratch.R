library(dplyr)
library(readxl)

##NOTE this script is not used for any figures in the paper
imm_meta <- read_xlsx("PatientID.xlsx") 

sample_meta <- read.delim("./SampleOverview_11-19-2025_10-47-53_PM.tsv")

head(imm_meta)

head(sample_meta)

sample_name_tag <- sample_meta |> select(sample_name, sample_tags)

head(sample_name_tag)

combined <- dplyr::bind_rows(sample_name_tag, imm_meta)

write.csv(combined, "./combined_metadata_draft.csv")

combined_2 <- read_csv("./combined_metadata_draft.csv")

head(combined_2)

library(dplyr)

combined_3 <- combined_2 %>%
  group_by(Dataset) %>%
  mutate(PatientID = if_else(Dataset == "Atezo",
                             paste0("A", dense_rank(as.numeric(PatientID))),
                             PatientID))

write.csv(combined_3, "./combined_metadata_3.csv")