library(immunarch)
library(ggplot2)
library(tidyverse)
library(stringr)
library(dplyr)
library(readxl)
library(readr)


imm_data <- repLoad("./Combined_data//")

meta <- read_xlsx("PatientID.xlsx") 
write.table(meta, file = "./Combined_data//metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)



imm_meta <- imm_data$meta
imm_meta$Sample <- names(imm_data$data)
meta_new <- meta

# Compute gini simpson
ginisimp <- repDiversity(imm_data$data, .method = "gini.simp")

# Merge with meta
clonality_annotated <- ginisimp %>%
  
  left_join(meta_new, by = "Sample") %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c("Pre-TX", "Post Rad", "4 Week", "16 Week")
    )
  )


gini_response <- clonality_annotated %>%
  mutate(Response_group = ifelse(Response == "PD", "Progression", "Non-Progression"))

# Make a new grouping column
gini_response <- gini_response %>%
  filter(!is.na(Response_group)) %>%  # Remove rows with NA in Response
  mutate(Response_group = case_when(
    Response %in% c("CR", "PR") ~ "Responder",
    Response %in% c("PD", "SD") ~ "Non-responder"
  ))


# Subset to Pre-TX and remove NAs
gini_pretx_clean <- gini_response %>%
  filter(Timepoint == "Pre-TX", !is.na(Response_group))

# Plot
ggplot(gini_pretx_clean, aes(x = Response_group, y = Value, fill = Response_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(
    aes(color = as.factor(PatientID)),
    width = 0.15,
    size = 2
  ) +
  labs(
    title = "Gini-Simpson Index by Response Group (Pre-TX)",
    x = "Response Group",
    y = "Gini-Simpson Value",
    color = "Patient ID"
  ) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "right")

#Wilcox test
wilcox_pretx <- wilcox.test(Value ~ Response_group, data = gini_pretx_clean)
print(wilcox_pretx)

library(ggpubr)

p1 <- ggplot(gini_pretx_clean, aes(x = Response_group, y = Value, fill = Response_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = as.factor(PatientID)), width = 0.15, size = 4) +
  labs(
    title = "Gini-Simpson Index by Response Group (Pre-TX)",
    x = "Response Group",
    y = "Gini-Simpson Value",
    color = "Patient ID"
  ) +
  theme_classic2(base_size = 20) +
  theme(legend.position = "right") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",    # <-- shows p-value instead of stars
    label.y = max(gini_pretx_clean$Value) * 1.001,
    size = 8,
    bracket.size = 1.5,
    vjust = .05,
    comparisons = list(c("Responder", "Non-responder"))
  )
p1
ggsave(
  "Gini_Pre_tx.png", 
  plot = p1, 
  width = 10, height = 8, 
  dpi = 800,
  bg = "transparent" # 600 DPI is publication quality
)