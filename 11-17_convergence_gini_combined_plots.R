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
# Flatten + clean
combined <- bind_rows(purrr::imap(imm_data_filt$data, function(df, sname) {
  df %>%
    mutate(
      Sample = sname,
      PatientID = imm_meta$PatientID[imm_meta$Sample == sname],
      Timepoint = imm_meta$Timepoint[imm_meta$Sample == sname]
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

# Merge with response info
convergence_with_response <- convergence_summary %>%
  left_join(imm_meta %>% select(PatientID, Response), by = "PatientID") %>%
  mutate(Response = na_if(Response, "NA")) %>%  
  mutate(Response_Group = case_when(
    Response %in% c("CR", "PR") ~ "Responder",
    Response %in% c("SD", "PD") ~ "Non-responder"
  )) %>%
  distinct(PatientID, num_convergent_clonotypes, .keep_all = TRUE) %>%
  filter(!is.na(Response_Group)) 


# Helper function for plots
make_boxplot <- function(df, xvar, yvar, title, ylab) {
  ggplot(df, aes_string(x = xvar, y = yvar, fill = xvar)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, size = .8) +
    geom_jitter(aes(color = as.factor(PatientID)), width = 0.15, size = 4) +
    labs(title = title, x = "Response Group", y = ylab, color = "Patient ID") +
    theme_classic2(base_size = 22) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      size = 14,
      bracket.size = 1.5,
      comparisons = list(c("Responder", "Non-responder")),
      hide.ns = TRUE,
      vjust = .5   # push the stars slightly downward
    )
}


# Figures
p1 <- make_boxplot(gini_response,
                   "Response_Group", "Value",
                   "Gini-Simpson Index by Response",
                   "Gini-Simpson Value")

p2 <- make_boxplot(convergence_with_response,
                   "Response_Group", "num_convergent_clonotypes",
                   "Convergent Clonotypes by Response",
                   "Number of Convergent Clonotypes")

# Combine with patchwork
p1 + p2

combined_plot <- p1 + p2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right") 
combined_plot

ggsave(
  "combined_plot.png", 
  plot = combined_plot, 
  width = 15, height = 10, 
  dpi = 800  # 600 DPI is publication quality
)
###############################
# Longitudinal line plots


# Gini-Simpson over time
p_line_gini <- ggplot(gini_response, 
                      aes(x = Timepoint, y = Value, group = PatientID, color = Response_Group)) +
  geom_line(aes(linetype = Response_Group), size = 1.2, alpha = 0.8) +
  geom_point(size = 4) +
  scale_x_discrete(limits = c("Pre-TX", "Post Rad", "4 Week", "16 Week")) +
  labs(title = "Longitudinal Gini-Simpson Diversity", 
       y = "Gini-Simpson Index", x = "Timepoint") +
  theme_classic(base_size = 20)

# Convergence over time
p_line_conv <- ggplot(convergence_with_response, 
                      aes(x = Timepoint, y = num_convergent_clonotypes, group = PatientID, color = Response_Group)) +
  geom_line(aes(linetype = Response_Group), size = 1.2, alpha = 0.8) +
  geom_point(size = 4) +
  scale_x_discrete(limits = c("Pre-TX", "Post Rad", "4 Week", "16 Week")) +
  labs(title = "Longitudinal Convergent Clonotypes", 
       y = "# Convergent Clonotypes", x = "Timepoint") +
  theme_classic(base_size = 20)

# Combine
p_line_gini + p_line_conv


