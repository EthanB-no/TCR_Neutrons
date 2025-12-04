library(ggplot2)
library(dplyr)
library(readr)
library(survminer)
library(survival)


##NOTE this script is note used for any figures in the paper
# Load data
clinical_data <- read_csv("./swimmer_plot.csv")

# Prepare the data
clinical_data <- clinical_data %>%
  mutate(
    PatientID = paste0("P", sprintf("%03d", `Patient No.`)),
    BestResponse = `Best Response per iRECIST`,
    BestResponse = ifelse(is.na(BestResponse) | BestResponse == "", "Unknown", BestResponse),
    BestResponse = factor(BestResponse, levels = c("CR", "PR", "SD", "PD", "Unknown")),
    status = factor(as.character(status), levels = c("0", "1")),
    progression_status = factor(progression_status, levels = c("censored", "progressed"))
  ) %>%
  arrange(OS_months) %>%
  mutate(PatientID = factor(PatientID, levels = PatientID))

# Plot
p1 <- ggplot(clinical_data, aes(y = PatientID)) +
  
  # OS segment (gray background)
  geom_segment(aes(x = 0, xend = OS_months, y = PatientID, yend = PatientID),
               color = "gray85", size = 8) +
  
  # PFS segment (colored foreground)
  geom_segment(aes(x = 0, xend = PFS_months, y = PatientID, yend = PatientID, color = BestResponse),
               size = 8) +
  
  # Only show PFS progression markers where patients progressed
  geom_point(
    data = clinical_data %>% filter(progression_status == "progressed"),
    aes(x = PFS_months, y = PatientID),
    shape = 21, fill = "black", color = "white", size = 3
  ) +
  
  # OS event/censoring marker at OS_months
  geom_point(
    aes(x = OS_months, y = PatientID, shape = status),
    size = 3, color = "black"
  ) +
  
  # Manual color for iRECIST
  scale_color_manual(values = c(
    "CR" = "#6a3d9a",  # purple
    "PR" = "#1f78b4",  # blue
    "SD" = "#33a02c",  # green
    "PD" = "#e31a1c",  # red
    "Unknown" = "orange"
  )) +
  
  # Manual shape for OS censoring/event
  scale_shape_manual(values = c("0" = 1, "1" = 16),
                     labels = c("Censored", "Event"),
                     name = "Progression from Baseline Status") +
  
  # Aesthetic polish
  labs(
    x = "Months",
    y = "Patient ID",
    color = "iRECIST Response",
    title = "Response and Survival Timeline"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  guides(shape = guide_legend(override.aes = list(color = "black")))

ggsave(
  "Swimmmer_plot_new.png", 
  plot = print(p1), 
  width = 12, height = 8, 
  dpi = 800,
  bg = "white" # 600 DPI is publication quality
)

####################################
#making kaplan meier curve
clinical_data <- clinical_data %>%
  mutate(
    # Convert factors to numeric binary event indicators
    os_event = as.numeric(as.character(status)),  # 1 = dead, 0 = censored
    prog_event = ifelse(progression_status == "progressed", 1, 0),
    
    # Composite PFS event: progression or death
    pfs_event = ifelse(os_event == 1 | prog_event == 1, 1, 0)
  )


# Create Surv object and fit
surv_obj <- Surv(time = clinical_data$PFS_months, event = clinical_data$pfs_event)

# No grouping variable example
fit <- survfit(surv_obj ~ 1)

# Plot
#fit_strat <- survfit(surv_obj ~ BestResponse, data = clinical_data)

p1 <- ggsurvplot(fit, data = clinical_data,
           risk.table = TRUE,
           pval = FALSE,
           xlab = "Time (Months)",
           ylab = "Progression-Free Survival",
           title =  "Kaplan-Meier Curve: Progression Free Survival",
           conf.int = FALSE,
           legend = "none",
           censor = TRUE)

p1

ggsave(
  "Swimmmer_plot_new.png", 
  plot = print(p1), 
  width = 12, height = 8, 
  dpi = 800  # 600 DPI is publication quality
)
ggexport(p1, filename = "Swimmer_plot_new.png", width = 12000, height = 8000, res = 800)
############################



