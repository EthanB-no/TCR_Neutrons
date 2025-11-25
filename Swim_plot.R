library(ggplot2)
library(dplyr)
library(readr)

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
clinical_data <- clinical_data %>% 
  filter(!is.na(PatientID))
# Plot
ggplot(clinical_data, aes(y = PatientID)) +
  
  # OS segment (gray background)
  geom_segment(aes(x = 0, xend = OS_months, y = PatientID, yend = PatientID),
               color = "gray85", size = 2.5) +
  
  # PFS segment (colored foreground)
  geom_segment(aes(x = 0, xend = PFS_months, y = PatientID, yend = PatientID, color = BestResponse),
               size = 3.5) +
  
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
                     name = "OS Status") +
  
  # Aesthetic polish
  labs(
    x = "Months",
    y = "Patient ID",
    color = "iRECIST Response",
    title = "Response and Survival Timeline"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.y = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  guides(shape = guide_legend(override.aes = list(color = "black")))

#########

timeline_plot <- ggplot(clinical_data, aes(y = PatientID)) +
  
  # OS segment (gray background line)
  geom_segment(aes(x = 0, xend = OS_months, y = PatientID, yend = PatientID),
               color = "gray85", size = 8) +
  
  # PFS segment (colored foreground)
  geom_segment(aes(x = 0, xend = PFS_months, y = PatientID, yend = PatientID, color = BestResponse),
               size = 10) +
  
  # PFS progression markers
  geom_point(
    data = clinical_data %>% filter(progression_status == "progressed"),
    aes(x = PFS_months, y = PatientID),
    shape = 21, fill = "black", color = "white", size = 5
  ) +
  
  # OS event/censoring marker
  geom_point(
    aes(x = OS_months, y = PatientID, shape = status),
    size = 5, color = "black"
  ) +
  
  # Manual color for iRECIST
  scale_color_manual(values = c(
    "CR" = "#6a3d9a",  
    "PR" = "#1f78b4",  
    "SD" = "#33a02c",  
    "PD" = "#e31a1c",  
    "Unknown" = "orange"
  )) +
  
  # Manual shape for OS censoring/event
  scale_shape_manual(values = c("0" = 1, "1" = 16),
                     labels = c("Censored", "Event"),
                     name = "OS Status") +
  
  labs(
    x = "Months",
    y = "Patient ID",
    color = "iRECIST Response",
    title = "Response and Survival Timeline"
  ) +
  theme_minimal(base_size = 28) +
  theme(
    axis.title = element_text(size = 30, face = "bold"),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 26, face = "bold"),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 34, face = "bold", hjust = 0.5),
    
    # 👇 force background to white
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  guides(
    shape = guide_legend(override.aes = list(color = "black", size = 5)),
    color = guide_legend(override.aes = list(size = 5))
  ) 

timeline_plot + 
  coord_cartesian(xlim = c(0, 20))



# Save high-res
ggsave(
  "response_survival_timeline.png", 
  plot = timeline_plot + 
    coord_cartesian(xlim = c(0, 20)), 
  width = 16, height = 10, 
  dpi = 900
)
