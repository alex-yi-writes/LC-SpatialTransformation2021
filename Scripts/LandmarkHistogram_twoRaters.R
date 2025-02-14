#-------------------------------------------------------------------------------------------------#
#----------------------- multiple plots with both raters' ratings --------------------------------#
#-------------------------------------------------------------------------------------------------#
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(foreign)
library(R.matlab)
library(gridExtra)
library(png)
library(cowplot)
library(magick)

# Set working directory and load data
setwd("/your/path/to/working/directory/")
landmarks_Rater1 <- readMat('/your/path/to/Distance_EPIlandmarks_Rater1.mat', 
                            to.data.frame = TRUE)
landmarks_Rater2 <- readMat('/your/path/to/Distance_EPIlandmarks_Rater2.mat', 
                            to.data.frame = TRUE)
ED_Rater1 <- landmarks_Rater1$Distance.export
ED_Rater2 <- landmarks_Rater2$Distance.export

# Hard-code sample size
sampleSize <- nrow(ED_Rater2)

# Calculate means and SDs for each ROI for each rater
# Periaqueductal Grey (column 4)
m_TBS_R1 <- mean(ED_Rater1[,4])
sd_TBS_R1 <- sd(ED_Rater1[,4])
m_TBS_R2 <- mean(ED_Rater2[,4])
sd_TBS_R2 <- sd(ED_Rater2[,4])

# Outline Brainstem (L) (column 5)
m_OBSL_R1 <- mean(ED_Rater1[,5])
sd_OBSL_R1 <- sd(ED_Rater1[,5])
m_OBSL_R2 <- mean(ED_Rater2[,5])
sd_OBSL_R2 <- sd(ED_Rater2[,5])

# Outline Brainstem (R) (column 6)
m_OBSR_R1 <- mean(ED_Rater1[,6])
sd_OBSR_R1 <- sd(ED_Rater1[,6])
m_OBSR_R2 <- mean(ED_Rater2[,6])
sd_OBSR_R2 <- sd(ED_Rater2[,6])

# 4th Ventricle Border (L) (column 7)
m_LCL_R1 <- mean(ED_Rater1[,7])
sd_LCL_R1 <- sd(ED_Rater1[,7])
m_LCL_R2 <- mean(ED_Rater2[,7])
sd_LCL_R2 <- sd(ED_Rater2[,7])

# 4th Ventricle Border (R) (column 8)
m_LCR_R1 <- mean(ED_Rater1[,8])
sd_LCR_R1 <- sd(ED_Rater1[,8])
m_LCR_R2 <- mean(ED_Rater2[,8])
sd_LCR_R2 <- sd(ED_Rater2[,8])

# Perifastigial Sulcus (column 9)
m_BBS_R1 <- mean(ED_Rater1[,9])
sd_BBS_R1 <- sd(ED_Rater1[,9])
m_BBS_R2 <- mean(ED_Rater2[,9])
sd_BBS_R2 <- sd(ED_Rater2[,9])

# Create data frames for each ROI (both raters)
data_TBS <- data.frame(
  Rater = c(rep("Rater1", sampleSize), rep("Rater2", sampleSize)),
  Distance = c(ED_Rater1[,4], ED_Rater2[,4])
)
p3 <- ggplot(data_TBS, aes(x = Distance, fill = Rater)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 0.2) +
  scale_fill_manual(values = c("#ff0000", "#0000ff")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold")) +
  xlim(-0.1, 4.5) + ylim(0, 50) +
  labs(title = "Periaqueductal Grey", x = "mm", y = "Number of scans") +
  theme(legend.position = "none")

data_OBSL <- data.frame(
  Rater = c(rep("Rater1", sampleSize), rep("Rater2", sampleSize)),
  Distance = c(ED_Rater1[,5], ED_Rater2[,5])
)
p4 <- ggplot(data_OBSL, aes(x = Distance, fill = Rater)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 0.2) +
  scale_fill_manual(values = c("#ff0000", "#0000ff")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold")) +
  xlim(-0.1, 4.5) + ylim(0, 50) +
  labs(title = "Outline Brainstem (L)", x = "mm", y = "Number of scans") +
  theme(legend.position = "none")

data_OBSR <- data.frame(
  Rater = c(rep("Rater1", sampleSize), rep("Rater2", sampleSize)),
  Distance = c(ED_Rater1[,6], ED_Rater2[,6])
)
p5 <- ggplot(data_OBSR, aes(x = Distance, fill = Rater)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 0.2) +
  scale_fill_manual(values = c("#ff0000", "#0000ff")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold")) +
  xlim(-0.1, 4.5) + ylim(0, 50) +
  labs(title = "Outline Brainstem (R)", x = "mm", y = "Number of scans") +
  theme(legend.position = c(0.9, 0.25), legend.direction = "vertical",
        legend.background = element_blank(), legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

data_LCL <- data.frame(
  Rater = c(rep("Rater1", sampleSize), rep("Rater2", sampleSize)),
  Distance = c(ED_Rater1[,7], ED_Rater2[,7])
)
p6 <- ggplot(data_LCL, aes(x = Distance, fill = Rater)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 0.2) +
  scale_fill_manual(values = c("#ff0000", "#0000ff")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold")) +
  xlim(-0.1, 4.5) + ylim(0, 50) +
  labs(title = "4th Ventricle Border (L)", x = "mm", y = "Number of scans") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 2.5, color = "red", size = 1.2)

data_LCR <- data.frame(
  Rater = c(rep("Rater1", sampleSize), rep("Rater2", sampleSize)),
  Distance = c(ED_Rater1[,8], ED_Rater2[,8])
)
p7 <- ggplot(data_LCR, aes(x = Distance, fill = Rater)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 0.2) +
  scale_fill_manual(values = c("#ff0000", "#0000ff")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold")) +
  xlim(-0.1, 4.5) + ylim(0, 50) +
  labs(title = "4th Ventricle Border (R)", x = "mm", y = "Number of scans") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 2.5, color = "red", size = 1.2)

data_BBS <- data.frame(
  Rater = c(rep("Rater1", sampleSize), rep("Rater2", sampleSize)),
  Distance = c(ED_Rater1[,9], ED_Rater2[,9])
)
p8 <- ggplot(data_BBS, aes(x = Distance, fill = Rater)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 0.2) +
  scale_fill_manual(values = c("#ff0000", "#0000ff")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold")) +
  xlim(-0.1, 4.5) + ylim(0, 50) +
  labs(title = "Perifastigial Sulcus", x = "mm", y = "Number of scans") +
  theme(legend.position = c(0.9, 0.25), legend.direction = "vertical",
        legend.background = element_blank(), legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

# Add inset images with separate stats labels for each rater
p3_inset <- ggdraw(p3) +
  draw_image("/your/path/to/Landmark_Insets/PeriaqueductalGrey.png", 
              x = 0.73, y = 0.58, width = 0.3, height = 0.3) +
  draw_label(sprintf("Mean = %.2f, SD = %.2f", m_TBS_R1, sd_TBS_R1),
             x = 0.79, y = 0.71 + 0.12, hjust = 1, size = 10, color = "#ff0000") +
  draw_label(sprintf("Mean = %.2f, SD = %.2f", m_TBS_R2, sd_TBS_R2),
             x = 0.79, y = 0.71 + 0.07, hjust = 1, size = 10, color = "#0000ff")

p4_inset <- ggdraw(p4) +
  draw_image("/your/path/to/Landmark_Insets/OutlineBrainstem_left.png", 
             x = 0.73, y = 0.58, width = 0.3, height = 0.3) +
  draw_label(sprintf("Mean = %.2f, SD = %.2f", m_OBSL_R1, sd_OBSL_R1),
             x = 0.79, y = 0.71 + 0.12, hjust = 1, size = 10, color = "#ff0000") +
  draw_label(sprintf("Mean = %.2f, SD = %.2f", m_OBSL_R2, sd_OBSL_R2),
             x = 0.79, y = 0.71 + 0.07, hjust = 1, size = 10, color = "#0000ff")

p5_inset <- ggdraw(p5) +
  draw_image("/your/path/to/Landmark_Insets/OutlineBrainstem_right.png", 
             x = 0.73, y = 0.58, width = 0.3, height = 0.3) +
  draw_label(sprintf("Mean = %.2f, SD = %.2f", m_OBSR_R1, sd_OBSR_R1),
             x = 0.79, y = 0.71 + 0.12, hjust = 1, size = 10, color = "#ff0000") +
  draw_label(sprintf("Mean = %.2f, SD = %.2f", m_OBSR_R2, sd_OBSR_R2),
             x = 0.79, y = 0.71 + 0.07, hjust = 1, size = 10, color = "#0000ff")

p6_inset <- ggdraw(p6) +
  draw_image("/your/path/to/Landmark_Insets/LC_left.png", 
             x = 0.73, y = 0.58, width = 0.3, height = 0.3) +
  draw_label(sprintf("Mean = %.2f, SD = %.2f", m_LCL_R1, sd_LCL_R1),
            x = 0.79, y = 0.71 + 0.12, hjust = 1, size = 10, color = "#ff0000") +
  draw_label(sprintf("Mean = %.2f, SD = %.2f", m_LCL_R2, sd_LCL_R2),
             x = 0.79, y = 0.71 + 0.07, hjust = 1, size = 10, color = "#0000ff")

p7_inset <- ggdraw(p7) +
  draw_image("/your/path/to/Landmark_Insets/LC_right.png", 
             x = 0.73, y = 0.58, width = 0.3, height = 0.3) +
  draw_label(sprintf("Mean = %.2f, SD = %.2f", m_LCR_R1, sd_LCR_R1),
            x = 0.79, y = 0.71 + 0.12, hjust = 1, size = 10, color = "#ff0000") +
  draw_label(sprintf("Mean = %.2f, SD = %.2f", m_LCR_R2, sd_LCR_R2),
             x = 0.79, y = 0.71 + 0.07, hjust = 1, size = 10, color = "#0000ff")

p8_inset <- ggdraw(p8) +
  draw_image("/your/path/to/Landmark_Insets/PerifastigialSulcus.png", 
             x = 0.73, y = 0.58, width = 0.3, height = 0.3) +
  draw_label(sprintf("Mean = %.2f, SD = %.2f", m_BBS_R1, sd_BBS_R1),
            x = 0.79, y = 0.71 + 0.12, hjust = 1, size = 10, color = "#ff0000") +
  draw_label(sprintf("Mean = %.2f, SD = %.2f", m_BBS_R2, sd_BBS_R2),
             x = 0.79, y = 0.71 + 0.07, hjust = 1, size = 10, color = "#0000ff")

# Arrange all plots
final_plot <- plot_grid(p3_inset, p4_inset, 
                        p5_inset, p6_inset, p7_inset, p8_inset, nrow = 2)

final_plot
