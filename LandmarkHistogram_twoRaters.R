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

# Create data frames for each ROI (both raters)

# Periaqueductal Grey (Top BS)
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

# Outline Brainstem (L)
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

# Outline Brainstem (R)
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

# 4th Ventricle Border (L)
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

# 4th Ventricle Border (R)
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

# Perifastigial Sulcus (Bottom BS)
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

# Add inset images (using cowplot)
p3_inset <- ggdraw(p3) +
  draw_image("/Users/alex/Dropbox/paperwriting/coreg/figures/PeriaqueductalGrey.png", 
             x = 0.73, y = 0.58, width = 0.3, height = 0.3)
p4_inset <- ggdraw(p4) +
  draw_image("/Users/alex/Dropbox/paperwriting/coreg/figures/OutlineBrainstem_left.png", 
             x = 0.73, y = 0.58, width = 0.3, height = 0.3)
p5_inset <- ggdraw(p5) +
  draw_image("/Users/alex/Dropbox/paperwriting/coreg/figures/OutlineBrainstem_right.png", 
             x = 0.73, y = 0.58, width = 0.3, height = 0.3)
p6_inset <- ggdraw(p6) +
  draw_image("/Users/alex/Dropbox/paperwriting/coreg/figures/LC_left.png", 
             x = 0.73, y = 0.58, width = 0.3, height = 0.3)
p7_inset <- ggdraw(p7) +
  draw_image("/Users/alex/Dropbox/paperwriting/coreg/figures/LC_right.png", 
             x = 0.73, y = 0.58, width = 0.3, height = 0.3)
p8_inset <- ggdraw(p8) +
  draw_image("/Users/alex/Dropbox/paperwriting/coreg/figures/PerifastigialSulcus.png", 
             x = 0.73, y = 0.58, width = 0.3, height = 0.3)

# Arrange all plots
final_plot <- plot_grid(p3_inset, p4_inset, 
                        p5_inset, p6_inset, p7_inset, p8_inset, nrow = 2)

final_plot
