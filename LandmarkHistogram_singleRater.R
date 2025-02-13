#-------------------------------------------------------------------------------------------------#
#----------------------- multiple plots with single rater ratings --------------------------------#
#-----------------------            LIDOtemplates                  -------------------------------#
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
ED_Rater1 <- landmarks_YY$Distance.export
nSamples <- nrow(ED_Rater1)

# Create data frames for each ROI
data_TBS <- data.frame(Rater = rep("Rater1", nSamples), Distance = ED_Rater1[,2])
data_OBSL <- data.frame(Rater = rep("Rater1", nSamples), Distance = ED_Rater1[,3])
data_OBSR <- data.frame(Rater = rep("Rater1", nSamples), Distance = ED_Rater1[,4])
data_LCL  <- data.frame(Rater = rep("Rater1", nSamples), Distance = ED_Rater1[,5])
data_LCR <- data.frame(Rater = rep("Rater1", nSamples), Distance = ED_Rater1[,6])
data_BBS <- data.frame(Rater = rep("Rater1", nSamples), Distance = ED_Rater1[,7])

# Define common theme settings and limits
common_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
  axis.text.x = element_text(size = 20, face = "bold"),
  axis.text.y = element_text(size = 20, face = "bold"),
  axis.title.x = element_text(size = 20, face = "bold"),
  axis.title.y = element_text(size = 20, face = "bold")
)
xlims <- c(-0.1, 4.5)
ylims <- c(0, 45)

# Create base plots
p3 <- ggplot(data_TBS, aes(x = Distance, fill = Rater)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, binwidth = 0.2) +
  scale_fill_manual(values = c("#2d89e0")) +
  common_theme + xlim(xlims) + ylim(ylims) +
  labs(title = "Periaqueductal Grey", x = "mm", y = "Number of scans") +
  theme(legend.position = "none")

p4 <- ggplot(data_OBSL, aes(x = Distance, fill = Rater)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, binwidth = 0.2) +
  scale_fill_manual(values = c("#2d89e0")) +
  common_theme + xlim(xlims) + ylim(ylims) +
  labs(title = "Outline Brainstem (L)", x = "mm", y = "Number of scans") +
  theme(legend.position = "none")

p5 <- ggplot(data_OBSR, aes(x = Distance, fill = Rater)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, binwidth = 0.2) +
  scale_fill_manual(values = c("#2d89e0")) +
  common_theme + xlim(xlims) + ylim(ylims) +
  labs(title = "Outline Brainstem (R)", x = "mm", y = "Number of scans") +
  theme(legend.position = c(0.9, 0.25), legend.direction = "vertical", 
        legend.background = element_blank(), legend.text = element_text(face = "bold"), 
        legend.title = element_text(face = "bold"))

p6 <- ggplot(data_LCL, aes(x = Distance, fill = Rater)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, binwidth = 0.2) +
  scale_fill_manual(values = c("#2d89e0")) +
  common_theme + xlim(xlims) + ylim(ylims) +
  labs(title = "4th Ventricle Border (L)", x = "mm", y = "Number of scans") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 2.5, color = "red", size = 1.2)

p7 <- ggplot(data_LCR, aes(x = Distance, fill = Rater)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, binwidth = 0.2) +
  scale_fill_manual(values = c("#2d89e0")) +
  common_theme + xlim(xlims) + ylim(ylims) +
  labs(title = "4th Ventricle Border (R)", x = "mm", y = "Number of scans") +
  theme(legend.position = c(0.9, 0.25), legend.direction = "vertical", 
        legend.background = element_blank(), legend.text = element_text(face = "bold"), 
        legend.title = element_text(face = "bold")) +
  geom_vline(xintercept = 2.5, color = "red", size = 1.2)

p8 <- ggplot(data_BBS, aes(x = Distance, fill = Rater)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, binwidth = 0.2) +
  scale_fill_manual(values = c("#2d89e0")) +
  common_theme + xlim(xlims) + ylim(ylims) +
  labs(title = "Perifastigial Sulcus", x = "mm", y = "Number of scans") +
  theme(legend.position = "none")

# Add inset images (approx. 1/3 of plot height) using cowplot's ggdraw & draw_image
# Adjust x and y positions if needed. Replace the file names with your actual image paths.
p3_inset <- ggdraw(p3) +
  draw_image("/your/path/to/Landmark_Insets/PeriaqueductalGrey.png", x = 0.73, y = 0.58, width = 0.3, height = 0.3)

p4_inset <- ggdraw(p4) +
  draw_image("/your/path/to/Landmark_Insets/OutlineBrainstem_left.png", x = 0.73, y = 0.58, width = 0.3, height = 0.3)

p5_inset <- ggdraw(p5) +
  draw_image("/your/path/to/Landmark_Insets/OutlineBrainstem_right.png", x = 0.73, y = 0.58, width = 0.3, height = 0.3)

p6_inset <- ggdraw(p6) +
  draw_image("/your/path/to/Landmark_Insets/LC_left.png", x = 0.73, y = 0.58, width = 0.3, height = 0.3)

p7_inset <- ggdraw(p7) +
  draw_image("/your/path/to/Landmark_Insets/LC_right.png", x = 0.73, y = 0.58, width = 0.3, height = 0.3)

p8_inset <- ggdraw(p8) +
  draw_image("/your/path/to/Landmark_Insets/PerifastigialSulcus.png", x = 0.73, y = 0.58, width = 0.3, height = 0.3)

# Arrange all plots
plot_grid(p3_inset, p4_inset, p5_inset, p8_inset, p6_inset, p7_inset, nrow = 2)
