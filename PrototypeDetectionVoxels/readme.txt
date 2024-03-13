Prototype Detection Voxels

Rolling average plot(s) was made in RStudio with:

x <- data$TIME
y <- data$CH1

window_size <- 11
n <- length(x)
avg_x <- sapply(1:(n - window_size + 1), function(i) mean(x[i:(i + window_size - 1)]))
avg_y <- sapply(1:(n - window_size + 1), function(i) mean(y[i:(i + window_size - 1)]))

plot(avg_x, avg_y, ylab = "Voltage (V)", xlab = "Time (s)", main = "Average", type = "l")

Landau fitting was preformed with landauFitPeaks.cc

Histogram fits are from histogramFits.cc

calibration curves came from calCurve.cc

Weighted means for timing were done on RStudio with:
data = rbind(dataNa, dataCs, dataMn)
x = data$timing
x = x[x != "NaN"]
mean(x)
sd(x)

x <- data$timing
y <- data$peak

finite_indices <- is.finite(x)
x <- x[finite_indices]
y <- y[finite_indices]

weighted_mean <- weighted.mean(x, w = y)
weighted_sd <- sqrt(sum(y * (x - weighted_mean)^2) / sum(y))
cat("Weighted Mean:", weighted_mean, "\n")
cat("Weighted Standard Deviation:", weighted_sd, "\n")

Coating histograms from the oscilloscope require some training data which I won't add here. Email me (brookshartsock@gmail.com) or Jon (folkertsjon@gmail.com) if you really want it. These where then fit using the same program as before

Coating comparison was done on RStudio with: (help from GPT4)
library(ggplot2)
library(tidyr)
library(dplyr)

data_long <- pivot_longer(data, 
                          cols = -Energy, 
                          names_to = c("Substance", "Type"), 
                          names_pattern = "([A-Za-z]+)(Peaks|SD)",
                          values_to = "Value")

# Separate peaks and SD into different columns
data_peaks <- data_long %>% 
  filter(Type == "Peaks") %>% 
  select(-Type) %>% 
  rename(Peak = Value)

data_sd <- data_long %>% 
  filter(Type == "SD") %>% 
  select(-Type) %>% 
  rename(SD = Value)

# Join Peaks and SD back into a single data frame
data_joined <- left_join(data_peaks, data_sd, by = c("Energy", "Substance"))

ggplot(data_joined, aes(x = Energy, y = Peak, color = Substance, group = Substance)) + 
  geom_point(aes(shape = Substance)) +  # Different shapes for each substance
  geom_errorbar(aes(ymin = Peak - SD, ymax = Peak + SD), width = 0.2) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Substance)) + # Different line types for each best fit line
  theme_minimal(base_size = 12) + 
  theme(
    panel.grid.major = element_line(color = "darkgray"), # Change grid to dark gray
    panel.grid.minor = element_line(color = "darkgray")  # Minor grid lines
  ) +
  labs(
    x = "Energy (MeV)", 
    y = "Peak Voltage (mV)", 
    title = "Coating Comparison",
    shape = "Substance",  # Legend for shapes
    linetype = "Substance" # Legend for line types
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c(1, 2, 3, 4, 5, 6)) + # Manual shapes
  scale_linetype_manual(values = c(1, 2, 3, 4, 5, 6))

Some of these files take peakIn.txt as an argument. Sometimes it's useful to run these in batches too with batch.bash
