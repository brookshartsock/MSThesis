Tracking

The tracking data from Geant comes from .csvs from VLXe or hodoscope

For the linear regression algorithm it's linRegVLXe.cc or hodoscopeLinReg.cc

The histogram plots with error bars (it's beautiful isn't it?) is RStudio:
# Define the variables for x and y axes directly from the data
x_var <- data$energy  
y_var <- data$STDY    

# Define the number of bins and the range
num_bins <- 9
range_start <- 0
range_end <- 10

# Boolean variable to toggle log10 scale on y-axis
use_log_scale <- FALSE  # Set to FALSE to use linear scale

# Create bins for the x-variable
breaks <- seq(range_start, range_end, length.out = num_bins + 1)
bins <- cut(x_var, breaks = breaks, include.lowest = TRUE, right = FALSE)

# Calculate the mean and standard deviation of y-variable for each bin
m <- tapply(y_var, bins, mean)
sd <- tapply(y_var, bins, sd)

# Create a new data frame for plotting
df <- data.frame(mean.y = m, sd = sd, bin.start = breaks[-length(breaks)], bin.end = breaks[-1])

# Plot with shaded bins, title, and conditional y-axis scale
p <- ggplot(df) +
    geom_rect(aes(xmin = bin.start, xmax = bin.end, ymin = 0, ymax = mean.y), fill = "gray90", color = NA) +
    geom_segment(aes(x = bin.start, xend = bin.start, y = 0, yend = mean.y)) +
    geom_segment(aes(x = bin.end, xend = bin.end, y = 0, yend = mean.y)) +
    geom_errorbar(aes(x = (bin.start + bin.end) / 2, y = mean.y, 
                      ymin = mean.y - sd, 
                      ymax = mean.y + sd),
                  width = ((range_end-range_start)/num_bins)/2, color = "red") +
    geom_segment(aes(x = bin.start, xend = bin.end, y = mean.y, yend = mean.y)) +
    theme_minimal() +
    theme(panel.grid.major = element_line(color = "darkgray"), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white")) +
    xlab("Energy") +  
    ylab("Mean STDY") +
    ggtitle("Your Title Here")

# Apply log scale if the boolean variable is TRUE
if (use_log_scale) {
    p <- p + scale_y_log10()
}

# Print the plot
print(p)

To populate the lookup lists, first you need data from trackingMC or hodoTrackingMC
Use either populateVLXeLList.cc or populateHodoLList.cc respectively. you may want a new list with newList.cc

To track with the lookup list use lListVLXe.cc (you must have a list already) with data from VLXe

For DSDR tracking use DSDRLinRegVLXe.cc or hodoscopeLinReg.cc has an option to toggle this

2DHistgrams are made from 2DHistograms.cc

Analysis of algorithm porformance is mostly from RStudio though:
fdata = data[data$slopeY > -11 & data$slopeY < 15 & data$slopeZ > -11 & data$slopeZ < 15,]

normE = sqrt(1 + fdata$estSlopeY^2 + fdata$estSlopeZ^2)
normA = sqrt(1 + fdata$slopeY^2 + fdata$slopeZ^2)
dotP = acos((fdata$slopeY*fdata$estSlopeY)/(normA*normE) + (fdata$slopeZ*fdata$estSlopeZ)/(normA*normE) + 1/(normA*normE))

mean(dotP)
sd(dotP)

mean(fdata$estSlopeY-fdata$slopeY)
sd(fdata$estSlopeY-fdata$slopeY)

mean(fdata$estSlopeZ-fdata$slopeZ)
sd(fdata$estSlopeZ-fdata$slopeZ)
