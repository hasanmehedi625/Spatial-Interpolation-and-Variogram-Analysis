# Variogram creation
rm(list = ls())
if (!require(pacman)) {
  install.packages("pacman")
}

pacman::p_load(ggplot2, tidyr, dplyr,scales,gridExtra)

pkgs <- c("ggplot2","dplyr","tidyr","scales","gridExtra")

# Set working directory
#setwd("D:/MRK_Documents/Teaching/Computational Geology/Lab exercises/VariogramKrigging")
print(paste("Current working directory:", getwd()))

# Step 1: Import/read your data
data <- read.table('SampleData.dat', header = FALSE)

# Step 2: Calculate Pairwise Distances
coords <- data[, 1:2]  # Extract coordinates
dist_matrix <- as.matrix(dist(coords))  # Euclidean distance matrix

# Step 3: Calculate Squared Differences for the variatble "elevation"
elvtn <- data[, 3]

# Pairwise squared differences
elvtn_matrix <- outer(elvtn, elvtn, FUN = function(x, y) 0.5 * (x - y)^2)
print(elvtn_matrix)

# Step 4: Extract the upper triangular part
dist_upper_half <- dist_matrix[upper.tri(dist_matrix)]
elvtn_upper_half <- elvtn_matrix[upper.tri(elvtn_matrix)]

plot(dist_upper_half, elvtn_upper_half)

# Step 5: Calculate Semi-Variance for each distance bin
bins <- seq(0, max(dist_upper_half), by = 10)  # Distance bins
midpoints <- (bins[-length(bins)] + bins[-1]) / 2  # Midpoints of bins

semi_variance <- sapply(1:(length(bins) - 1), function(i) {
  in_bin <- which(dist_upper_half >= bins[i] & dist_upper_half < bins[i + 1])
  if (length(in_bin) > 0) {
    mean(elvtn_upper_half[in_bin])
  } else {
    NA
  }
})

semi_variance <- semi_variance[!is.na(semi_variance)]
midpoints <- midpoints[!is.na(semi_variance)]


plot(midpoints, semi_variance)
# Step 6: Define the Model Functions
exp_model <- function(h, nugget, sill, range) {
  nugget + (sill - nugget) * (1 - exp(-h / range))
}

gaussian_model <- function(h, nugget, sill, range) {
  nugget + sill * (1 - exp(-(h / range)^2))
}

spherical_model <- function(h, nugget, sill, range) {
  ifelse(h <= range, 
         nugget + sill * (1.5 * (h / range) - 0.5 * (h / range)^3), 
         nugget + sill)
}

# Model Parameters
nugget <- 0.0
sill <- 1.7
hrange <- 50
model_parameters=data.frame(nugget=nugget,sill=sill,range=hrange)
h_values <- seq(0, max(bins), length.out = 100)

exp_variogram <- exp_model(h_values, nugget, sill, hrange)
gaussian_variogram <- gaussian_model(h_values, nugget, sill, hrange)
spherical_variogram <- spherical_model(h_values, nugget, sill, hrange)

# Step 7: Plotting
data_plot <- data.frame(X = data[, 1], Y = data[, 2], Z = data[, 3])
variogram_Cloud_plot <- data.frame(Distance = dist_upper_half, Variance = elvtn_upper_half)
exp_semi_var_plot <- data.frame(Distance = midpoints, Variance = semi_variance)
model_plot <- data.frame(Distance = h_values, 
                         Exponential = exp_variogram, 
                         Gaussian = gaussian_variogram, 
                         Spherical = spherical_variogram)

# Plotting
p1 <- ggplot(data_plot, aes(x = X, y = Y, color = Z)) +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  labs(title = "Data Points", color = "Z Value") +
  theme_minimal()


p2 <- ggplot() +
  geom_point(data = variogram_Cloud_plot, aes(x = Distance, y = Variance, color = "Variogram cloud"), 
             size = 2, alpha = 0.6) +
  geom_point(data = exp_semi_var_plot, aes(x = Distance, y = Variance, color = "Experimental Variogram"), 
             size = 2, alpha = 0.6) +
  geom_line(data = model_plot, aes(x = Distance, y = Exponential, color = "Exponential Model"), 
            linetype = "dashed", size = 1) +
  geom_line(data = model_plot, aes(x = Distance, y = Gaussian, color = "Gaussian Model"), 
            linetype = "dotted", size = 1) +
  geom_line(data = model_plot, aes(x = Distance, y = Spherical, color = "Spherical Model"), 
            size = 1) +
  scale_color_manual(values = c("Variogram cloud" = "grey",
                                "Experimental Variogram" = "red",
                                "Exponential Model" = "blue", 
                                "Gaussian Model" = "cyan", 
                                "Spherical Model" = "green")) +
  labs(title = "Model Fitting", x = "Lag Distance", y = "Semi-variance", 
       color = "Legend") +
  theme_minimal()

grid.arrange(p1, p2, nrow = 1)

write.csv(model_parameters, "model_parameters.csv", row.names = FALSE)
