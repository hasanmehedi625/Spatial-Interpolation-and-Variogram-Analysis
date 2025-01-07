rm(list = ls())
if (!require(pacman)) {
  install.packages("pacman")
}

pacman::p_load(ggplot2, tidyr, dplyr,scales,gridExtra)

pkgs <- c("ggplot2","dplyr","tidyr","scales","gridExtra")
# Set working directory
#setwd("D:/MRK_Documents/Teaching/Computational Geology/Lab exercises/VariogramKrigging")
cat("Current working directory:", getwd(), "\n")

# Load the data
data <- read.table("SampleData.dat", header = FALSE)
x_known <- data[, 1]
y_known <- data[, 2]
z_known <- data[, 3]

# Define the grid
grid_x <- seq(0, 300, length.out = 200)
grid_y <- seq(0, 300, length.out = 200)
grid <- expand.grid(grid_x = grid_x, grid_y = grid_y)

# Variogram model (spherical model)
model_parameters=read.csv("model_parameters.csv", header = FALSE)
nugget <- 0.0
sill <- 1.7
hrange <- 50

variogram <- function(h, nugget, sill, range) {
  ifelse(h <= range,
         nugget + sill * (1.5 * (h / range) - 0.5 * (h / range)^3),
         nugget + sill)
}

# Initialize matrices for kriged values and variance
z_kriged <- matrix(0, nrow = length(grid_x), ncol = length(grid_y))
kriging_variance <- matrix(0, nrow = length(grid_x), ncol = length(grid_y))

# Kriging computations
for (i in seq_along(grid_x)) {
  for (j in seq_along(grid_y)) {
    # Current grid point
    x <- grid_x[i]
    y <- grid_y[j]
    
    # Compute distances to known points
    distances <- sqrt((x_known - x)^2 + (y_known - y)^2)
    
    # Variogram values for the current point
    gamma_k <- sapply(distances, variogram, nugget = nugget, sill = sill, range = hrange)
    
    # Variogram matrix for known points
    n <- length(x_known)
    gamma_matrix <- matrix(0, nrow = n + 1, ncol = n + 1)
    for (k in 1:n) {
      for (l in 1:n) {
        gamma_matrix[k, l] <- variogram(sqrt((x_known[k] - x_known[l])^2 + (y_known[k] - y_known[l])^2),
                                        nugget = nugget, sill = sill, range = hrange)
      }
    }
    gamma_matrix[n + 1, 1:n] <- 1
    gamma_matrix[1:n, n + 1] <- 1
    gamma_matrix[n + 1, n + 1] <- 0
    
    # Add Lagrange multiplier to gamma_k
    gamma_k <- c(gamma_k, 1)
    
    # Solve the kriging system
    weights <- solve(gamma_matrix, gamma_k)
    
    # Compute kriged estimate and variance
    z_kriged[i, j] <- sum(weights[1:n] * z_known)
    kriging_variance[i, j] <- sum(weights * gamma_k)
  }
}



# Kriged mean plot
kriged_data <- data.frame(
  x = rep(grid_x, each = length(grid_y)),
  y = rep(grid_y, times = length(grid_x)),
  z = as.vector(t(z_kriged))
)
variance_data <- data.frame(
  x = rep(grid_x, each = length(grid_y)),
  y = rep(grid_y, times = length(grid_x)),
  z = as.vector(t(kriging_variance))
)

p1=ggplot(kriged_data, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradientn(colours = terrain.colors(256), name = "Kriged Mean") +
  geom_point(data = data.frame(x = x_known, y = y_known),
                       aes(x = x, y = y), size = 2, shape = 21, fill = "black") +
  scale_color_gradientn(colours = terrain.colors(256), name = "Data Points") +
  labs(title = "Kriged Mean on Regular Grid", x = "X", y = "Y") +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )

# Kriged variance plot
p2=ggplot(variance_data, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradientn(colours = terrain.colors(256), name = "Kriged Variance") +
  geom_point(data = data.frame(x = x_known, y = y_known),
             aes(x = x, y = y), size = 2, shape = 21, fill = "black") +
  labs(title = "Kriged Variance on Regular Grid", x = "X", y = "Y") +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )
grid.arrange(p1, p2, nrow = 1)

# Save the plot with high resolution (300 dpi)
ggsave("krigged_mean.png", p1, dpi = 300, width = 10, height = 8, units = "in")
# Save the plot with high resolution (300 dpi)
ggsave("krigged_variance.png", p2, dpi = 300, width = 10, height = 8, units = "in")

