# Gavin's Code for comparing region 2 to 4
library(tidyverse)
o2 <- read_table("https://tofu.byu.edu/docs/files/stat666/datasets/oliver2a")
o4 <- read_table("https://tofu.byu.edu/docs/files/stat666/datasets/oliver4a")
set.seed(39)
o2_imp <- run_mvi(run_em(o2, max_steps = 7), num_sets = 1)$df_mvi
o4_imp <- run_mvi(run_em(o4, max_steps = 24), num_sets = 1)$df_mvi

o2_imp$region <- 2
o4_imp$region <- 4

o2_without_region <- o2_imp %>% select(-region)
o4_without_region <- o4_imp %>% select(-region)


#library(Hotelling)
# Perform Hotelling's two-sample T-squared test
# hotelling_result <- hotelling.test(x = o2_without_region, y = o4_without_region)


# Hotelling's two-sample T^2 test from first principles

n1 <- nrow(o2_without_region)
n2 <- nrow(o4_without_region)
p <- ncol(o2_without_region)

# means
mean1 <- colMeans(o2_without_region)
mean2 <- colMeans(o4_without_region)
mean_diff <- mean1 - mean2

# covariances
cov1 <- cov(o2_without_region)
cov2 <- cov(o4_without_region)
pooled_cov <- ((n1 - 1) * cov1 + (n2 - 1) * cov2) / (n1 + n2 - 2)

# Hotelling's T^2 statistic
t_squared <- (n1 * n2) / (n1 + n2) * t(mean_diff) %*% solve(pooled_cov) %*% mean_diff

# Convert the T^2 to an F-statistic
f_statistic <- (n1 + n2 - p - 1) * t_squared / (p * (n1 + n2 - 2))
df1 <- p
df2 <- n1 + n2 - p - 1

# p-value
p_value <- 1 - pf(f_statistic, df1, df2)

cat("Hotelling's T^2 Statistic:", as.numeric(t_squared), "\n")
cat("F-statistic:", as.numeric(f_statistic), "\n")
cat("Numerator degrees of freedom:", df1, "\n")
cat("Denominator degrees of freedom:", df2, "\n")
cat("P-value:", as.numeric(p_value), "\n")



# possible write-up
# To compare whether Regions 2 and 4 are the same, we used Hotelling's 2-Sample T test. We calculated a T^2 value of 189.7534 which converts to an F-Statistic of 21.87435. 
# Using the F distribution with numerator degrees of freedom of 8 and a denominator degrees of freedom of 83 results in a P-Value of <.0001.
# This is statistically significant and provides overwhelming evidence that the overall fatty-acid profiles of regions 2 and 4 differ. 
# The agronomist's belief that region 2 and region 4 olives have evolved to have essentially the same profile in terms of the eight fatty acids is not supported by the data. 
# The observed differences are unlikely to have come from random sampling variation alone.
# It is important to note that the two samples were drawn by different organizations with potentially different data collection procedures, chemical analysis tools, 
# and data censoring mechanisms. This combined with missing data that was imputed using MVI exposes limitations of this analysis. 
# We do not know if significant variation comes from the different data collection procedures.






# ---------- Discriminant Function Calculation ----------

# Discriminant vector (a)
a <- solve(pooled_cov) %*% mean_diff

# Normalize for interpretability (optional)
a_std <- a / sqrt(sum(a^2))

# Print the coefficients
disc_fun <- data.frame(Variable = colnames(o2_without_region),
                       Coefficient = as.numeric(a_std))
print(disc_fun)

# ---------- Discriminant Scores for each observation ----------
combined_data <- rbind(o2_imp, o4_imp)

# Select only fatty acids (same order)
X <- as.matrix(combined_data[, colnames(o2_without_region)])

# Compute discriminant scores
combined_data$D_score <- X %*% a

# Check mean discriminant score by region
aggregate(D_score ~ region, data = combined_data, mean)

# ---------- Visualization ----------
library(ggplot2)

ggplot(combined_data, aes(x = D_score, fill = factor(region))) +
  geom_histogram(alpha = 0.6, bins = 20, position = "identity") +
  labs(title = "Discriminant Function Scores by Region",
       x = "Discriminant Score",
       y = "Count",
       fill = "Region") +
  theme_minimal()






