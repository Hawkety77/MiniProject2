library(readr)
o2 <- read_table("https://tofu.byu.edu/docs/files/stat666/datasets/oliver2a")
o4 <- read_table("https://tofu.byu.edu/docs/files/stat666/datasets/oliver4a")
set.seed(39)
o2_imp <- run_mvi(run_em(o2, max_steps = 7), num_sets = 1)$df_mvi
o4_imp <- run_mvi(run_em(o4, max_steps = 24), num_sets = 1)$df_mvi

o2_imp$region <- 2
o4_imp$region <- 4

ofull <- rbind(o2_imp, o4_imp)
ofull$region <- as.factor(ofull$region)

# COMPUTE Wilk's Lambdas #

fit_full <- manova(
  cbind(palmitic, palmitoleic, stearic, oleic, linoleic, eicosanoic, linolenic, eicosenoic) ~ region,
  data = ofull
  )
fit_red <- manova(
  cbind(palmitic, palmitoleic, stearic, oleic, eicosanoic, eicosenoic) ~ region,
  data = ofull
  )
lambda_yz <- summary(fit_full, test="Wilks")$stats["region", "Wilks"]
lambda_y  <- summary(fit_red,  test="Wilks")$stats["region", "Wilks"]
lambda_z_bar_y <- lambda_yz / lambda_y

# COMPUTE F-statistic

n_tmt = length(unique(ofull$region))
n_cols <- ncol(ofull)-1
q = ncol(fit_red$fitted.values)
p = n_cols - q
n = nrow(ofull)

# Test S
v_h <- n_tmt-1
v_e <- n - n_tmt - q # -q number of vars in y

t = find_t(q, v_h)
w = v_e + v_h - .5 * (p + v_h + 1)
df1 = p * v_h
df2 = w*t - .5 * (p * v_h - 2)

Fs <- wilk_to_f(lambda_z_bar_y, t, df1, df2)
paste("F-stat for region effect:", Fs)
paste("df1:", df1, "df2:", df2)
paste("P-value:", 1-pf(Fs, df1, df2))

# Simple write-up:
# because p-value < .05, we reject the null hypothesis and conclude that
# linoleic and linoeic acids ARE important in contributing to the significance
# and separation of the two regions. We did this using a conditional Wilk's
# lambda statistic, converted into an F statistic, whereby we determine that the
# contribution of the two acids are beyond the information already available
# from the other 6 acids.

# Note, something helpful for computing disc. scores

# since there are only two groups, can do something like...
xbar1 <- apply(o2_imp[,-9],2,mean)
xbar2 <- apply(o4_imp[,-9],2,mean)
S1 <- var(o2_imp[,-9])
S2 <- var(o4_imp[,-9])
n1 <- nrow(o2_imp)
n2 <- nrow(o4_imp)
Spl <- ((n1-1)*S1 + (n2-1)*S2)/(n1+n2-2)

a <- solve(Spl) %*% (xbar2 - xbar1)
a
astar <- sqrt(diag(Spl)) * a
astar

# simple write-up:
# in conjunction with what we found with our conditional wilk's lambda test statistic
# we find that the discriminant functions separating regions also indicate some form of
# separation for linoleic and linolenic acids. In this case, we see that the standardized
# discriminant scores for each of them have high positive influence, indicating that observations
# with higher values of linoleic and linolenic acids are indicative for belonging to Region 2





