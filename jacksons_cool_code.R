
o2 <- read_table("https://tofu.byu.edu/docs/files/stat666/datasets/oliver2a")
o4 <- read_table("https://tofu.byu.edu/docs/files/stat666/datasets/oliver4a")

set.seed(39)

### EM Algorithm ###
num_steps = 50
em_o2 <- run_em(o2, max_steps = num_steps)
em_o2$mu
em_o2$sigma

plot(x = 1:em_o2$steps, y = em_o2$loglik, type = 'l')
plot(x = 1:em_o2$steps, y = em_o2$param_diffs[,1], type = 'l')
plot(x = 1:em_o2$steps, y = em_o2$param_diffs[,2], type = 'l')

# maximize liklihood, but stop at that step
final_step <- which.max(em_o2$loglik)
em_o2 <- run_em(o2, max_steps = final_step)
em_o2$df_imputed
# problem with imputed dataset is that it does not take into account random error, we need
# to create a distribution of values for each missing cell

### MVI Algorithm ###
num_sets <- 100 # number of datasets to create
mvi_o2 <- run_mvi(em_o2, num_sets = num_sets)
mvi_o2$df_dist # oof thats a mess wtf jackson

# check this out
num_miss <- nrow(mvi_o2$df_dist) # number of missing cells

par(mfrow = c(3, 3))
for (i in 1:num_miss) {
  imp_vals <- unlist(mvi_o2$df_dist[i, 1:num_sets])
  cell_ids <- (mvi_o2$df_dist[i, (num_sets+1) : (num_sets+2)])
  chosen_val <- mvi_o2$df_mvi[cell_ids$row[1], cell_ids$col[1]]
  hist(x = imp_vals, main = paste("row", cell_ids[1], 'col', cell_ids[2]))
  abline(v = chosen_val, col = "red", lwd = 3)
}
# NOTE: we created a normal distribution of each imputed value, then randomly selected
# from the distribution this imputed value.
# The red line indicates which value is in the final.
par(mfrow = c(1,1))

final_o2 <- mvi_o2$df_mvi

# if you want to change the final dataframe to any of the other 'num_sets', try this:
final_o2[em$which_na] <- unlist(mvi_o2$df_dist['set1'])


### Little's MCAR Test ###
# another option to determining mcar vs mar/mnar
naniar::mcar_test(o2)
naniar::mcar_test(o4)


### Testing if Linoleic and Linolenic acids are influential in distinguishing regions ###

# 🚧 under construction 🚧 #


