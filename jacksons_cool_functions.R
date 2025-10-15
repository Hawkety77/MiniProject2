
### EM Algorithm ###

run_em <- function(df, max_steps = 10, tol = 1e-4) {

  df <- as.matrix(df)
  n <- nrow(df)

  # initial estimates
  mu <- colMeans(df, na.rm = TRUE)
  sigma <- cov(df, use = "pairwise.complete.obs")
  param_diffs <- matrix(0, nrow = max_steps, ncol = 2)
  colnames(param_diffs) <- c("mu", "sigma")
  loglik_list <- numeric(max_steps)

  for (step in 1:max_steps) {
    df_filled <- df

    # --- E-step: impute missing values ---
    for (i in 1:n) {
      obs_idx <- which(!is.na(df[i, ]))
      mis_idx <- which(is.na(df[i, ]))
      if (length(mis_idx) > 0) {
        mu_obs <- mu[obs_idx]
        mu_mis <- mu[mis_idx]

        Sigma_oo <- sigma[obs_idx, obs_idx, drop = FALSE]
        Sigma_mo <- sigma[mis_idx, obs_idx, drop = FALSE]

        # regression coefficients
        B <- Sigma_mo %*% solve(Sigma_oo)

        x_obs <- df[i, obs_idx]
        df_filled[i, mis_idx] <- mu_mis + B %*% (x_obs - mu_obs)
      }
    }

    # --- M-step: update mu and sigma ---
    mu_new <- colMeans(df_filled)
    sigma_new <- cov(df_filled)

    # --- Compute observed-data log-likelihood ---
    loglik <- 0
    for (i in 1:n) {
      obs_idx <- which(!is.na(df[i, ]))
      x_obs <- df[i, obs_idx]
      mu_obs <- mu_new[obs_idx]
      Sigma_oo <- sigma_new[obs_idx, obs_idx, drop = FALSE]

      diff <- x_obs - mu_obs
      k <- length(obs_idx)

      term <- -0.5 * (
        log(det(Sigma_oo)) +
          t(diff) %*% solve(Sigma_oo) %*% diff +
          k * log(2 * pi)
      )
      loglik <- loglik + term
    }
    loglik_list[step] <- loglik

    # check convergence
    mu_diff <- max(abs(mu_new - mu))
    sigma_diff <- max(abs(sigma_new - sigma))
    param_diffs[step, ] <- c(mu_diff, sigma_diff)

    if (mu_diff < tol && sigma_diff < tol) {
      message("Converged at step ", step)
      break
    }

    mu <- mu_new
    sigma <- sigma_new
  }

  which_na = which(is.na(df), arr.ind = T)
  which_obs = which(!is.na(df), arr.ind = T)
  return(list(mu = mu, sigma = sigma,
              df_imputed = df_filled,
              steps = step,
              loglik = loglik_list[1:step],
              param_diffs = param_diffs[1:step,],
              which_na = data.frame(which_na),
              which_obs = data.frame(which_obs)
              ))

}

### MVI Algorithm ###

run_mvi <- function(em, num_sets = 10) {

  X_star <- em$df_imputed
  n <- nrow(X_star)
  d <- ncol(X_star)
  mu_star <- em$mu
  sigma_star <- em$sigma
  X_final <- X_star

  full_miss <- nrow(em$which_na)
  X_miss_dist <- data.frame(matrix(0, ncol = num_sets + 2, nrow = full_miss))
  X_miss_dist[,((num_sets+1):(num_sets+2))] <- em$which_na
  colnames(X_miss_dist) <- c(paste0("set", 1:num_sets), 'row', 'col')

  for (set in 1:num_sets) {
    for (i in 1:n) {
      obs_idx <- em$which_obs[,2][em$which_obs[,1] == i]
      mis_idx <- em$which_na[,2][em$which_na[,1] == i]
      n_miss <- length(mis_idx)
      if (n_miss > 0) {
        mu_obs <- mu_star[obs_idx]
        mu_mis <- mu_star[mis_idx]

        Sigma_oo <- sigma_star[obs_idx, obs_idx, drop = FALSE]
        Sigma_mo <- sigma_star[mis_idx, obs_idx, drop = FALSE]
        B <- Sigma_mo %*% solve(Sigma_oo)

        # compute error term
        miss <- 1:d %in% mis_idx
        E_Sigma <- sigma_star[miss,miss] - sigma_star[miss,!miss] %*% solve(sigma_star[!miss,!miss]) %*% sigma_star[!miss,miss]
        E <- MASS::mvrnorm(1, rep(0, n_miss), E_Sigma)

        # compute new imputation based on error term
        x_obs <- X_star[i, obs_idx]
        X_final[i, mis_idx] <- mu_mis + B %*% (x_obs - mu_obs) + E

        # set dist
        get_mat_rows <- X_miss_dist[X_miss_dist[num_sets+1] == i,]
        get_mat_rows[,set] <- X_final[i, mis_idx]
        X_miss_dist[X_miss_dist[num_sets+1] == i,set] <- get_mat_rows[,set]
      }
    }
  }



  full_error = X_final - X_star

  return(list(
    diff_mat = full_error,
    df_mvi = data.frame(X_final),
    df_dist = data.frame(X_miss_dist)
    ))

}

