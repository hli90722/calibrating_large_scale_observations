# Check the vector beta has length 2
check_len_two <- function(beta) {
  stopifnot(length(beta) == 2)
}


# Returns a list of K indices dividing indices 1:n at random into
# K folds as equally sized as possible
get_fold_inds <- function(n, K, seed=NULL, type="foldid") {
  if (K == 1) {
    inds <- list()
    inds[[1]] <- 1:n
    if (type == "foldid") {
      return(rep(1, n))
    } else {
      return(inds)
    }
  }
  set.seed(seed)
  perm <- sample(n)
  cut_ind <- cut(1:n, breaks=K, labels=FALSE)
  if (type == "standard") {
    return(lapply(1:K, function(k) perm[which(cut_ind==k)]))
  } else if (type == "foldid") {
    return(cut_ind[perm])
  } else {
    stop("type must be one of 'standard' or 'foldid'")
  }
}

# Gets the index in vec that starts the first run of num consecutive non-decreases
# If no such index exists, returns length(vec)
get_early_stopping_index <- function(vec, num) {
  n <- length(vec)
  num <- ceiling(num)
  if (num > n) {
    stop("Need to pick num no larger than length(vec)")
  }
  if (n == 1) {
    return(1)
  }
  n_non_decreases <- 0
  for (i in 1:(n-1)) {
    curr_diff <- vec[(i+1)] - vec[i]
    if (curr_diff >= 0) {
      n_non_decreases <- n_non_decreases + 1
    } else {
      n_non_decreases <- 0
    }
    if (n_non_decreases == num) {
      return(i - num + 1)
    }
  }
  return(n)
}

# Adds a column of 1s to mat
add_intercept <- function(mat) {
  mat <- as.matrix(mat)
  return(cbind(rep(1, nrow(mat)), mat))
}

get_summands <- function(all_df, type, adj_hat, normalize) {
  if (!(all(c("s", "y", "mu_hat") %in% names(all_df)))) {
    stop("Not all of the following columns are in all_df: s, y, mu_hat")
  }
  if (type == "eff" && is.null(adj_hat)) {
    stop("adj_hat cannot be null with type = 'eff'")
  }
  
  ## Compute aipsw summands
  rho_hat <- mean(all_df$s)
  
  ## Compute collab summands
  if (normalize) {
    if (type != "aipsw") {
      stop("normalize=TRUE not supported for type != 'aipsw'")
    }
    N <- nrow(all_df)
    summands <-  N * with(all_df, s*(1-pi_hat)/pi_hat*(y-mu_hat)) / sum(with(all_df, s*(1-pi_hat)/pi_hat)) + with(all_df, (1-s)*mu_hat) / (1 - rho_hat)
  } else {
    if (type == "aipsw") {
      summands <- with(all_df,  s*(1-pi_hat)/pi_hat*(y-mu_hat) + (1-s)*mu_hat) / (1 - rho_hat)
    } else if (type == "aipsw_g") {
      summands <- with(all_df,  s*(1-g_hat)/g_hat*(y-mu_hat) + (1-s)*mu_hat) / (1 - rho_hat)
    } else if (type == "collab") {
      summands <-  with(all_df, s*(1-g_hat)/g_hat*(y-mu_hat) + (1-s)*mu_hat) / (1 - rho_hat) + with(all_df, (s-pi_hat) * (mu_hat-k_hat)) / (1 - rho_hat)
    } else if (type == "bar") {
      A_mat <- matrix(c(rho_hat, with(all_df, mean(s*Delta_hat)), with(all_df, mean(s*Delta_hat)), with(all_df, mean(s*Delta_hat^2))), nrow=2)
      b_vec <- c(1-rho_hat, with(all_df, mean((1-s)*Delta_hat)))
      coefs <- solve(A_mat) %*% b_vec
      all_df$alpha_hat <- coefs[1] + coefs[2] * all_df$Delta_hat
      summands <- with(all_df, s*alpha_hat*(y-mu_hat) + (1-s)*mu_hat) / (1 - rho_hat) + with(all_df, (s-pi_hat) * (mu_hat-k_hat)) / (1 - rho_hat)
    } else if (type == "eff") {
      adj <- (adj_hat * all_df$del_hat + all_df$s) * with(all_df, h_hat * (y-m_hat))
      summands <- with(all_df, s*(1-pi_hat)/pi_hat*(y-mu_hat) + (1-s)*mu_hat) / (1 - rho_hat) - adj
    } else if (type == "eff_g") {
      adj <- (adj_hat * all_df$del_hat + all_df$s) * with(all_df, h_hat * (y-m_hat))
      summands <- with(all_df,  s*(1-g_hat)/g_hat*(y-mu_hat) + (1-s)*mu_hat) / (1 - rho_hat) - adj
    } else if (type == "eff_collab") {
      adj <- (adj_hat * all_df$del_hat + all_df$s) * with(all_df, h_hat * (y-m_hat))
      summands <-  with(all_df, s*(1-g_hat)/g_hat*(y-mu_hat) + (1-s)*mu_hat) / (1 - rho_hat) + with(all_df, (s-pi_hat) * (mu_hat-k_hat)) / (1 - rho_hat) - adj
    } else if (type == "or") {
      summands <- with(all_df, (1-s)*mu_hat) / (1 - rho_hat)
    }
  }
  return(summands)
}

# Estimate influence function
estimate_if <- function(all_df, type, adj_hat) {
  summands <- get_summands(all_df, type, adj_hat, normalize=FALSE)
  rho_hat <- mean(all_df$s)
  if (type != "or") {
    adj <- 0
  } else {
    Mat <- matrix(c(rho_hat, with(all_df, mean(s*Delta_hat)),
                    with(all_df, mean(s*Delta_hat)), with(all_df, mean(s*Delta_hat^2))),
                  nrow=2)
    adj <- with(all_df, s*(y-mu_hat)) * as.numeric(with(all_df, cbind(1, Delta_hat) %*% solve(Mat) %*% c(1-rho_hat, mean((1-s)*Delta_hat)))) / (1-rho_hat)
  }
  return(summands + adj - (1-all_df$s) * mean(summands) / (1-rho_hat)) 
}

## Computes estimator and CI's
compute_estimator_and_ci <- function(all_df, type, adj_hat=NULL, normalize=FALSE, conf=0.95) {
  summands <- get_summands(all_df, type, adj_hat, normalize)
  
  ## Compute point estimate
  theta_hat <- mean(summands)
  
  ## Compute CI's
  N <- nrow(all_df)
  se <- sqrt(mean(estimate_if(all_df, type, adj_hat)^2) / N)
  ci_lower <- theta_hat - qnorm((1+conf)/2) * se
  ci_upper <- theta_hat + qnorm((1+conf)/2) * se
  
  return(list(theta_hat=theta_hat, ci_lower=ci_lower, ci_upper=ci_upper))
}
