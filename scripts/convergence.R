library(rhdf5)

source("scripts/theoretical.R")
source("scripts/utils.R")

# Generate a random genotype covariance matrix

#' Haar-distributed random orthogonal matrix
#'
#' @param dim Dimensions of the matrix
#'
#' @return Randomly generated orthogonal matrix with respect to Haar measure
rand_orthogonal <- function(dim) {
  qr.Q(qr(matrix(rnorm(dim^2), dim)))
}

rand_cov_G0 <- function(n_loci) {
  U_G0 <- rand_orthogonal(n_loci)
  U_G0 %*% diag(runif(n_loci)) %*% t(U_G0)
}

rand_cov_e <- function(n_traits) {
  U_e <- rand_orthogonal(n_traits)
  U_e %*% diag(runif(n_traits)) %*% t(U_e)
}

#' Generate a random trait effect size matrix
#'
#' @param n_loci Number of genetic loci
#' @param n_traits Number of traits
#'
#' @return A matrix of dimension (n_loci, n_traits)
rand_effects <- function(n_loci, n_traits) {
  idx_all <- seq_len(n_loci)

  effs_mat <- matrix(0, nrow = n_loci, ncol = n_traits)
  effs_size <- floor(n_loci / n_traits)

  for (j in seq_len(n_traits)) {
    effs_idx <- sample(idx_all, size = effs_size)
    effs_mat[effs_idx, j] <- rnorm(
      n = effs_size,
      sd = sqrt(1 / effs_size)
    )
  }

  effs_mat
}

rand_cor_m <- function(n_traits) {
  U_m <- rand_orthogonal(n_traits)
  V_m <- rand_orthogonal(n_traits)
  svals <- runif(n_traits) |> sort()

  U_m %*% diag(svals) %*% t(V_m)
}

#' Function for the fixed point evaluation of the genotype covariance recursion
#'
#' @param vech_G Half-vectorisation of the current genotype covariance
#' @param cov_G0 Initial (panmictic) genotype covariance matrix
#' @param cov_e Environmental component covariance matrix
#' @param cor_m Assorting matrix
#' @param lambda Effect size matria
#'
#' @return Evaluation of the fixed point function for nleqslv
fpe_function <- function(vech_G, cov_G0, cov_e, cor_m, lambda) {
  cov_G <- unvech(vech_G, nrow(cov_G0))
  cov_g <- t(lambda) %*% cov_G %*% lambda
  cov_X <- cov_g + cov_e
  mul_G <- cov_G %*% lambda

  cor_m_sym <- 0.5 * (cor_m + t(cor_m))

  vech((cov_G0 + mul_G %*% (cor_m_sym %c/c% cov_X) %*% t(mul_G)) - cov_G)
}

#' Heuristic check of convergence
#'
#' @param n_loci Number of genetic loci
#' @param n_traits Number of traits
#' @param n_checks Number of checks to perform
#'
#' @return list containing convergence data
fpe_check <- function(n_loci, n_traits, n_checks) {
  cov_G0 <- rand_cov_G0(n_loci)
  cov_e <- rand_cov_e(n_traits)
  cor_m <- rand_cor_m(n_traits)
  lambda <- rand_effects(n_loci, n_traits)

  cov_g <- t(lambda) %*% cov_G0 %*% lambda
  theta <-
    max(eigen(cov_g, symmetric = TRUE)$values) /
    min(eigen(cov_e, symmetric = TRUE)$values)

  cov_G_inf <- limiting_cov_G(cov_G0, cov_e, lambda, cor_m, 1e4)[[1e4]]

  criterion <- norm(cov_G0, type = "2") * norm(cor_m, type = "2") * theta

  solutions <- vector("list", n_checks)
  dist_inf <- vector("list", n_checks)
  n_converged <- 0

  for (j in seq_len(n_checks)) {
    sol <- nleqslv(
      vech(rand_cov_G0(n_loci)),
      fpe_function,
      cov_G0 = cov_G0,
      cov_e = cov_e,
      cor_m = cor_m,
      lambda = lambda
    )

    if (sol$termcd == 1) {
      n_converged <- n_converged + 1
      solutions[[j]] <- unvech(sol$x, n_loci)
    }
  }

  prop_converged <- n_converged / n_checks
  solutions <- solutions[seq_len(n_converged)]

  max_dist_pred <- 0
  for (j in seq_len(n_converged)) {
    dist <- norm(solutions[[j]] - cov_G_inf, type = "F")
    max_dist_pred <- max(max_dist_pred, dist)
  }

  max_dist <- 0
  for (j in seq_len(n_converged - 1)) {
    for (k in (j + 1):n_converged) {
      dist <- norm(solutions[[j]] - solutions[[k]], type = "F")
      max_dist <- max(max_dist, dist)
    }
  }

  list(
    cov_G0 = cov_G0,
    cov_e = cov_e,
    cor_m = cor_m,
    lambda = lambda,
    criterion = criterion,
    n_converged = n_converged,
    prop_converged = prop_converged,
    max_dist_pred = max_dist_pred,
    max_dist = max_dist
  )
}

#' Heuristic check of convergence across several configurations
#'
#' @param n_loci Number of genetic loci
#' @param n_traits Number of traits
#' @param n_setups Number of random setups to evaluate
#' @param n_checks Number of checks to perform
#'
#' @return data.frame containing convergence data across setups
fpe_checks <- function(n_loci, n_traits, n_setups, n_checks) {
  results <- parallel::mclapply(
    seq_len(n_setups),
    function(j) {
      check <- fpe_check(n_loci, n_traits, n_checks)

      c(
        param_id = j,
        criterion = check$criterion,
        prop_converged = check$prop_converged,
        max_dist_pred = check$max_dist_pred,
        max_dist = check$max_dist
      )
    },
    mc.cores = parallel::detectCores() - 1
  )

  results <- do.call(rbind, results) |> as.data.frame()

  results
}

convergence <- fpe_checks(25, 5, 1000, 50)
