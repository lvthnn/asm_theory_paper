library(Matrix)
library(nleqslv)

source("scripts/utils.R")

#' Matrix operator montone exponentiation function
#'
#' @param mat The Hermitian matrix to raise to a power
#' @param pwr The power to raise the matrix to
#'
#' @return Matrix obtained by raising the spectrum of mat to the supplied power
matrix_power <- function(mat, pwr) {
  eig <- eigen(mat, symmetric = TRUE)
  eig$vectors %*% diag(eig$values^pwr) %*% t(eig$vectors)
}

#' Conjugate multiplication: mat2^{1/2} * mat1 * mat2^{1/2}
#'
#' @param mat1 Inner matrix
#' @param mat2 Outer matrix
#'
#' @return Matrix obtained from conjugate multiplication of mat1 by mat2
`%c*c%` <- function(mat1, mat2) {
  mat2_sqrt_inv <- matrix_power(mat2, 0.5)
  mat2_sqrt_inv %*% mat1 %*% mat2_sqrt_inv
}

#' Conjugate matrix division: mat2^{-1/2} * mat1 * mat2^{-1/2}
#'
#' @param mat1 Inner matrix
#' @param mat2 Outer matrix
#'
#' @return Matrix obtained from conjugate division of mat1 by mat2
`%c/c%` <- function(mat1, mat2) {
  mat2_sqrt <- matrix_power(mat2, -0.5)
  mat2_sqrt %*% mat1 %*% mat2_sqrt
}

#' Simulate a step from the genotype covariance recursion
#'
#' @param cov_G Current genotype covariance matrix
#' @param cov_G0 Initial (panmictic) genotype covariance matrix
#' @param cov_e Environmental component covariance matrix
#' @param lambda Effect matrix
#' @param cor_m Mate correlation (assorting) matrix)
#'
#' @return Forward iterate of the genotype covariance recursion
cov_G_step <- function(cov_G, cov_G0, cov_e, lambda, cor_m) {
  # Symmetrise the assorting matrix
  cor_m <- 0.5 * (cor_m + t(cor_m))

  # Compute genetic and total trait covariance, and regression matrix
  cov_g <- t(lambda) %*% cov_G %*% lambda
  cov_X <- cov_g + cov_e
  mul_G <- cov_G %*% lambda

  (0.5 * cov_G0) +
    (0.5 * cov_G) +
    (0.5 * mul_G %*% (cor_m %c/c% cov_X) %*% t(mul_G))
}

#' Simulate a genotype covariance trajectory under assortment
#'
#' @param cov_G0 Initial (panmictic) genotype covariance matrix
#' @param cov_e Environmental component covariance matrix
#' @param effects Effect matrix
#' @param cor_m Mate correlation (assorting) matrix)
#'
#' @return Forward iterate of the genotype covariance recursion
limiting_cov_G <- function(cov_G0, cov_e, effects, cor_m, n_gen) {
  vcov_G <- vector("list", n_gen)
  vcov_G[[1]] <- cov_G0

  for (t in (1 + seq_len(n_gen - 1))) {
    vcov_G[[t]] <- cov_G_step(
      vcov_G[[t - 1]],
      cov_G0,
      cov_e,
      effects,
      cor_m
    )
  }

  vcov_G
}

#' Simulate genetic component covariance trajectory
#'
#' @param cov_g Current genetic component covariance matrix
#' @param cov_g_0 Initial (panmictic) genetic component covariance matrix
#' @param cov_e Environmental component covariance matrix
#' @param cor_m Assorting matrix
#'
#' @return Forward iterate of the genetic component covariance recursion
cov_g_step <- function(cov_g, cov_g0, cov_e, cor_m) {
  # Symmetrise the assorting matrix
  cor_m <- 0.5 * (cor_m + t(cor_m))

  # Total trait covariance
  cov_X <- cov_g + cov_e

  (0.5 * cov_g0) +
    (0.5 * cov_g) +
    0.5 * (cov_g %*% (cor_m %c/c% cov_X) %*% cov_g)
}

#' Simulate the genetic component covariance recursion
#'
#' @param cov_g_0 The initial genetic component covariance matrix
#' @param cov_e The environmental component covariance matrix
#' @param cor_m The symmetrised mate correlation matrix
#' @param n_gen Number of generations to simulate
#'
#' @return A list of iterates from the genetic component covariance recursion
limiting_cov_g <- function(cov_g0, cov_e, cor_m, n_gen) {
  vcov_g <- vector("list", n_gen)
  vcov_g[[1]] <- cov_g0

  for (t in (1 + seq_len(n_gen - 1))) {
    vcov_g[[t]] <- cov_g_step(
      vcov_g[[t - 1]],
      cov_g0,
      cov_e,
      cor_m
    )
  }

  vcov_g
}

#' Commute the limiting heritability of a non-commuting system
#'
#' @param cov_g0 Initial genetic component covariance matrix
#' @param cov_e Environmental component covariance matrix
#' @param cor_m Assorting matrix
#' @param n_gen Number of generations simulated
#'
#' @return Limiting heritability matrix in non-commuting case
lim_h2_noncommuting <- function(
  cov_g0,
  cov_e,
  cor_m,
  n_gen
) {
  p <- nrow(cov_g0)

  # Compute mate correlation quantities
  eig_m <- eigen(cor_m)
  eigvec_m <- eig_m$vectors
  eigval_m <- eig_m$values
  sgn_m <- eigvec_m %*% diag(sign(eigval_m)) %*% t(eigvec_m)
  abs_m <- eigvec_m %*% diag(abs(eigval_m)) %*% t(eigvec_m)

  # Get the limiting genetic covariance matrix
  vcov_g <- limiting_cov_g(
    cov_g0,
    cov_e,
    cor_m,
    n_gen
  )

  cov_g_inf <- vcov_g[[n_gen]]
  cov_X_inf <- cov_g_inf + cov_e

  # Compute implicit term
  theta <- cov_g0 %c/c% cov_X_inf

  # Check feasibility in the indefinite case
  if (!isSymmetric(sgn_m %*% (theta %c*c% abs_m))) {
    stop("Discriminant not symmetric; can not evaluate closed form.")
  }

  I <- diag(p)

  hg2_inf <- 0.5 *
    sgn_m %*%
      (I -
        matrix_power(
          I - 4 * sgn_m %*% (theta %c*c% abs_m),
          0.5
        )) %c/c%
      abs_m

  cov_g_inf <- hg2_inf %c*c% cov_X_inf

  list(
    cov_g_inf = cov_g_inf,
    cov_X_inf = cov_X_inf,
    hg2_inf = hg2_inf,
    h2_inf = diag(cov_g_inf) / diag(cov_X_inf)
  )
}


#' Commute the limiting heritability of a non-commuting system using nleqslv
#'
#' @param cov_g0 Initial genetic component covariance matrix
#' @param cov_e Environmental component covariance matrix
#' @param cor_m Assorting matrix
#'
#' @return Limiting heritability matrix in non-commuting case
lim_h2_noncommuting_nleqslv <- function(
  cov_g0,
  cov_e,
  cor_m
) {
  # Symmetrise the assorting matrix
  cor_m <- 0.5 * (cor_m + t(cor_m))

  hg2_fpe <- function(vech_g, cov_g0, cov_e, cor_m) {
    cov_g <- unvech(vech_g, nrow(cov_g0))
    cov_X <- cov_g + cov_e
    hg2 <- cov_g %c/c% cov_X
    theta <- cov_g0 %c/c% cov_X

    vech(hg2 %*% cor_m %*% hg2 - hg2 + theta)
  }

  solution <- nleqslv(
    vech(cov_g0),
    hg2_fpe,
    cov_g0 = cov_g0,
    cov_e = cov_e,
    cor_m = cor_m
  )

  if (solution$termcd != 1) {
    warning("nleqslv did not converge; check solver output")
  }

  cov_g_inf <- unvech(solution$x, nrow(cov_g0))
  cov_X_inf <- cov_g_inf + cov_e
  hg2_inf <- cov_g_inf %c/c% cov_X_inf
  h2_inf <- diag(cov_g_inf) / diag(cov_X_inf)

  list(
    cov_g_inf = cov_g_inf,
    cov_X_inf = cov_X_inf,
    hg2_inf = hg2_inf,
    h2_inf = h2_inf
  )
}

#' Limiting per-trait heritability in the commuting case
#'
#' @param h2_0 Initial trait heritability
#' @param cor_m Eigenvalue of assorting matrix
#' @param tol Zero tolerance of assorting eigenvalue
#'
#' @return Limiting heritability of the supplied eigencomponent
lim_h2_scalar <- function(
  h2_0,
  cor_m,
  tol = 1e-12
) {
  ifelse(
    abs(cor_m) < tol,
    h2_0,
    (1 - sqrt(1 - 4 * cor_m * h2_0 * (1 - h2_0))) / (2 * cor_m * (1 - h2_0))
  )
}

#' Compute the limiting heritability of a commuting system
#'
#' @param cov_g0 Initial genetic component covariance
#' @param cov_e Environmental component covariance
#' @param cor_m Assorting matrix
#'
#' @return Limiting heritability matrix of the commuting system
lim_h2_commuting <- function(
  cov_g0,
  cov_e,
  cor_m
) {
  eig_g0 <- eigen(cov_g0, symmetric = TRUE)
  eigvec_g0 <- eig_g0$vectors

  # Assumes cov_g0, cov_e, and cor_m have the same eigenbasis
  eigval_g0 <- eig_g0$values
  eigval_e <- diag(t(eigvec_g0) %*% cov_e %*% eigvec_g0)
  eigval_m <- diag(t(eigvec_g0) %*% cor_m %*% eigvec_g0)

  hg2_0 <- eigval_g0 / (eigval_g0 + eigval_e)
  hg2_inf <- limiting_h2_scalar(hg2_0, eigval_m)

  eigval_xinf <- eigval_e / (1 - hg2_inf)
  eigval_ginf <- hg2_inf * eigval_xinf
  cov_ginf <- eigvec_g0 %*% diag(eigval_ginf) %*% t(eigvec_g0)

  list(
    hg2_inf = hg2_inf,
    h2_inf = diag(cov_ginf) / (diag(cov_ginf + cov_e))
  )
}

#' Compute the kinship matrix K
#'
#' @param cov_g_inf Equilibrium genetic component covariance matrix
#' @param cov_X_inf Equilibrium trait covariance matrix
#' @param cor_m Assorting matrix
#'
#' @return The kinship matrix
kinship_K_matrix <- function(cov_g_inf, cov_X_inf, cor_m) {
  0.5 *
    (diag(nrow(cov_g_inf)) +
      cov_g_inf %*% (cor_m %c/c% cov_X_inf))
}

#' Compute the equilibrium genetic component covariance between relatives
#'
#' @param cov_g0 Initial (panmictic) genetic component covariance
#' @param cov_e Environmental component covariance
#' @param cor_m Assorting matrix
#' @param n_gen Number of generations of forward iterations
kinship_cov_g <- function(cov_g0, cov_e, cor_m, n_gen) {
  cov_g_inf <- limiting_cov_g(
    cov_g0,
    cov_e,
    cor_m,
    n_gen
  )[[n_gen]]

  cov_X_inf <- cov_g_inf + cov_e

  K <- kinship_K_matrix(cov_g_inf, cov_X_inf, cor_m)
  K2 <- K %*% K
  K3 <- K %*% K2

  mrca_pair <- diag(nrow(cov_g0)) + cor_m %c/c% cov_X_inf
  mrca_pair <- 0.5 *
    (cov_g_inf +
      cov_g_inf %*% (cor_m %c/c% cov_X_inf) %*% cov_g_inf)

  list(
    parent = K %*% cov_g_inf,
    grandparent = K2 %*% cov_g_inf,
    half_sibling = K %*% cov_g_inf %*% t(K),
    full_sibling = mrca_pair,
    first_cousin = K %*% mrca_pair %*% t(K),
    second_cousin = K2 %*% mrca_pair %*% t(K2),
    third_cousin = K3 %*% mrca_pair %*% t(K3)
  )
}
