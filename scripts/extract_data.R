if (!requireNamespace("rhdf5", quietly = TRUE)) {
  stop(
    "Package rhdf5 required to open results files. Install via BioConductor."
  )
}

read_params <- function(h5file) {
  var_gen <- drop(rhdf5::h5read(h5file, "/params/phenome/var_gen"))
  var_env <- drop(rhdf5::h5read(h5file, "/params/phenome/var_env"))
  gen_cor <- rhdf5::h5read(h5file, "/params/phenome/gen_cor")
  env_cor <- rhdf5::h5read(h5file, "/params/phenome/env_cor")

  sd_gen <- sqrt(var_gen)
  sd_env <- sqrt(var_env)

  list(
    cov_g0 = diag(sd_gen) %*% gen_cor %*% diag(sd_gen),
    cov_e = diag(sd_env) %*% env_cor %*% diag(sd_env),
    cor_m = rhdf5::h5read(h5file, "/params/mating/mate_cor"),
    var_g = var_gen,
    var_e = var_env,
    cor_g0 = gen_cor,
    cor_e = env_cor,
    n_gen = as.integer(rhdf5::h5read(h5file, "/params/simulation/n_gens")),
    n_ind = as.integer(rhdf5::h5read(h5file, "/params/simulation/n_ind")),
    trait_names = rhdf5::h5read(h5file, "/params/phenome/names")
  )
}

read_dataspace <- function(h5file, path) {
  data <- rhdf5::h5read(h5file, path)
  # attr <- rhdf5::h5read(h5file, path)
  # dimnames(data) <- list(attr$row_labels, attr$col_labels)

  data
}

read_scenario <- function(h5file) {
  params <- read_params(h5file)

  estimators <- list(
    h2 = "pheno_h2",
    cov_g = "pheno_genetic_genetic_cov",
    cov_x = "pheno_total_total_cov",
    cor_m = "mate_total_cor",
    cor_m_obj = "mate-objective",
    sib_cov = "sibling_genetic_cov",
    cousin_1_cov = "cousin_1_genetic_cov",
    cousin_2_cov = "cousin_2_genetic_cov",
    cousin_3_cov = "cousin_3_genetic_cov",
    ancestor_1_cov = "ancestor_1_genetic_cov",
    ancestor_2_cov = "ancestor_2_genetic_cov"
  )

  raw <- lapply(estimators, function(e) {
    read_dataspace(h5file, paste0("/raw/", e))
  })
  mean <- lapply(estimators, function(e) {
    read_dataspace(h5file, paste0("/summary/", e, "/mean"))
  })
  stderr <- lapply(estimators, function(e) {
    read_dataspace(h5file, paste0("/summary/", e, "/stderr"))
  })
  confint_lo <- lapply(estimators, function(e) {
    read_dataspace(h5file, paste0("/summary/", e, "/confint_lo"))
  })
  confint_hi <- lapply(estimators, function(e) {
    read_dataspace(h5file, paste0("/summary/", e, "/confint_hi"))
  })

  list(
    params = params,
    raw = raw,
    mean = mean,
    stderr = stderr,
    confint_lo = confint_lo,
    confint_hi = confint_hi
  )
}
