source("scripts/extract_data.R")
source("scripts/theoretical.R")
source("scripts/utils.R")

library(ggplot2)
library(ggpalettes)
library(patchwork)

set_knoop_theme()

fig_h2_commuting <- function() {
  positive <- read_scenario("results/results_asm_positive_commuting.h5")
  negative <- read_scenario("results/results_asm_negative_commuting.h5")
  indefinite <- read_scenario("results/results_asm_indefinite_commuting.h5")

  make_h2_df <- function(ds, name) {
    h2_mean <- ds$mean$h2 |> drop()
    h2_stderr <- ds$stderr$h2 |> drop()
    h2_confint_lo <- ds$confint_lo$h2 |> drop()
    h2_confint_hi <- ds$confint_hi$h2 |> drop()

    lim_h2 <- lim_h2_commuting(
      ds$params$cov_g0,
      ds$params$cov_e,
      ds$params$cor_m
    )

    dimnames(h2_mean) <-
      dimnames(h2_stderr) <-
        dimnames(h2_confint_lo) <-
          dimnames(h2_confint_hi) <- list(paste0("trait", 1:5), 1:50)

    h2_mean <- data.frame(as.table(h2_mean))
    h2_stderr <- data.frame(as.table(h2_stderr))
    h2_confint_lo <- data.frame(as.table(h2_confint_lo))
    h2_confint_hi <- data.frame(as.table(h2_confint_hi))

    colnames(h2_mean) <- c("trait", "gen", "mean")
    colnames(h2_stderr) <- c("trait", "gen", "stderr")
    colnames(h2_confint_lo) <- c("trait", "gen", "confint_lo")
    colnames(h2_confint_hi) <- c("trait", "gen", "confint_hi")

    data.frame(
      name = name,
      trait = h2_mean$trait,
      gen = as.integer(h2_mean$gen),
      mean = h2_mean$mean,
      stderr = h2_stderr$stderr,
      confint_lo = h2_confint_lo$confint_lo,
      confint_hi = h2_confint_hi$confint_hi,
      h2_inf = lim_h2$h2_inf[as.integer(h2_mean$trait)] ,
      hg2_inf = lim_h2$hg2_inf[as.integer(h2_mean$trait)] 
    )
  }

  df_h2 <- do.call(
    rbind,
    Map(
      make_h2_df,
      list(positive, negative, indefinite),
      c("Positive", "Negative", "Indefinite")
    )
  )

  df_h2$name <- factor(
    df_h2$name,
    levels = c("Positive", "Indefinite", "Negative")
  )

  ggplot(
    df_h2,
    aes(
      x = gen,
      y = mean,
      ymin = confint_lo,
      ymax = confint_hi
    )
  ) +
    geom_hline(
      aes(
        yintercept = h2_inf,
        colour = trait
      ),
      linetype = "dashed",
      alpha = 0.5
    ) +
    geom_line(aes(colour = trait)) +
    geom_ribbon(aes(fill = trait), alpha = 0.25) +
    scale_colour_clinical() +
    scale_fill_clinical() +
    facet_wrap(~name) +
    labs(
      x = "Generation",
      y = "Heritability"
    ) +
    theme(
      legend.position = "none"
    )
}

fig_cov_g_noncommuting <- function() {
  positive <- read_scenario("results/results_asm_positive_noncommuting.h5")
  negative <- read_scenario("results/results_asm_negative_noncommuting.h5")
  indefinite <- read_scenario("results/results_asm_indefinite_noncommuting.h5")

  make_cov_df <- function(ds, name) {
    mean <- drop(ds$mean$cov_g)[,, 50]
    confint_lo <- drop(ds$confint_lo$cov_g)[,, 50]
    confint_hi <- drop(ds$confint_hi$cov_g)[,, 50]

    pred <- lim_cov_g(
      ds$params$cov_g0,
      ds$params$cov_e,
      ds$params$cor_m,
      1e4
    )[[1e4]]

    idx <- which(upper.tri(mean, diag = TRUE), arr.ind = TRUE)
    data.frame(
      name = name,
      type = ifelse(idx[, 1] == idx[, 2], "Variance", "Covariance"),
      predicted = pred[idx],
      sim_mean = mean[idx],
      sim_confint_lo = confint_lo[idx],
      sim_confint_hi = confint_hi[idx]
    )
  }

  df_cov <- do.call(
    rbind,
    Map(
      make_cov_df,
      list(positive, negative, indefinite),
      c("Positive", "Negative", "Indefinite")
    )
  )

  df_cov$name <- factor(
    df_cov$name,
    levels = c("Positive", "Indefinite", "Negative")
  )

  fig_identity <- ggplot(
    df_cov,
    aes(x = predicted, y = sim_mean, shape = type, colour = type)
  ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.25) +
    geom_point() +
    scale_colour_clinical() +
    facet_wrap(~name) +
    labs(
      x = "Predicted covariance",
      y = "Simulated covariance",
      colour = "Entry type",
      shape = "Entry type"
    )

  fig_residual <- ggplot(
    df_cov,
    aes(
      x = predicted,
      y = sim_mean - predicted,
      ymin = sim_confint_lo - predicted,
      ymax = sim_confint_hi - predicted,
      shape = type
    )
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.25) +
    geom_errorbar(
      aes(colour = type),
      alpha = 0.5,
      width = 0.00625,
      show.legend = FALSE
    ) +
    geom_point(aes(colour = type)) +
    scale_colour_clinical() +
    facet_wrap(~name) +
    scale_y_continuous(
      transform = scales::transform_pseudo_log(base = 10),
      labels = scales::label_log()
    ) +
    labs(
      x = "Predicted covariance",
      y = "Residual",
      colour = "Entry type",
      shape = "Entry type"
    ) +
    theme(
      strip.text.x = element_blank(),
      plot.margin = margin_part(t = 5.5)
    )

  (fig_identity / fig_residual) +
    plot_layout(guides = "collect", axes = "collect") &
    theme(legend.position = "bottom")
}

fig_kinship <- function() {
  positive <- read_scenario("results/results_asm_positive_noncommuting.h5")
  negative <- read_scenario("results/results_asm_negative_noncommuting.h5")
  indefinite <- read_scenario("results/results_asm_indefinite_noncommuting.h5")

  make_kinship_df <- function(ds, name) {
    n_gen <- ds$params$n_gen

    pred <- kinship_covariance(
      ds$params$cov_g0,
      ds$params$cov_e,
      ds$params$cor_m,
      1e4
    )

    mapping <- list(
      parent = "ancestor_1_cov",
      grandparent = "ancestor_2_cov",
      full_sibling = "sib_cov",
      first_cousin = "cousin_1_cov",
      second_cousin = "cousin_2_cov",
      third_cousin = "cousin_3_cov"
    )

    idx <- which(upper.tri(pred$parent, diag = TRUE), arr.ind = TRUE)

    do.call(
      rbind,
      Map(
        function(pred_name, sim_name) {
          sim_mean <- drop(ds$mean[[sim_name]])[,, n_gen]
          sim_lo <- drop(ds$confint_lo[[sim_name]])[,, n_gen]
          sim_hi <- drop(ds$confint_hi[[sim_name]])[,, n_gen]

          data.frame(
            name = name,
            relationship = pred_name,
            type = ifelse(idx[, 1] == idx[, 2], "Variance", "Covariance"),
            predicted = pred[[pred_name]][idx],
            sim_mean = sim_mean[idx],
            sim_confint_lo = sim_lo[idx],
            sim_confint_hi = sim_hi[idx]
          )
        },
        names(mapping),
        mapping
      )
    )
  }

  df_kinship <- do.call(
    rbind,
    Map(
      make_kinship_df,
      list(positive, negative, indefinite),
      c("Positive", "Negative", "Indefinite")
    )
  )

  df_kinship$relationship <- factor(
    df_kinship$relationship,
    levels = c(
      "parent",
      "grandparent",
      "full_sibling",
      "first_cousin",
      "second_cousin",
      "third_cousin"
    ),
    labels = c(
      "Parent",
      "Grandparent",
      "Full sibling",
      "First cousin",
      "Second cousin",
      "Third cousin"
    )
  )

  df_kinship$name <- factor(
    df_kinship$name,
    levels = c("Positive", "Indefinite", "Negative")
  )

  fig_identity <- ggplot(
    df_kinship,
    aes(
      x = predicted,
      y = sim_mean,
      ymin = sim_confint_lo,
      ymax = sim_confint_hi,
      colour = name,
      shape = name
    )
  ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.25) +
    geom_point() +
    facet_wrap(~relationship, scales = "free") +
    scale_colour_clinical() +
    labs(
      x = "Predicted covariance",
      y = "Simulated covariance",
      colour = "Assortment type",
      shape = "Assortment type"
    ) +
    theme(legend.position = "bottom")

  fig_residual <- ggplot(
    df_kinship,
    aes(
      x = predicted,
      y = sim_mean - predicted,
      ymin = sim_confint_lo - predicted,
      ymax = sim_confint_hi - predicted,
      colour = name,
      shape = name
    )
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.25) +
    geom_errorbar(width = 0, alpha = 0.25, show.legend = FALSE) +
    geom_point() +
    scale_colour_clinical() +
    scale_x_continuous(
      breaks = \(lims) c(0, signif(lims[2] * 0.8, 1)),
      labels = scales::label_number(drop0trailing = TRUE)
    ) +
    scale_y_continuous(
      transform = scales::pseudo_log_trans(base = 10),
      labels = scales::label_log()
    ) +
    labs(
      x = "Predicted covariance",
      y = "Residual",
      colour = "Assortment type",
      shape = "Assortment type"
    ) +
    facet_grid(name ~ relationship, scales = "free") +
    theme(
      strip.text.y = element_blank(),
      legend.position = "bottom",
      plot.margin = margin_part(t = 12)
    )

  (fig_identity / fig_residual) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

ggsave(
  fig_h2_commuting(),
  filename = "img/h2_commuting.pdf",
  width = 6.5,
  height = 2.5,
  device = cairo_pdf
)

ggsave(
  fig_cov_g_noncommuting(),
  filename = "img/cov_g_noncommuting.pdf",
  width = 6.5,
  height = 4.5,
  device = cairo_pdf
)

ggsave(
  fig_kinship(),
  filename = "img/kinship.pdf",
  width = 6.5,
  height = 8.5,
  device = cairo_pdf
)
