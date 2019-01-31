#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Simulation comparing the effect of treating each subject as their own control
## vs. having separate subjects as controls.

## ---- setup ----
library("plyr")
library("dplyr")
library("lme4")
library("reshape2")
library("ggplot2")
library("readr")

theme_set(theme_bw())
min_theme <- theme_update(
  panel.border = element_blank(),
  panel.grid = element_blank(),
  panel.spacing = unit(0, "line"),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  axis.text = element_text(size = 10),
  axis.title = element_text(size = 12),
  strip.background = element_blank(),
  strip.text = element_text(size = 10),
  legend.key = element_blank()
)

## ---- functions ----
#' Get Covariate Matrix for Simulations
#'
#' The X matrix in this simulation experiment includes a column of times, an
#' indicator of whether those times are in a prespecified window, and an
#' indicator for subject ID. The rows of this matrix index all combinations of
#' subjects x times.
#'
#' @param n_times [integer] The number of time points included in the
#'   simulation.
#' @param n_subjects [integer] The number of subjects included in the
#'   simulation.
#' @param window_size [integer] The number of time points to be considered for
#'   the "treatment effect" window.
#' @return x [data_frame] A data_frame indexing combinations of subjects and
#'   timepoints, with columns for time, subject, and a treatment window
#'   indicator.
simulate_covariates <- function(n_times, n_subjects, window_size) {
  times <- seq(-n_times / 2, n_times / 2)
  subjects <- seq_len(n_subjects)
  effect_window <- seq_len(window_size) - 1
  x <- as_data_frame(expand.grid(time = times, subject = subjects))
  x$window_indic <- 1 * (x$time %in% effect_window)
  x
}

#' Simulate a response vector, according to different designs
#'
#' To compare the estimation of treatment effects accoridng to different
#' experimental designs, we need to generate the actual observed y-values across
#' subjects in different conditions. This generates those y's, according to
#' either an "internal" subjects-as-their-own-controls design or an "external"
#' separate-subjects-as-treatments design. Here, the treatment effect is the
#' same across all subjects who are treated, but different subjects have a
#' different baseline.
#'
#' @param x [data.frame] Covariates
#' @param control_type [character] If "internal", all subjects are given
#'   treatment during the treatment window. Otherwise, half of patients are
#'   never given the treatment.
#' @param effect_size [numeric] The offset in the response in the treatment
#'   window, applied to all treatment subjects.
#' @param meas_sigma [numeric >= 0] The internal variability / noise a subject
#'   over time.
#' @param subject_sigma [numeric >= 0] The variation in baselines across
#'   subjects.
#' @return y [list] A list with two components
#'   $y [numeric vector] The actual measurement value
#'   $treat_subjects [binary vector] A binary indicator vector of whether the
#'   subject corresponding to the specified measurement index ever receives
#'   treatment.
simulate_response <- function(x, control_type, effect_size, meas_sigma, subject_sigma) {
  n_subjects <- length(unique(x$subject))
  y <- rnorm(n_subjects, 0, subject_sigma)[x$subject]

  if (control_type == "internal") {
    treat_subjects <- rep(TRUE, nrow(x))
  } else {
    treat_subjects <- x$subject %in% seq_len(n_subjects / 2)
  }

  cur_x <- x %>%
    filter(treat_subjects) %>%
    .[["window_indic"]]

  y[treat_subjects] <- y[treat_subjects] + effect_size * cur_x
  list(
    "treat_subjects" = 1 * treat_subjects,
    "y" = y + rnorm(length(y), 0, meas_sigma)
  )
}

## ---- simulate ----
## simulation parameters
n_times <- 20
window_size <- 5
meas_sigma <- 1
n_subjects <- 8

n_sim <- 10
effect_sizes <- round(seq(0, 2, length.out = 30), 3)
subject_sigmas <- round(seq(0, 5, length.out = 12), 3)

file.remove("x.tsv")
file.remove("y.tsv")
file.remove("sim.tsv")

## Generate x-values (subject IDs and times) and write to file
x <- simulate_covariates(n_times, n_subjects, window_size)
write.table(
  x,
  file = "x.tsv",
  sep = "\t",
  append = TRUE,
  col.names = FALSE
)

counter <- 1
n_total <- 4 * n_sim * length(effect_sizes) * length(subject_sigmas)
for (i in seq_len(n_sim)) {
  for (j in seq_along(effect_sizes)) {
    for (k in seq_along(subject_sigmas)) {
      for (control_type in c("internal", "external")) {
        for (inference_type in c("Mixed Effects", "Fixed Effect")) {

          if (counter %% 10 == 0) {
            cat(sprintf("simulating %s / %s \n", counter, n_total))
          }

          ## Simulate and write y to file, according to the experimental design
          y <- simulate_response(
            x,
            control_type,
            effect_sizes[j],
            meas_sigma,
            subject_sigmas[k]
          )

          write.table(
            cbind(counter, y$treat_subjects, y$y),
            file = "y.tsv", sep = "\t",
            append = TRUE,
            col.names = FALSE
          )

          ## Get t-statistics and write to file
          if (inference_type == "Mixed Effects") {
            fit <- lmer(y ~ window_indic + (1 | subject), data = cbind(x, y = y$y))
          } else {
            fit <- lm(
              y ~ window_indic,
              data = cbind(x, window_indic = x$window_indic * y$treat_subjects, y = y$y)
            )
          }

          t_value <- summary(fit)$coefficients[, "t value"]["window_indic"]
          cat(
            sprintf(
              "%d\t%d\t%f\t%f\t%s\t%s\t%f\n",
              counter, i, effect_sizes[j],
              subject_sigmas[k], control_type,
              inference_type, t_value
            ),
            file = "sim.tsv",
            append = TRUE
          )
          counter <- counter + 1
        }
      }
    }
  }
}

## ---- visualize-effects ----
sim <- read_tsv(
  "sim.tsv",
  col_names = c(
    "counter",
    "i",
    "effect_size",
    "subject_sigma",
    "control_type",
    "inference_type",
    "t_value"
  )
)

p <- ggplot(sim) +
  geom_point(
    aes(x = effect_size, y = t_value, col = control_type, shape = inference_type),
    alpha = 0.6, size = 0.75
  ) +
  scale_color_manual(values = c("#8ec269", "#7169c2")) +
  facet_wrap(~ subject_sigma, ncol = 4) +
  theme(
    strip.text.x = element_blank(),
    panel.border = element_rect(fill = "transparent", size = 1)
  ) +
  scale_x_continuous(breaks = c(0.5, 1.5)) +
  guides("col" = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  guides("shape" = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(
    x = "True Effect Size",
    y = "t-value",
    col = "Control Type",
    shape = "Inference"
  )

ggsave("S6-simulation_control_results.png", p, width = 5, height = 3)
ggsave("S6-simulation_control_results.eps", p, width = 5, height = 3)

## ---- visualize-inputs ----
x <- read_tsv(
  "x.tsv",
  col_names = c("ix", "time", "subject", "window_indic")
)
y <- read_tsv(
  "y.tsv",
  col_names = c("ix", "counter", "treat_subjects", "y")
)

combined <- x %>%
  left_join(y) %>%
  left_join(sim)

combined$window_indic <- ifelse(combined$window_indic, "ImmPost", "Pre")
combined$window_indic[combined$time > window_size] <- "Post"
combined$window_indic <- factor(
  combined$window_indic,
  levels = c("Pre", "ImmPost", "Post")
)

p <- ggplot(combined %>%
       filter(
         i == 1,
         effect_size %in% effect_sizes[c(TRUE, FALSE, FALSE, FALSE)],
         subject_sigma %in% subject_sigmas[c(TRUE, FALSE, FALSE)]
       )) +
  geom_line(
    aes(x = time, y = y, col = as.factor(window_indic), group = interaction(subject, counter)),
    size = 0.4, alpha = 0.7
  ) +
  scale_color_manual(values = c("#440154FF", "#DB4551", "#25858EFF")) +
  scale_x_continuous(breaks = c(-5, 0, 5)) +
  scale_y_continuous(breaks = c(-6, 0, 6)) +
  facet_grid(subject_sigma ~ effect_size) +
  labs(col = "Perturbation") +
  guides("col" = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme(
    axis.text = element_text(size = 8),
    panel.border = element_rect(fill = "transparent", size = 1)
  )
ggsave("S5-simulation_control_data_subset.png", p, width = 5.5, height = 3)
ggsave("S5-simulation_control_data_subset.eps", p, width = 5.5, height = 3)
