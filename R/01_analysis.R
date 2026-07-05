# =============================================================================
#  01_analysis.R
#  Hierarchical Bayesian dynamic borrowing across the two molecular strata of
#  the ALASCCA trial (Martling et al., N Engl J Med 2025;393:1051-1064).
#
#  Inputs are the PUBLISHED summary-level subgroup hazard ratios and their
#  95% confidence intervals only; no individual patient data are used, so the
#  analysis is reproducible directly from the trial report.
#
#  Model (per endpoint):
#      y_i     ~ Normal(theta_i, se_i^2)     log-HR observed in stratum i
#      theta_i ~ Normal(mu, tau^2)           i = A, B  (the borrowing layer)
#      mu      ~ vague ;  tau ~ Half-Normal(scale)
#  tau is the data-driven amount of borrowing between strata.
#
#  Writes: results/results_table.csv, results/tau_sensitivity.csv
#  Run from the repository root:  Rscript R/01_analysis.R
# =============================================================================
suppressMessages({ library(bayesmeta); library(dplyr) })
set.seed(20260705)
dir.create("results", showWarnings = FALSE)

# Recover se(log-HR) from a published 95% CI on the HR scale.
se_from_ci <- function(lo, hi) (log(hi) - log(lo)) / (2 * 1.959964)
tau_prior  <- function(scale) function(t) dhalfnormal(t, scale = scale)
# P(HR < thr) from a Normal(mean, sd) posterior on the log-HR scale.
P_lt <- function(mean, sd, thr) pnorm(log(thr), mean, sd)

# Published ALASCCA subgroup results (Table used as summary-level input).
endpoints <- list(
  DFS = list(
    y  = c(A = log(0.61), B = log(0.51)),
    se = c(A = se_from_ci(0.34, 1.08), B = se_from_ci(0.29, 0.88))),
  Recurrence = list(
    y  = c(A = log(0.49), B = log(0.42)),
    se = c(A = se_from_ci(0.24, 0.98), B = se_from_ci(0.21, 0.83)))
)

# ---- 1. Standalone vs dynamic borrowing (tau ~ Half-Normal(0.5)) ------------
rows <- list(); tau_rows <- list()
for (ep in names(endpoints)) {
  y  <- endpoints[[ep]]$y
  se <- endpoints[[ep]]$se
  fit <- bayesmeta(y = y, sigma = se, labels = c("A", "B"),
                   tau.prior = tau_prior(0.5))

  tau_rows[[ep]] <- data.frame(
    endpoint = ep,
    tau_median = fit$summary["median", "tau"],
    tau_lower  = fit$summary["95% lower", "tau"],
    tau_upper  = fit$summary["95% upper", "tau"],
    # standard-error distance between the two observed log-HRs
    z_AB = abs(y["A"] - y["B"]) / sqrt(se["A"]^2 + se["B"]^2))

  mu_m <- fit$summary["median", "mu"]; mu_s <- fit$summary["sd", "mu"]
  rows[[length(rows) + 1]] <- data.frame(
    endpoint = ep, stratum = "Pooled (mu)", model = "Dynamic borrowing",
    hr = exp(mu_m), lo = exp(fit$summary["95% lower", "mu"]),
    hi = exp(fit$summary["95% upper", "mu"]),
    P_lt1 = P_lt(mu_m, mu_s, 1.0), P_lt060 = P_lt(mu_m, mu_s, 0.60))

  for (g in c("A", "B")) {
    # Standalone: each stratum alone (flat prior == frequentist estimate).
    rows[[length(rows) + 1]] <- data.frame(
      endpoint = ep, stratum = g, model = "Standalone",
      hr = exp(y[g]), lo = exp(y[g] - 1.959964 * se[g]),
      hi = exp(y[g] + 1.959964 * se[g]),
      P_lt1 = P_lt(y[g], se[g], 1.0), P_lt060 = P_lt(y[g], se[g], 0.60))
    # Dynamic borrowing: shrunken stratum estimate theta_g.
    th_m <- fit$theta["median", g]
    th_s <- (fit$theta["95% upper", g] - fit$theta["95% lower", g]) / (2 * 1.959964)
    rows[[length(rows) + 1]] <- data.frame(
      endpoint = ep, stratum = g, model = "Dynamic borrowing",
      hr = exp(fit$theta["median", g]), lo = exp(fit$theta["95% lower", g]),
      hi = exp(fit$theta["95% upper", g]),
      P_lt1 = P_lt(th_m, th_s, 1.0), P_lt060 = P_lt(th_m, th_s, 0.60))
  }
}
results <- bind_rows(rows) %>%
  mutate(across(c(hr, lo, hi, P_lt1, P_lt060), ~round(., 3))) %>%
  arrange(endpoint, factor(stratum, levels = c("A", "B", "Pooled (mu)")), model)
heterogeneity <- bind_rows(tau_rows)

write.csv(results, "results/results_table.csv", row.names = FALSE)
cat("== Table 1: standalone vs dynamic borrowing ==\n"); print(results, row.names = FALSE)
cat("\n== Between-stratum heterogeneity tau (Half-Normal(0.5) prior) ==\n")
print(heterogeneity, row.names = FALSE, digits = 3)

# ---- 2. Sensitivity of the pooled effect to the tau-prior scale -------------
sens <- list()
for (ep in names(endpoints)) for (sc in c(0.25, 0.50, 1.00)) {
  y <- endpoints[[ep]]$y; se <- endpoints[[ep]]$se
  fit <- bayesmeta(y = y, sigma = se, labels = c("A", "B"), tau.prior = tau_prior(sc))
  mu_m <- fit$summary["median", "mu"]; mu_s <- fit$summary["sd", "mu"]
  thA_m <- fit$theta["median", "A"]
  thA_s <- (fit$theta["95% upper", "A"] - fit$theta["95% lower", "A"]) / (2 * 1.959964)
  sens[[length(sens) + 1]] <- data.frame(
    endpoint = ep, tau_prior_scale = sc,
    tau_median = fit$summary["median", "tau"],
    mu_HR = exp(mu_m), mu_lo = exp(fit$summary["95% lower", "mu"]),
    mu_hi = exp(fit$summary["95% upper", "mu"]), P_mu_lt1 = P_lt(mu_m, mu_s, 1.0),
    A_borrowed_HR = exp(thA_m), A_P_lt1 = P_lt(thA_m, thA_s, 1.0))
}
sensitivity <- bind_rows(sens) %>% mutate(across(where(is.numeric), ~round(., 3)))
write.csv(sensitivity, "results/tau_sensitivity.csv", row.names = FALSE)
cat("\n== tau-prior sensitivity ==\n"); print(sensitivity, row.names = FALSE)
cat("\nDone. Tables written to results/.\n")
