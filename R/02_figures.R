# =============================================================================
#  02_figures.R
#  Figure 2 (forest: standalone vs dynamic borrowing) and
#  Figure 3 (borrowing is data-driven: tau and amount borrowed vs discordance).
#
#  Writes: figures/fig2_forest.png, figures/fig3_dynamic.png
#  Run from the repository root:  Rscript R/02_figures.R
# =============================================================================
suppressMessages({ library(bayesmeta); library(ggplot2); library(patchwork); library(dplyr) })
set.seed(20260705)
dir.create("figures", showWarnings = FALSE)

TEALD <- "#0C4F63"; TEAL <- "#1796A3"; SKY <- "#BFE3E6"
GREY  <- "#9AA0A6"; RUST <- "#A23B2E"; INK <- "#233238"
se_from_ci <- function(lo, hi) (log(hi) - log(lo)) / (2 * 1.959964)
tau_prior  <- function(scale) function(t) dhalfnormal(t, scale = scale)

endpoints <- list(
  DFS = list(y = c(A = log(0.61), B = log(0.51)),
             se = c(A = se_from_ci(0.34, 1.08), B = se_from_ci(0.29, 0.88))),
  Recurrence = list(y = c(A = log(0.49), B = log(0.42)),
             se = c(A = se_from_ci(0.24, 0.98), B = se_from_ci(0.21, 0.83)))
)

# ---- Figure 2: forest, standalone vs dynamic borrowing ----------------------
rows <- list(); mu <- list()
for (ep in names(endpoints)) {
  y <- endpoints[[ep]]$y; se <- endpoints[[ep]]$se
  fit <- bayesmeta(y = y, sigma = se, labels = c("A", "B"), tau.prior = tau_prior(0.5))
  mu[[ep]] <- exp(fit$summary["median", "mu"])
  for (g in c("A", "B")) {
    rows[[length(rows) + 1]] <- data.frame(endpoint = ep, stratum = g, type = "Standalone (ALASCCA)",
      hr = exp(y[g]), lo = exp(y[g] - 1.96 * se[g]), hi = exp(y[g] + 1.96 * se[g]))
    rows[[length(rows) + 1]] <- data.frame(endpoint = ep, stratum = g, type = "Dynamic borrowing",
      hr = exp(fit$theta["median", g]), lo = exp(fit$theta["95% lower", g]), hi = exp(fit$theta["95% upper", g]))
  }
}
df <- bind_rows(rows)
df$endpoint <- factor(df$endpoint, levels = c("DFS", "Recurrence"),
                      labels = c("Disease-free survival", "Recurrence (time to recurrence)"))
df$type <- factor(df$type, levels = c("Standalone (ALASCCA)", "Dynamic borrowing"))
df$ypos <- ifelse(df$stratum == "A", 2, 1) + ifelse(df$type == "Dynamic borrowing", 0.16, -0.16)
df$label <- sprintf("%.2f (%.2f-%.2f)", df$hr, df$lo, df$hi)
muband <- data.frame(endpoint = factor(levels(df$endpoint), levels = levels(df$endpoint)),
                     mu = unlist(mu))

fig2 <- ggplot(df, aes(hr, ypos, color = type)) +
  geom_vline(xintercept = 1, color = "grey55", linewidth = 0.4) +
  geom_vline(data = muband, aes(xintercept = mu), linetype = "22", color = TEALD, linewidth = 0.5) +
  geom_text(data = muband, aes(x = mu, y = 2.75, label = sprintf("pooled mu = %.2f", mu)),
            inherit.aes = FALSE, color = TEALD, size = 3.2, vjust = 0) +
  geom_segment(aes(x = lo, xend = hi, y = ypos, yend = ypos), linewidth = 1.0) +
  geom_point(aes(shape = type), size = 3.1, fill = "white", stroke = 1.1) +
  geom_text(aes(x = hi, label = label), hjust = -0.12, size = 2.75, show.legend = FALSE) +
  scale_color_manual(values = c("Standalone (ALASCCA)" = GREY, "Dynamic borrowing" = TEAL), name = NULL) +
  scale_shape_manual(values = c("Standalone (ALASCCA)" = 21, "Dynamic borrowing" = 19), name = NULL) +
  scale_x_continuous(trans = "log", breaks = c(0.2, 0.3, 0.5, 0.7, 1.0, 1.4), limits = c(0.2, 2.2)) +
  scale_y_continuous(breaks = c(1, 2),
                     labels = c("Group B\n(other PI3K variants)", "Group A\n(PIK3CA exon 9/20)"),
                     limits = c(0.6, 3.0)) +
  facet_wrap(~endpoint, ncol = 1) +
  labs(title = "Aspirin vs placebo: each stratum borrows strength from the other",
       subtitle = "Open circles = each stratum alone (ALASCCA). Filled circles = hierarchical dynamic borrowing.",
       x = "Hazard ratio (log scale) - values <1 favour aspirin", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold", color = TEALD, hjust = 0),
        plot.title = element_text(face = "bold", size = 12.5, color = TEALD),
        plot.subtitle = element_text(size = 9.3, color = INK),
        axis.text.y = element_text(color = INK, size = 9))
ggsave("figures/fig2_forest.png", fig2, width = 7.4, height = 6.0, dpi = 300, bg = "white")

# ---- Figure 3: borrowing is data-driven -------------------------------------
# Hold Group A fixed at HR 0.60; move Group B from concordant to opposite,
# at two precision levels; record tau and the fraction borrowed into A.
sweep_one <- function(seA, seB, label) {
  grid <- seq(log(0.45), log(1.6), length.out = 25)
  bind_rows(lapply(grid, function(lb) {
    fit <- bayesmeta(y = c(A = log(0.60), B = lb), sigma = c(A = seA, B = seB),
                     labels = c("A", "B"), tau.prior = tau_prior(0.5))
    data.frame(precision = label, hrB = exp(lb),
               tau = fit$summary["median", "tau"],
               borrow = 100 * max(1 - fit$theta["sd", "A"] / seA, 0))
  }))
}
dyn <- bind_rows(
  sweep_one(0.30, 0.30, "Imprecise data\n(ALASCCA-like, SE~=0.30)"),
  sweep_one(0.12, 0.12, "Precise data\n(large trials, SE~=0.12)"))
dyn_long <- bind_rows(
  transmute(dyn, precision, hrB, metric = "Amount borrowed into Group A (%)", value = borrow),
  transmute(dyn, precision, hrB, metric = "Between-stratum heterogeneity  tau", value = tau))
dyn_long$metric <- factor(dyn_long$metric,
  levels = c("Amount borrowed into Group A (%)", "Between-stratum heterogeneity  tau"))

fig3 <- ggplot(dyn_long, aes(hrB, value, color = precision)) +
  annotate("rect", xmin = 0.47, xmax = 0.55, ymin = -Inf, ymax = Inf, fill = SKY, alpha = 0.35) +
  annotate("text", x = 0.51, y = Inf, label = "ALASCCA\nGroup B", vjust = 1.3, size = 2.7,
           color = TEALD, lineheight = 0.9) +
  geom_vline(xintercept = 0.60, color = "grey70", linetype = "33") +
  annotate("text", x = 0.60, y = -Inf, label = "A", vjust = -0.6, size = 3, color = "grey40") +
  geom_line(linewidth = 1.05) +
  scale_color_manual(values = c("Imprecise data\n(ALASCCA-like, SE~=0.30)" = TEAL,
                                "Precise data\n(large trials, SE~=0.12)" = RUST), name = NULL) +
  scale_x_continuous(trans = "log", breaks = c(0.5, 0.6, 0.8, 1.0, 1.3, 1.6)) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  labs(title = "Borrowing is data-driven, not assumed",
       subtitle = "Group A fixed at HR 0.60; Group B moved from concordant (left) to opposite (right).",
       x = "Hazard ratio of Group B (discordance with Group A increases to the right)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", legend.text = element_text(size = 8.3),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", color = TEALD, hjust = 0),
        plot.title = element_text(face = "bold", size = 12.5, color = TEALD),
        plot.subtitle = element_text(size = 9.3, color = INK))
ggsave("figures/fig3_dynamic.png", fig3, width = 7.4, height = 6.2, dpi = 300, bg = "white")

cat("Done. Figures written to figures/ (fig2_forest.png, fig3_dynamic.png).\n")
