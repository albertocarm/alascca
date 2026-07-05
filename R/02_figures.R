# =============================================================================
#  02_figures.R
#  Data figures for the re-analysis, using the JCO colour palette (ggsci).
#  Figures are written as vector PDF (for the manuscript) and PNG (preview).
#  Explanatory text is kept in the manuscript captions, not inside the panels.
#
#    Figure 2  forest: standalone (ALASCCA) vs dynamic borrowing
#    Figure 3  adaptive behaviour: tau and amount borrowed vs discordance
#    Figure 4  quantile dot plot of the posterior probability of benefit
#
#  Run from the repository root:  Rscript R/02_figures.R
# =============================================================================
suppressMessages({ library(bayesmeta); library(ggplot2); library(ggsci); library(dplyr) })
set.seed(20260705)
dir.create("figures", showWarnings = FALSE)

jco  <- pal_jco("default")(10)
BLUE <- jco[1]   # #0073C2
GOLD <- jco[2]   # #EFC000
GREY <- jco[3]   # #868686
RED  <- jco[4]   # #CD534C
NAVY <- jco[6]   # #003C67
INK  <- "#000000"

se_from_ci <- function(lo, hi) (log(hi) - log(lo)) / (2 * 1.959964)
tau_prior  <- function(scale) function(t) dhalfnormal(t, scale = scale)
save_fig <- function(name, plot, w, h) {
  ggsave(file.path("figures", paste0(name, ".pdf")), plot, width = w, height = h, device = cairo_pdf)
  ggsave(file.path("figures", paste0(name, ".png")), plot, width = w, height = h, dpi = 300, bg = "white")
}

endpoints <- list(
  DFS = list(y = c(A = log(0.61), B = log(0.51)),
             se = c(A = se_from_ci(0.34, 1.08), B = se_from_ci(0.29, 0.88))),
  Recurrence = list(y = c(A = log(0.49), B = log(0.42)),
             se = c(A = se_from_ci(0.24, 0.98), B = se_from_ci(0.21, 0.83)))
)
fit_ep <- function(ep) bayesmeta(y = endpoints[[ep]]$y, sigma = endpoints[[ep]]$se,
                                 labels = c("A", "B"), tau.prior = tau_prior(0.5))

base_theme <- theme_bw(base_size = 11) +
  theme(text = element_text(colour = INK), axis.text = element_text(colour = INK),
        panel.grid.minor = element_blank(), legend.position = "top",
        legend.title = element_blank(), strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(colour = INK, face = "bold"))

# ---- Figure 2: forest, standalone vs dynamic borrowing ----------------------
rows <- list(); mu <- list()
for (ep in names(endpoints)) {
  fit <- fit_ep(ep); y <- endpoints[[ep]]$y; se <- endpoints[[ep]]$se
  mu[[ep]] <- exp(fit$summary["median", "mu"])
  for (g in c("A", "B")) {
    rows[[length(rows) + 1]] <- data.frame(endpoint = ep, stratum = g, type = "Standalone",
      hr = exp(y[g]), lo = exp(y[g] - 1.96 * se[g]), hi = exp(y[g] + 1.96 * se[g]))
    rows[[length(rows) + 1]] <- data.frame(endpoint = ep, stratum = g, type = "Dynamic borrowing",
      hr = exp(fit$theta["median", g]), lo = exp(fit$theta["95% lower", g]), hi = exp(fit$theta["95% upper", g]))
  }
}
df <- bind_rows(rows)
df$endpoint <- factor(df$endpoint, levels = c("DFS", "Recurrence"),
                      labels = c("Disease free survival", "Recurrence"))
df$type <- factor(df$type, levels = c("Standalone", "Dynamic borrowing"))
df$ypos <- ifelse(df$stratum == "A", 2, 1) + ifelse(df$type == "Dynamic borrowing", 0.17, -0.17)
df$label <- sprintf("%.2f (%.2f to %.2f)", df$hr, df$lo, df$hi)
muband <- data.frame(endpoint = factor(levels(df$endpoint), levels = levels(df$endpoint)), mu = unlist(mu))

fig2 <- ggplot(df, aes(hr, ypos, colour = type)) +
  geom_vline(xintercept = 1, colour = "grey50", linewidth = 0.4) +
  geom_vline(data = muband, aes(xintercept = mu), linetype = "22", colour = NAVY, linewidth = 0.5) +
  geom_segment(aes(x = lo, xend = hi, y = ypos, yend = ypos), linewidth = 0.9) +
  geom_point(aes(shape = type), size = 2.8, fill = "white", stroke = 1.0) +
  geom_text(aes(x = hi, label = label), hjust = -0.1, size = 2.7, colour = INK, show.legend = FALSE) +
  scale_colour_manual(values = c("Standalone" = GREY, "Dynamic borrowing" = BLUE)) +
  scale_shape_manual(values = c("Standalone" = 21, "Dynamic borrowing" = 19)) +
  scale_x_continuous(trans = "log", breaks = c(0.2, 0.3, 0.5, 0.7, 1.0, 1.4), limits = c(0.2, 2.4)) +
  scale_y_continuous(breaks = c(1, 2), labels = c("Group B", "Group A"), limits = c(0.65, 2.5)) +
  facet_wrap(~endpoint, ncol = 1) +
  labs(x = "Hazard ratio (log scale)", y = NULL) + base_theme
save_fig("figure_2_forest", fig2, 7.0, 5.2)

# ---- Figure 3: adaptive behaviour of the borrowing --------------------------
sweep_one <- function(seA, seB, label) {
  grid <- seq(log(0.45), log(1.6), length.out = 25)
  bind_rows(lapply(grid, function(lb) {
    fit <- bayesmeta(y = c(A = log(0.60), B = lb), sigma = c(A = seA, B = seB),
                     labels = c("A", "B"), tau.prior = tau_prior(0.5))
    data.frame(precision = label, hrB = exp(lb), tau = fit$summary["median", "tau"],
               borrow = 100 * max(1 - fit$theta["sd", "A"] / seA, 0))
  }))
}
dyn <- bind_rows(
  sweep_one(0.30, 0.30, "Imprecise data (SE 0.30)"),
  sweep_one(0.12, 0.12, "Precise data (SE 0.12)"))
dyn_long <- bind_rows(
  transmute(dyn, precision, hrB, panel = "Amount borrowed into Group A (%)", value = borrow),
  transmute(dyn, precision, hrB, panel = "Between stratum heterogeneity, tau", value = tau))
dyn_long$panel <- factor(dyn_long$panel,
  levels = c("Amount borrowed into Group A (%)", "Between stratum heterogeneity, tau"))

fig3 <- ggplot(dyn_long, aes(hrB, value, colour = precision)) +
  annotate("rect", xmin = 0.47, xmax = 0.55, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
  geom_vline(xintercept = 0.60, colour = "grey60", linetype = "33") +
  geom_line(linewidth = 1.0) +
  scale_colour_manual(values = c("Imprecise data (SE 0.30)" = BLUE, "Precise data (SE 0.12)" = RED)) +
  scale_x_continuous(trans = "log", breaks = c(0.5, 0.6, 0.8, 1.0, 1.3, 1.6)) +
  facet_wrap(~panel, scales = "free_y", ncol = 1) +
  labs(x = "Hazard ratio of Group B", y = NULL) + base_theme
save_fig("figure_3_adaptive", fig3, 7.0, 5.4)

# ---- Figure 4: quantile dot plot (each dot = 1% of posterior probability) ----
# Group A, disease free survival: analysed alone vs with dynamic borrowing.
fitD <- fit_ep("DFS")
mean_std <- endpoints$DFS$y["A"]; sd_std <- endpoints$DFS$se["A"]
mean_bor <- fitD$theta["median", "A"]
sd_bor   <- (fitD$theta["95% upper", "A"] - fitD$theta["95% lower", "A"]) / (2 * 1.959964)
qd <- function(m, s, lab) data.frame(model = lab, loghr = qnorm(ppoints(100), m, s))
dots <- bind_rows(qd(mean_std, sd_std, "Group A alone"),
                  qd(mean_bor, sd_bor, "Group A with dynamic borrowing"))
dots$model <- factor(dots$model, levels = c("Group A alone", "Group A with dynamic borrowing"))
dots$benefit <- ifelse(dots$loghr < 0, "Aspirin benefit (HR < 1)", "No benefit (HR >= 1)")
brks <- c(0.3, 0.4, 0.5, 0.7, 1.0, 1.4)

fig4 <- ggplot(dots, aes(x = loghr, fill = benefit)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.5) +
  geom_dotplot(method = "histodot", binwidth = 0.05, dotsize = 0.9, stackratio = 1.05, colour = NA) +
  scale_fill_manual(values = c("Aspirin benefit (HR < 1)" = BLUE, "No benefit (HR >= 1)" = RED)) +
  scale_x_continuous(breaks = log(brks), labels = brks) +
  facet_wrap(~model, ncol = 1) +
  labs(x = "Hazard ratio for disease free survival (log scale)", y = NULL) +
  base_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     panel.grid.major.y = element_blank())
save_fig("figure_4_quantile", fig4, 7.0, 4.8)

cat("Done. Figures 2 to 4 written to figures/ as PDF and PNG.\n")
