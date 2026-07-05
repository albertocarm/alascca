# =============================================================================
#  03_figure_concept.R
#  Figure 1: Kruschke style diagram of the hierarchical borrowing model.
#  Reads top to bottom: priors on mu and tau, the shared distribution that the
#  two stratum effects are drawn from, and the likelihood that links them to the
#  observed log hazard ratios. Colours follow the JCO palette (ggsci).
#
#  Run from the repository root:  Rscript R/03_figure_concept.R
# =============================================================================
suppressMessages({ library(ggplot2); library(patchwork); library(ggsci) })
dir.create("figures", showWarnings = FALSE)

jco  <- pal_jco("default")(10)
BLUE <- jco[1]; FILL <- "#BcD9EC"; INK <- "#000000"
dhn  <- function(x, s) ifelse(x < 0, 0, sqrt(2/pi)/s * exp(-x^2/(2*s^2)))
arr  <- arrow(angle = 20, length = unit(0.26, "cm"), type = "closed")
th   <- theme_void() +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 10, colour = INK),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold", colour = INK),
        plot.margin = margin(2, 2, 2, 2))
theme_set(th)

# priors ----------------------------------------------------------------------
p_mu <- ggplot(data.frame(x = seq(-6, 6, .02)), aes(x, dnorm(x, 0, 3))) +
  geom_area(fill = FILL, colour = INK, linewidth = 0.4) +
  annotate("text", x = 0, y = 0.055, label = "vague", size = 3.2, colour = INK) +
  labs(title = "mu", subtitle = "mean effect")
p_tau <- ggplot(data.frame(x = seq(0, 2, .005)), aes(x, dhn(x, 0.5))) +
  geom_area(fill = FILL, colour = INK, linewidth = 0.4) +
  annotate("text", x = 0.95, y = 0.62, label = "Half-Normal(0.5)", size = 3.0, colour = INK) +
  labs(title = "tau", subtitle = "heterogeneity between strata")

# mu and tau feed the shared distribution ------------------------------------
p_conv <- ggplot() + xlim(0, 1) + ylim(0, 1) +
  annotate("segment", x = 0.25, xend = 0.5, y = 0.95, yend = 0.1, linewidth = 0.4, colour = INK) +
  annotate("segment", x = 0.75, xend = 0.5, y = 0.95, yend = 0.1, linewidth = 0.4, colour = INK) +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.1, yend = 0, linewidth = 0.5, colour = INK, arrow = arr) +
  annotate("text", x = 0.5, y = 0.5, label = "'~'", parse = TRUE, size = 7, colour = INK)

p_theta <- ggplot() + xlim(0, 1) + ylim(0, 1) +
  annotate("text", x = 0.5, y = 0.62, parse = TRUE, size = 5, fontface = "bold", colour = INK,
           label = "theta[A]*','~theta[B] ~ '~' ~ Normal(mu, tau^2)") +
  annotate("text", x = 0.5, y = 0.24, size = 3.4, colour = INK,
           label = "true effect of each stratum, shared through mu")

p_arrow <- ggplot() + xlim(0, 1) + ylim(0, 1.1) +
  annotate("segment", x = 0.5, xend = 0.5, y = 1, yend = 0, linewidth = 0.5, colour = INK, arrow = arr) +
  annotate("text", x = 0.5, y = 0.5, label = "likelihood", size = 3.2, vjust = -0.7, colour = INK)

# likelihood and data ---------------------------------------------------------
p_lik <- ggplot(data.frame(x = seq(-3, 3, .02)), aes(x, dnorm(x, -0.6, 0.3))) +
  geom_area(fill = FILL, colour = INK, linewidth = 0.4) +
  annotate("text", x = -0.6, y = 0.45, parse = TRUE, size = 3.0, colour = INK,
           label = "Normal(theta[j], sigma[j]^2)") +
  labs(subtitle = "observed log hazard ratio")
p_data <- ggplot() + xlim(0, 1) + ylim(0, 1) +
  annotate("text", x = 0.5, y = 0.5, parse = TRUE, size = 4.2, colour = INK,
           label = "y[A]*','~y[B]~': data (strata A and B)'")

fig1 <- ((p_mu | p_tau) / p_conv / p_theta / p_arrow / p_lik / p_data) +
  plot_layout(heights = c(2, 1.4, 0.8, 0.7, 2, 0.7))

ggsave("figures/figure_1_kruschke.pdf", fig1, width = 5.2, height = 7.2, device = cairo_pdf)
ggsave("figures/figure_1_kruschke.png", fig1, width = 5.2, height = 7.2, dpi = 300, bg = "white")
cat("Done. Figure 1 (Kruschke diagram) written to figures/.\n")
