# =============================================================================
#  03_figure_concept.R
#  Figure 1: the three ways to analyse two molecular strata
#  (split, dynamic borrowing, lump). Kept deliberately spare; the interpretation
#  is given in the manuscript caption. Colours follow the JCO palette (ggsci).
#
#  Run from the repository root:  Rscript R/03_figure_concept.R
# =============================================================================
suppressMessages({ library(ggplot2); library(ggsci) })
dir.create("figures", showWarnings = FALSE)

jco  <- pal_jco("default")(10)
BLUE <- jco[1]; GREY <- jco[3]; DGOLD <- "#8F7700"; INK <- "#000000"
save_both <- function(name, plot, w, h) {
  ggsave(file.path("figures", paste0(name, ".pdf")), plot, width = w, height = h, device = cairo_pdf)
  ggsave(file.path("figures", paste0(name, ".png")), plot, width = w, height = h, dpi = 300, bg = "white")
}

# Map a hazard ratio to a column's local x range (log scale).
xmap <- function(hr, cx, hw = 11, lohr = 0.25, hihr = 1.7)
  cx - hw + (log(hr) - log(lohr)) / (log(hihr) - log(lohr)) * 2 * hw

cols <- data.frame(
  cx   = c(17, 50, 83),
  head = c("SPLIT", "DYNAMIC BORROWING", "LUMP"),
  sub  = c("each stratum on its own", "hierarchical partial pooling", "one combined estimate"),
  col  = c(GREY, BLUE, DGOLD),
  hr = c(0.61, 0.58, 0.56), lo = c(0.34, 0.36, 0.37), hi = c(1.08, 0.94, 0.83),
  tag  = c("wide interval, crosses 1", "narrower, excludes 1", "narrow, but assumes A = B"),
  stringsAsFactors = FALSE)

yF <- 30
p <- ggplot() + xlim(0, 100) + ylim(4, 56) + coord_cartesian(clip = "off") + theme_void() +
  theme(text = element_text(colour = INK))

# spectrum arrow
p <- p +
  annotate("segment", x = 10, xend = 90, y = 50, yend = 50, linewidth = 0.9, colour = INK,
           arrow = grid::arrow(ends = "both", length = unit(0.20, "cm"), type = "closed")) +
  annotate("text", x = 10, y = 52.2, label = "share nothing", hjust = 0, size = 3.0, colour = GREY) +
  annotate("text", x = 90, y = 52.2, label = "share everything", hjust = 1, size = 3.0, colour = DGOLD)

# three strategies with a mini forest of the Group A disease free survival estimate
for (i in 1:3) {
  cx <- cols$cx[i]; x1 <- xmap(1, cx)
  p <- p +
    annotate("text", x = cx, y = 44.5, label = cols$head[i], size = 3.8, fontface = "bold", colour = cols$col[i]) +
    annotate("text", x = cx, y = 41.6, label = cols$sub[i], size = 2.7, colour = INK) +
    annotate("segment", x = cx - 12, xend = cx + 12, y = yF - 4, yend = yF - 4, linewidth = 0.3, colour = "grey70") +
    annotate("segment", x = x1, xend = x1, y = yF - 3.6, yend = yF + 3.2, linetype = "22", linewidth = 0.35, colour = "grey55") +
    annotate("text", x = x1, y = yF - 5.4, label = "1", size = 2.5, colour = "grey45") +
    annotate("text", x = xmap(0.3, cx), y = yF - 5.4, label = "0.3", size = 2.5, colour = "grey65") +
    annotate("segment", x = xmap(cols$lo[i], cx), xend = xmap(cols$hi[i], cx), y = yF, yend = yF,
             linewidth = 1.3, colour = cols$col[i]) +
    annotate("point", x = xmap(cols$hr[i], cx), y = yF, size = 3.0, colour = cols$col[i]) +
    annotate("text", x = cx, y = yF - 8.4, label = cols$tag[i], size = 2.6, colour = INK)
}
p <- p +
  annotate("text", x = 50, y = 36.5, label = "Group A disease free survival estimate under each strategy",
           size = 2.9, fontface = "italic", colour = INK) +
  annotate("text", x = 50, y = 9.5,
           label = "theta_A, theta_B ~ Normal(mu, tau^2);  tau is estimated from the data and sets how much the strata share",
           size = 2.7, colour = INK)

save_both("figure_1_concept", p, 8.2, 4.6)
cat("Done. Figure 1 written to figures/ as PDF and PNG.\n")
