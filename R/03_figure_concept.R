# =============================================================================
#  03_figure_concept.R
#  Figure 1 (teaching schematic): split -> dynamic borrowing -> lump.
#  Illustrates the power/specificity trade-off and the hierarchical dial (tau)
#  using the ALASCCA Group-A disease-free-survival estimate.
#
#  Writes: figures/fig1_concept.png
#  Run from the repository root:  Rscript R/03_figure_concept.R
# =============================================================================
suppressMessages({ library(ggplot2) })
dir.create("figures", showWarnings = FALSE)

TEALD <- "#0C4F63"; TEAL <- "#1796A3"; GREY <- "#9AA0A6"
AMBER <- "#C6771E"; INK <- "#233238"; PAPER <- "#F4F9F9"

# Map a hazard ratio to a column's local x-range (log scale).
xmap <- function(hr, cx, hw = 12, lohr = 0.25, hihr = 1.7)
  cx - hw + (log(hr) - log(lohr)) / (log(hihr) - log(lohr)) * 2 * hw

cols <- data.frame(
  cx   = c(18, 50, 82),
  head = c("SPLIT", "DYNAMIC BORROWING", "LUMP"),
  sub  = c("analyse each stratum\nseparately (no pooling)",
           "hierarchical partial\npooling (Bayes)",
           "one combined analysis\n(full pooling)"),
  hcol = c(GREY, TEAL, AMBER),
  hr = c(0.61, 0.58, 0.56), lo = c(0.34, 0.36, 0.37), hi = c(1.08, 0.94, 0.83),
  note = c("Underpowered: Group A\nDFS 95% CI includes 1\n(HR 0.61, 0.34-1.08)",
           "Keeps power AND specificity:\nGroup A DFS now excludes 1\n(HR 0.58, 0.36-0.94)",
           "Maximum power, but assumes\nA = B; a true difference\nwould be erased (HR 0.56)"),
  stringsAsFactors = FALSE)

yF <- 34
p <- ggplot() + xlim(0, 100) + ylim(2, 60) + coord_cartesian(clip = "off") + theme_void()

# Top spectrum arrow.
p <- p +
  annotate("segment", x = 8, xend = 92, y = 54, yend = 54, linewidth = 1.1, color = INK,
           arrow = grid::arrow(ends = "both", length = unit(0.22, "cm"), type = "closed")) +
  annotate("text", x = 8, y = 56.4, label = "share NOTHING", hjust = 0, size = 3.1,
           color = GREY, fontface = "bold") +
  annotate("text", x = 92, y = 56.4, label = "share EVERYTHING", hjust = 1, size = 3.1,
           color = AMBER, fontface = "bold") +
  annotate("text", x = 50, y = 56.6, label = "tau  is learned from the data", hjust = 0.5,
           size = 3.0, color = TEAL, fontface = "italic") +
  annotate("segment", x = 50, xend = 50, y = 53.4, yend = 51.2, linewidth = 0.5, color = TEAL,
           arrow = grid::arrow(length = unit(0.15, "cm"), type = "closed"))

# Three strategy columns.
for (i in 1:3) {
  cx <- cols$cx[i]; x1 <- xmap(1, cx)
  p <- p +
    annotate("tile", x = cx, y = 48.5, width = 27, height = 4.6, fill = cols$hcol[i], alpha = 0.16) +
    annotate("text", x = cx, y = 49.4, label = cols$head[i], size = 3.7, fontface = "bold", color = cols$hcol[i]) +
    annotate("text", x = cx, y = 46.7, label = cols$sub[i], size = 2.6, color = INK, lineheight = 0.9) +
    annotate("segment", x = cx - 13, xend = cx + 13, y = yF - 4.4, yend = yF - 4.4, linewidth = 0.35, color = "grey75") +
    annotate("segment", x = x1, xend = x1, y = yF - 4.0, yend = yF + 3.2, linetype = "22", linewidth = 0.4, color = "grey55") +
    annotate("text", x = x1, y = yF - 5.6, label = "HR 1", size = 2.3, color = "grey45") +
    annotate("text", x = xmap(0.3, cx), y = yF - 5.6, label = "0.3", size = 2.3, color = "grey65") +
    annotate("segment", x = xmap(cols$lo[i], cx), xend = xmap(cols$hi[i], cx), y = yF, yend = yF,
             linewidth = 1.5, color = cols$hcol[i]) +
    annotate("point", x = xmap(cols$hr[i], cx), y = yF, size = 3.4, color = cols$hcol[i]) +
    annotate("text", x = cx, y = yF - 9.6, label = cols$note[i], size = 2.5, color = INK, lineheight = 0.95)
}
p <- p + annotate("text", x = 50, y = yF + 6.6,
  label = "What each strategy does to the Group A disease-free-survival estimate",
  size = 2.9, fontface = "italic", color = TEALD)

# Bottom hierarchical-model strip.
p <- p +
  annotate("tile", x = 50, y = 11, width = 92, height = 13.5, fill = PAPER, color = "grey80") +
  annotate("text", x = 50, y = 15.8, label = "The hierarchical model behind it",
           size = 3.3, fontface = "bold", color = TEALD) +
  annotate("text", x = 50, y = 12.4,
    label = "theta_A , theta_B  ~  Normal( mu , tau^2 )   -   each stratum's true log-HR is drawn around a common effect mu",
    size = 2.9, color = INK) +
  annotate("text", x = 50, y = 8.8,
    label = "tau -> 0 :  strata forced together (LUMP)        tau large :  strata kept apart (SPLIT)",
    size = 2.8, color = INK) +
  annotate("text", x = 50, y = 5.4,
    label = "Concordant strata (small tau) give strong borrowing;  discordant strata with firm data (large tau) switch borrowing off.",
    size = 2.7, fontface = "italic", color = TEAL)

ggsave("figures/fig1_concept.png", p, width = 8.6, height = 5.4, dpi = 300, bg = "white")
cat("Done. Figure written to figures/fig1_concept.png.\n")
