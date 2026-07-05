# ALASCCA reanalysis: Bayesian dynamic borrowing

Code for a hierarchical Bayesian reanalysis of the two molecular strata of the
ALASCCA trial (adjuvant aspirin in colorectal cancer with PI3K pathway
alterations; Martling et al., N Engl J Med 2025;393:1051-1064). It accompanies a
short methodological paper submitted to Clinical and Translational Oncology.

## The question

ALASCCA analysed its two prespecified strata separately: group A, the canonical
PIK3CA hotspot mutations in exon 9 or 20, and group B, other predicted activating
alterations in PIK3CA, PIK3R1 or PTEN. The two strata gave almost identical hazard
ratios, but because each was analysed on its own the group A disease free survival
estimate was not significant (0.61, 95% CI 0.34 to 1.08), while group B was.
Analysing the strata separately keeps biological noise down but loses statistical
power; pooling them into one analysis recovers the power but assumes the strata
behave identically.

## The reanalysis

Rather than choose, we let the strata share information according to how similar
they actually turn out to be. Each stratum effect is treated as a draw from a
common distribution, and a single heterogeneity parameter tau, estimated from the
data, controls how much they are pulled together:

    y_i     ~ Normal(theta_i, se_i^2)
    theta_i ~ Normal(mu, tau^2),   i = A, B
    tau     ~ Half-Normal(0.5)

When tau is small the strata share strength and shrink toward the common effect;
when it is large they stay apart. On the ALASCCA data tau is small (posterior
median around 0.24 for disease free survival), so group A borrows from the
concordant group B. Its disease free survival estimate moves to 0.58 (95% credible
interval 0.36 to 0.94), which no longer includes one, and the posterior probability
that aspirin helps rises from 95 to 99 percent. Recurrence is reinforced in the
same way. A sensitivity analysis included here shows the borrowing shutting itself
off when strata are made to diverge with precise data, so it cannot manufacture an
agreement that is not there.

Everything runs from the published summary hazard ratios, without individual
patient data, so the results can be reproduced from the trial report alone.

## Reproducing it

You need R 4.3 or later with a few packages:

    install.packages(c("bayesmeta", "ggplot2", "ggsci", "dplyr"))

Then, from the repository root:

    Rscript R/01_analysis.R          # estimates and tau-prior sensitivity -> results/
    Rscript R/02_figures.R           # forest, adaptive borrowing, quantile dot plot -> figures/
    Rscript R/03_figure_concept.R    # schematic of split / borrow / lump -> figures/

Tables are written to `results/` and figures to `figures/` as PDF and PNG. The
results were produced under R 4.3.2 with bayesmeta 3.5, ggplot2 4.0.0, ggsci 4.0.0
and dplyr 1.1.4.

## Source and license

Martling A, Hed Myrberg I, Nilbert M, et al. Low-dose aspirin for PI3K-altered
localized colorectal cancer. N Engl J Med. 2025;393(11):1051-1064.

Released under the MIT License (see `LICENSE`).
