## Retrospective Analysis of BASA

dir <- paste0(here::here("results/"), "lower.b0/sim_513/")
year <- 1
nyr <- 15
n.peels <- 5
annual.estimate <- matrix(NA, 43+nyr, n.peels)


for(y in (nyr-n.peels):15){
    est.ssb <- read_csv(paste0(dir, "/year_", y, "/model/mcmc_out/PFRBiomass.csv"), col_names=FALSE, show_col_types = FALSE) %>%
                  summarise(
                    across(
                      everything(),
                      median
                    )
                  )
    # assessment.biomass <- accumulate.assessment.posteriors(nyr.sim, total.sims, seeds, trial)
    # assessment.biomass.quants <- apply(assessment.biomass[[1]], 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
    # annual.assessment.biomass <- apply(assessment.biomass[[1]], c(2,3), median)
    est.ssb <- t(as.matrix(est.ssb))
    rownames(est.ssb) <- NULL
    for(j in 1:length(est.ssb)){
        annual.estimate[j, y] <- est.ssb[j, 1]
    }
}

final.year.est <- annual.estimate[44:nrow(annual.estimate), ncol(annual.estimate)]
peel.estimates <- apply(as.matrix(annual.estimate), 2, function(x) x[max(which(!is.na(x)))])
bias <- (peel.estimates - final.year.est)/final.year.est
mohn.rho <- mean(bias)
mohn.rho

plot(1980:2037, rep(NA, 58), ylim=c(0, 120000))
rownames(annual.estimate) <- 1980:2037
colors <- palette(rainbow(15))
for(p in 1:15){
    data <- annual.estimate[,p]
    lines(1980:2037, data, col=colors[p], lwd=3)
}
points(2023:2037, peel.estimates, col=colors, pch=19, lwd=10)