# basa_run.r
# Created by John Trochta
# Date created:  06/08/2019
# Summary:
# This script runs the Bayesian ASA model using the No-U-turn (NUTS) MCMC sampler to obtain posteriors.
# The adnuts package developed by Cole Monnahan to run NUTS with ADMB (hence adnuts) is used.
# I highly recommend users should read the following for more details on NUTS application to stock assessments:
#   Cole C Monnahan, Trevor A Branch, James T Thorson, Ian J Stewart, Cody S Szuwalski, 
#         Overcoming long Bayesian run times in integrated fisheries stock assessments, 
#         ICES Journal of Marine Science, Volume 76, Issue 6, November-December 2019, 
#         Pages 1477–1488, https://doi-org.eres.qnl.qa/10.1093/icesjms/fsz059

# Some annotated guidance on diagnostic checks is provided below:
# 1) If divergent transitions>0, increase target acceptance rate to reduce step size
#    e.g. control = list(adapt_delta = 0.95)
# 2) IF extreme global correlations, pass dense matrix estimated from previous run
#    e.g. control-list(metric=M) where M is matrix in untransformed space
#    for ADMB models, use MLE covairance with control=list(metric="mle")

# With this adnuts, most important diagnostics are the:
#   1) ESS (accounts for autocorrelation)-500 ESS is sufficient for most quantities
#   2) Potential Scale reduction (R hat)-R hat fails if >1.1
#   3) No max tree depths exceeded (<12)
#   4) 0% divergences

# For divergence diagnoses and resolutions:
# https://discourse.mc-stan.org/t/divergent-transitions-a-primer/17099

#################################################################
library(data.table)
library(tidyverse)
library(adnuts)
library(snowfall)
library(rstan)
library(r4ss)

function.dir <- here::here("R/utils/")
source(file=paste0(function.dir, "fun_read_dat.R"))
source(file=paste0(function.dir, "calculate_ess.R"))
source(file=paste0(function.dir, "init_admb_params.R"))

#################################################################
OS <- "MAC"

# BE SURE TO CHECK YOUR DIRECTORY

run.basa.adnuts <- function(model.dir, seed, n.iter=350, n.warmup=100, max.duration=1000){

    setwd(model.dir)

    basa.root <- "~/Desktop/Projects/basa/model/"
    admb.model <- paste0(basa.root, "PWS_ASA")
    if(file.exists("PWS_ASA.tpl")){
      file.remove("PWS_ASA.tpl")
      file.remove("PWS_ASA(par).ctl")
    }
    file.symlink(paste0(admb.model, ".tpl"), ".")
    file.symlink(paste0(basa.root, "PWS_ASA(par).ctl"), ".")
    
    system("admb -s PWS_ASA", ignore.stdout = TRUE)
    

    PWS.ASA.data <- data_reader("PWS_ASA.dat")
    PWS.disease.data <- data_reader("PWS_ASA_disease.dat")

    nyr.fit <- PWS.ASA.data[[1]] 
    nyr.tot <- PWS.ASA.data[[2]] 
    nage <- 10                   # number of age classes (age 0 to 9+) (only care about 3-9+)

    # Read in measured age comps
    seine.age.comp <- PWS.ASA.data[[20]]
    spawn.age.comp <- PWS.ASA.data[[21]]
    vhsv.age.comp <- PWS.disease.data[[1]]
    ich.age.comp <- PWS.disease.data[[5]]

    # Read in the actual sample sizes
    # The following commands skips over lines of the file to read in tables 
    # so change the number skipped if file is modified  
    seine.samp.size <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=1,                 nrows=nyr.tot)
    spawn.samp.size <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=1+1*(nyr.tot),   nrows=nyr.tot)
    vhsv.samp.size  <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=1+2*(nyr.tot)+1,   nrows=nyr.tot)
    ich.samp.size   <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=1+3*(nyr.tot)+2, nrows=nyr.tot)


    # Create empty matrices to fill estimated ESS and age comps
    seine.ess.its <- matrix(0, nyr.tot, 1)  # Matrix to hold all iterations of the routine
    seine.ess.its <- seine.samp.size        # fill in the first column with the recorded sample size 

    spawn.ess.its <- matrix(0, nyr.tot, 1)
    spawn.ess.its <- spawn.samp.size

    vhsv.ess.its  <- matrix(0, nyr.tot, 1)
    vhsv.ess.its  <- vhsv.samp.size
    ich.ess.its   <- matrix(0, nyr.tot, 1)
    ich.ess.its   <- ich.samp.size

    ess.its <- listN(seine.ess.its, spawn.ess.its, vhsv.ess.its, ich.ess.its)

    # Change phases of the ESS in phases file to use the PWS_ASA(ESS_estimate)
    phases <- readLines("PWS_ASA(ESS).ctl", -1)
    phases[5] <- 1
    writeLines(phases, "PWS_ASA(ESS).ctl")

    seine.ess <- seine.samp.size
    spawn.ess <- spawn.samp.size
    vhsv.ess  <- vhsv.samp.size
    ich.ess   <- ich.samp.size

    age.comps <- list(seine=seine.age.comp, spawn=spawn.age.comp, vshv=vhsv.age.comp, ich=ich.age.comp)
    start.ess <- list(seine=seine.ess, spawn=spawn.ess, vhvs=vhsv.ess, ich=ich.ess)
    samp.size <- list(seine=seine.samp.size, spawn=spawn.samp.size, vhsv=vhsv.samp.size, ich=ich.samp.size)

    calc.ess <- calculate.ess(age.comps, start.ess, samp.size, ess.its, nyr.fit)

    # Turn of the phases for the ESS calculation so it no longer recalculates ESS in future model runs
    # Now write the converged ESS to a ctl file to be used for model runs
    write.table(
      rbind(
        "# PWS age comp effective sample sizes", paste0("# (", date(), ")"), " ",
        "# Determines which ctl file for the age comps and ESS to use (1 uses ESS control to be iteratively estimated)", -1, " ",
        "# Seine ESS",    calc.ess$seine.ess, " ", 
        "# Spawn ESS",    calc.ess$spawn.ess, " ", 
        "# Sero ESS",     calc.ess$vhsv.ess, " ", 
        "# Ich prev ESS", calc.ess$ich.ess
      ),
      file = "PWS_ASA(ESS).ctl", 
      append = F, 
      sep = " ",
      row.names=FALSE, 
      col.names=FALSE, 
      quote=F
    )
      
    ######################################################
    # Create reps x starting par vectors, and run NUTS
    setwd(model.dir)
    reps <- 4
    #set.seed(8558)
    set.seed(1120)
    seeds <- sample(1:1e4, size=reps)
    #system("admb -s PWS_ASA")
    system("./PWS_ASA -pinwrite -hbf 1", ignore.stdout = TRUE)

    inits <- init.admb.params(reps)

    # ADMB command for running nuts
    # PWS_ASA -nox -noest -nohess -maxfn 0 -nuts -mcmc 2000 -warmup 500 -chain 1 -mcseed 8682524 -max_treedepth 12 -adapt_delta 0.8 -adapt_mass -mcpin init.pin

    # Pilot run to check 

    if(dir.exists("mcmc_out")){
        system("rm mcmc_out/*.csv")
    }

    start.time <- Sys.time()
    fit.1 <- sample_nuts(model='./PWS_ASA',
                         path=model.dir,
                         init=inits, seeds=seeds, chains=reps, cores=reps,
                         iter=n.iter,
                         warmup=n.warmup,
                         duration=max.duration,
                         mceval=TRUE,
                         control=list(adapt_delta=0.9, metric="mle")
                    )
    end.time <- Sys.time()
    print(end.time - start.time)

    # Extracts NUTS stats (energy, leapfrog transitions,etc)
    mon <- monitor(fit.1$samples, warmup=fit.1$warmup, print=FALSE)
    x <- extract_sampler_params(fit.1)

    # Quick check for divergences & Gelman-Ruben statistic
    n.divergences <- sum(x$divergent__)/nrow(x)
    r.hat <- max(mon[, "Rhat"])<=1.1

    # Write summary of parameter posteriors (medians, percentiles, etc)
    write.csv(mon, file="mcmc_out/table_par_posterior_summary.csv")

    # Write all MCMC samples of the parameters
    mcmc.samps <- data.frame(matrix(fit.1$samples, nrow=dim(fit.1$samples)[3], byrow=TRUE))
    names(mcmc.samps) <- fit.1$par_names
    write.csv(mcmc.samps, file="mcmc_out/iterations.csv", row.names=FALSE)

    ## Examine the slowest mixing parameters
    # slow <- names(sort(mon[,"n_eff"]))[1:8]
    # pairs_admb(fit=fit, pars=slow)

    # Create summary file of NUTS/MCMC diagnostics
    sum.dia <- data.frame(divergences.from.extract.function=sum(x$divergent__)/nrow(x),
                          min.ESS=min(mon[, "n_eff"]),
                          which.min.ESS=names(which.min(mon[, "n_eff"])),
                          max.Rhat=max(mon[, "Rhat"]),
                          which.max.Rhat=names(which.max(mon[, "Rhat"])),
                          time.elapsed=end.time-start.time)
    write.table(sum.dia, file="mcmc_out/table_convergence_diagnostics.csv", sep=",", row.names=FALSE)
    saveRDS(fit.1, file="mcmc_out/NUTS_fit.RDS")

    return(
      list(
          divergences=n.divergences,
          converged=r.hat
      )
    )

}


