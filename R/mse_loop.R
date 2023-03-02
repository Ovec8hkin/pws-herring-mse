  # mse_loop.r
  # Created by John Trochta
  # Date created:  07/27/2020
  # Summary:
  # Runs the closed-loop simulation of an MSE for PWS herring (1 replicate)

source(file=paste0(here::here("R/utils/"),            "fun_read_dat.R"))
source(file=paste0(here::here("R/operating_model/"),  "fun_fish.R"))
source(file=paste0(here::here("R/operating_model/"),  "fun_operm.R"))
source(file=paste0(here::here("R/operating_model/"),  "fun_obsm.R"))
source(file=paste0(here::here("R/utils/"),            "fun_write_dat.R"))
source(file=paste0(here::here("R/estimation_model/"), "run_basa.R"))
source(file=paste0(here::here("R/utils/",             "hindcast.R")))

files.sources = list.files(here::here("R/operating_model/control_rules"), full.names=TRUE)
sapply(files.sources, source)

library(tidyverse)

#nyr.sim <- 10
nage <- 10

listN <- function(...){
  anonList <- list(...)
  names(anonList) <- as.character(substitute(list(...)))[-1]
  anonList
}

# Initialize variables in operating model (fun_operm.R)
initialize.popdyn.variables <- function(nyr.sim){

  rownames <- 1:nyr.sim
  colnames <- 0:(nage-1)

  survival.summer         <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  survival.winter         <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  maturity                <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  prefish.spawn.biomass   <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  seine.catch             <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  gillnet.catch           <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  pound.catch             <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  foodbait.catch          <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  n.spawners              <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  spawn.biomass.age.comp  <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  spawn.biomass           <- rep(0, nyr.sim)
  seine.biomass           <- rep(0, nyr.sim)
  true.nya                <- matrix(0, nyr.sim+1, nage, dimnames=list(c(rownames, nyr.sim+1), colnames))

  survey.indices          <- list(
      mdm                     <- rep(0, nyr.sim),
      egg                     <- rep(0, nyr.sim),
      PWSSC.hydro             <- rep(0, nyr.sim),
      ADFG.hydro              <- rep(0, nyr.sim),
      juv.schools             <- rep(0, nyr.sim),
      spawn.age.comp          <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames)),
      seine.age.comp          <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames)),
      vhsv.antibody           <- matrix(0, nyr.sim, nage*2, dimnames = list(rownames, rep(colnames, each=2))),
      ich.antibody            <- matrix(0, nyr.sim, nage*2, dimnames = list(rownames, rep(colnames, each=2)))  
  )
  names(survey.indices) <- c("mdm", "egg", "PWSSC.hydro", "ADFG.hydro", "juv.schools", "spawn.age.comp", "seine.age.comp", "vhsv.antibody", "ich.antibody")

  annual.age0.devs        <- rep(0, nyr.sim)

  pop_dyn <- listN(survival.summer, survival.winter, maturity, prefish.spawn.biomass, 
                seine.catch, gillnet.catch, pound.catch, foodbait.catch,
                n.spawners, spawn.biomass.age.comp, true.nya, 
                survey.indices,
                annual.age0.devs)

  return(pop_dyn)
  
}

set.initial.conditions <- function(dir, pop_dyn, mat, waa, sim.seed, start.year=2021){ 
  # Takes a random sample of the NYA distribution from year 0
  # (either then actual stock assessment or the hindcast model).
  set.seed(sim.seed)
  nyr <- start.year-1980+1
  init.nya <- read_csv(paste0(dir, "/year_0/model/mcmc_out/Num_at_age.csv"), col_names=FALSE, show_col_types = FALSE) %>%
                select_at(((10*nyr)-9):(10*nyr)) %>%
                slice_sample(n=1)

  pop_dyn$true.nya[1, ] <- as.numeric(init.nya) # Should these be rounded to integers?
  # pop_dyn$prefish.spawn.biomass[1, ] <- mat*as.numeric(init.nya)*waa
  
  return(pop_dyn)
}

# Reads in iterations.csv and takes a random samples
# of the parameter vector. Used for getting a random
# set of parameter estimates in order to have full
# level of uncertainty.
set.parameters <- function(dir, sim.seed){
  set.seed(sim.seed)
  iter.data <- read_csv(paste0(dir, "/year_0/model/mcmc_out/iterations.csv"), show_col_types=FALSE) %>% 
                na.omit() %>%
                select(!matches("[0-9]]$")) %>%
                slice_sample(n=1)
  return(as.list(iter.data))
}

# Computed median weight-at-age from raw WAA data
# over the years 2012-2022. WAA is assumed stable
# on an annual basis.
calculate.waa <- function(dat.files){
  waa.data <- dat.files$PWS_ASA.dat[[4]]
  true.waa <- apply(waa.data[31:41, ], 2, median)
  return(true.waa)
}

# Computed median fecundity from raw fecundity data
# over the years collected. Fecundity is assumed stable
# on an annual basis.
calculate.fec <- function(dat.files){
  fec.data <- dat.files$PWS_ASA.dat[[5]]
  true.fec <- apply(fec.data, 2, median, na.rm=TRUE)
  return(true.fec)
}

get.assessment.estimates <- function(dir, sim.year){
    est.ssb <- read_csv(paste0(dir, "/year_", sim.year, "/model/mcmc_out/PFRBiomass.csv"), col_names=FALSE, show_col_types = FALSE) %>%
                    select(last_col()) %>%
                    summarise(
                      across(
                        everything(),
                        quantile,
                        probs=c(0.025, 0.25, 0.5, 0.75, 0.975)
                      )
                    )

    est.nya <- read_csv(paste0(dir, "/year_", sim.year, "/model/mcmc_out/Num_at_age.csv"), col_names=FALSE, show_col_types = FALSE) %>% 
                select_at((ncol(.)-9):ncol(.)) %>% 
                summarise(
                  across(
                    everything(), 
                    quantile, 
                    probs=c(0.025, 0.25, 0.5, 0.75, 0.975)
                  )
                )

    return(listN(est.ssb, est.nya))
}

create.model.dir <- function(directory, year){
    model.dir <- paste0(directory, "/year_", year, "/model")
    if(dir.exists(model.dir)){
        #unlink(model.dir, recursive = TRUE, force=TRUE)
        unlink(model.dir, recursive=TRUE)
    }
    dir.create(model.dir, recursive = TRUE)
    return(paste0(model.dir, "/"))
}

compute.proportion.big.fish <- function(nya, big.fish.threshold=100){

    dist.means = c(NA, NA, NA, 63.7, 84.3, 120.0, 112.0, 127.0, 128.0, 142.0)
    dist.sds   = c(NA, NA, NA, 17.6, 22.2, 23.2,  18.2,  28.5,  37.5,  21.2)

    nya.integer = as.integer(nya)

    all.fish.weights <- c()

    for(i in 4:10){
        m = dist.means[i]
        s = dist.sds[i]
        n = nya.integer[i]
        all.fish.weights <- c(all.fish.weights, rnorm(n, m, s))
    }

    return(sum(all.fish.weights > big.fish.threshold)/length(all.fish.weights))

}

generate.recruitment.deviates <- function(nyr.sim, sim.seed){
  set.seed(sim.seed)
  max.regime.length <- 15
  
  devs <- rep(NA, nyr.sim)
  sigmas <- rep(NA, nyr.sim)
  
  devs[1:max.regime.length] <- rnorm(max.regime.length, 0.0645, 1.35)
  sigmas[1:max.regime.length] <- rep(1.20, max.regime.length)
  
  high.regime <- TRUE
  for(y in 1:(nyr.sim-max.regime.length)){
    if(y %% max.regime.length == 0) high.regime <- !high.regime
    
    dev <- ifelse(high.regime == 1, rnorm(1, -1.270, 1.05), rnorm(1, 0.0645, 1.35))
    sig <- ifelse(high.regime == 1, 1.05, 1.35)
    
    devs[y+max.regime.length] <- dev
    sigmas[y+max.regime.length] <- sig
    
  }
  
  return(list(devs=devs, sigmas=sigmas))
  
}

run.simulation <- function(hcr.options, nyr.sim, sim.seed=NA, write=NA, 
                           start.year=1, stop.year=NA, init.start.year=2021,
                           assessment=TRUE, hindcast=TRUE, cr.name="test"){
    # print(write)
    if(is.na(write)){
      write <- paste0(here::here("results"), "/test")
    }

    # Copy over current stock assessment to use for starting values
    model.0.dir <- create.model.dir(write, 0)
    if(is.na(hindcast) || !hindcast){
        file.remove(model.0.dir) # This is a hack to remove the internal model/ directory
        basa.root.dir <- "~/Desktop/Projects/basa/model/"
        file.symlink(basa.root.dir, paste0(write, "/year_0/"))
    }else{
        print("Running hindcast")
        hindcast(model.0.dir, sim.seed)
    }
    

    # Read in data files JUST ONCE, store then write within code
    dat.files <- read.data.files(paste0(write, "/year_", start.year-1, "/model/"))

    # These are the observed (not effective) sample sizes 
    sample.sizes <- list(
      spac = 1500,
      seac = 500
    )

    # Read these out of the mcmc_out/*.csv files
    true.waa <- calculate.waa(dat.files)
    true.fec <- calculate.fec(dat.files)
    
    # Assuming knife-edged selectivity at age-3
    # Could replace with maturity ogive?
    fish.selectivity <- matrix(1, nrow=4, ncol=10)
    fish.selectivity[, 1:3] <- 0 # Selectivity 0 for fish age 0-2

    log.mean.age0 <- read_csv(file.path(model.0.dir, "mcmc_out", "Num_at_age.csv"), col_names = FALSE, show_col_types = FALSE) %>%
      select(seq(1, 430, 10)) %>%
      rename_with(~ as.character(1980:2022), everything()) %>%
      summarise(
        across(
          as.character(1980:2022), 
          \(x) mean(x, na.rm=TRUE)
        )
      ) %>%
      as.matrix %>%
      as.vector %>%
      mean %>%
      log

    # Read in other parameters
    params <- read.par.file("~/Desktop/Projects/basa/model/PWS_ASA.par")
    par.samples <- set.parameters(write, sim.seed)
    params[names(par.samples)] <- par.samples
    params$log_MeanAge0     <- log.mean.age0
    params$female.spawners  <- tail(dat.files$PWS_ASA.dat$perc.female, 1)
    params$pk               <- dat.files$PWS_ASA.dat$pk
    params$waa              <- true.waa
    params$fec              <- true.fec
    params$selectivity      <- fish.selectivity
    params$catch.sd         <- 0.1

    maturity <- calc.maturity(params$mat_par_1, params$mat_par_2)

    # Initialize projection estimates (e.g. from forecast model)
    pop_dyn <- initialize.popdyn.variables(nyr.sim)
    pop_dyn <- set.initial.conditions(write, pop_dyn, maturity, params$waa, sim.seed, init.start.year-1)

    # Generate new age-0 recruitment deviates.
    # Devs are pulled from a regime-based recruitment function.
    recruitment <- generate.recruitment.deviates(nyr.sim, sim.seed)
    pop_dyn$annual.age0.devs <- recruitment$devs
    params$sigma_age0devs <- recruitment$sigmas
    #pop_dyn$annual.age0.devs <- rnorm(nyr.sim, mean=-0.35, sd=0.9) # change this so that the fishery doesn't immediately recover
    

    # Start loop
    control.rule <- rep(0, nyr.sim)

    # Need previous median assessment biomasses to compute HRs from gradient rule
    prev.ass.biomass <- read_csv(paste0(model.0.dir, "mcmc_out/PFRBiomass.csv"), col_names=FALSE, show_col_types = FALSE)%>% 
                          summarise( 
                            across( 
                              everything(), 
                              quantile, 
                              probs=c(0.025, 0.25, 0.5, 0.75, 0.975) 
                            ) 
                          ) 
    nyr <- ncol(prev.ass.biomass)
    ass.biomass <- data.frame(matrix(0, nrow=nyr+nyr.sim, ncol=6))
    colnames(ass.biomass) <- c("Year", "Biomass2.5", "Biomass25", "Biomass50", "Biomass75", "Biomass97.5")
    ass.biomass$Year <- 1980:(1980+nyr+nyr.sim-1)
    ass.biomass[1:nyr, 2:6] <- t(prev.ass.biomass)

    # Read data from previous, interrupted, model run into the appropiate pop_dyn variables
    # to facilitate restarting the run at any year necesarry. 
    if(start.year != 1){
        # This handles the population dynamics matrices
        tryCatch({
            vars <- names(pop_dyn)
            vars <- vars[! vars %in% c("survey.indices")]
            fnames <- apply(as.matrix(vars), 1, str_replace_all, pattern="[.]", replacement="_")
            file.paths <- apply(as.matrix(fnames), 1, function(f) paste0(write, "year_", start.year-1, "/results/", f, ".csv"))
            for(i in 1:length(vars)){
                data <- as.matrix(read.csv(file.paths[i]))
                dims <- dim(data)
                if(dims[2] > 2){
                    pop_dyn[[vars[i]]][1:dims[1],1:(dims[2]-1)] <- data[,2:ncol(data)]
                }else{
                    pop_dyn[[vars[i]]][1:dims[1]] <- data[,2:ncol(data)]
                }
            }

            # These handle the three special case matrics: annual.age0.devs, harvest_rate, and assessment_biomass
            age0.dev.data <- read.csv(paste0(write, "year_", start.year-1, "/results/annual_age0_devs.csv"))
            pop_dyn$annual.age0.devs[1:nrow(age0.dev.data)] <- age0.dev.data[,2]

            harvest.rate.data <- read.csv(paste0(write, "year_", start.year-1, "/results/harvest.csv"))
            control.rule[1:nrow(harvest.rate.data)] <- harvest.rate.data[,2]

            ass.biomass.data <- read.csv(paste0(write, "year_", start.year-1, "/results/assessment_biomass.csv"))
            ass.biomass[1:nrow(ass.biomass.data), 2:ncol(ass.biomass)] <- ass.biomass.data[,3:ncol(ass.biomass.data)]

        }, error=function(e){
            start.year <- 1
        })

    }

    if(is.na(stop.year)){
      stop.year <- nyr.sim
    }

    for(y in start.year:stop.year){  
        model.dir <- create.model.dir(write, y)
        setwd(model.dir)

        true.ssb <- sum(maturity*pop_dyn$true.nya[y, ]*params$waa)
        true.nya <- pop_dyn$true.nya[y, ]

        # If we want to run the full assessment, then use the assessment estimates in the HCR
        # calculation. Otherwise, use the deterministic biomass. This is useful for debugging
        # the operating model quickly (no need to run the whole assessment).
        hcr.name <- hcr.options$type

        if(assessment){
            ass.est <- get.assessment.estimates(write, y-1)
            ass.biomass[nyr+y, 2:6] <- t(ass.est$est.ssb)
            hcr.params <- list(
              curr.biomass = as.numeric(ass.est$est.ssb[3, ]),
              age.structure = as.numeric(ass.est$est.nya[3, ]),
              rel.biomass = as.numeric(ass.biomass[nyr+y, 4]/ass.biomass[nyr+y-3, 4]),
              weight.prop = compute.proportion.big.fish(as.numeric(ass.est$est.nya[3, ])),
              options = hcr.options
            )
        }else{
            hcr.params <- list(
              curr.biomass = true.ssb,
              age.structure = true.nya,
              rel.biomass = true.ssb/sum(pop_dyn$prefish.spawn.biomass[y-3, ]),
              options = hcr.options
            )
        }

        tar_hr <- do.call(hcr.name, c(hcr.params))

        control.rule[y] <- tar_hr

        # Execute fishery following management recommendation
        catch.at.age <- fun_fish(tar_hr, true.ssb, true.nya, params$waa, params$selectivity, params$catch.sd)

        # Run operating model
        pop_dyn <- fun_operm(y, pop_dyn$true.nya[y, ], catch.at.age, params, pop_dyn, sim.seed, project=TRUE)

        # Generate observations with error
        obs_w_err <- fun_obsm(pop_dyn$survey.indices, dat.files$PWS_ASA.dat$waa, dat.files$PWS_ASA.dat$fecundity, dat.files$PWS_ASA.dat$perc.female, 2.17, sample.sizes, y, sim.seed)

        wd <- getwd()
        
        # Write results so we can restart a failed run from specified year
        if(!is.na(write)){
            setwd(write)
            results.dir <- paste0(write, "/", "year_", y, "/results/")
            if(dir.exists(results.dir)){
              unlink(results.dir, recursive = TRUE)
            }
            dir.create(results.dir, recursive = TRUE)
            setwd(results.dir)

            files <- apply(as.matrix(names(pop_dyn)), 1, str_replace_all, pattern="[.]", replacement="_")
            lapply(seq_along(pop_dyn), function(i){
              write.csv(pop_dyn[[i]], paste0(files[i], ".csv"), row.names = TRUE)
            })

            files <- apply(as.matrix(names(pop_dyn$survey.indices)), 1, str_replace_all, pattern="[.]", replacement="_")
            lapply(seq_along(pop_dyn$survey.indices), function(i){
              write.csv(pop_dyn$survey.indices[[i]], paste0(files[i], ".csv"), row.names = TRUE)
            })

            write.csv(control.rule, "harvest.csv")
            write.csv(ass.biomass, "assessment_biomass.csv")

        }
        setwd(wd)

        # Write all new data to .dat file (and modify covariate and agecomp_samp_sizes.txt files)
        fun_write_dat(dat.files, catch.at.age, obs_w_err, params, sample.sizes, y)


        # Run BASA
        #  - allow option to skip years of assessment
        if(assessment){

            # Write all new data to .dat file (and modify covariate and agecomp_samp_sizes.txt files)
            #fun_write_dat(dat.files, catch.at.age, obs_w_err, params, sample.sizes, y)

            # This is a wrapper around the run.basa function that reruns the assessment procedured
            # if NUTS mysteriously fails or times out. This merely helps with making sure that the
            # simulations run to completion without random NUTS errors interrupting.

            desired.samples <- 1000   # this should be the number of iterations per chains * the number of chains (4)
            max.duration <- 20        # max runtime duration for the NUTS sample in minutes
            iters <- 1

            repeat{
                if(iters > 1){
                  Sys.sleep(5)
                }

                # print(paste("Trying with max duration of", max.duration, "minutes."))
                convergence.diags <- run.basa.adnuts(model.dir, sim.seed, n.iter=1500, n.warmup=500, max.duration = max.duration, n.chains=1)   # This is the important calculation

                # Check to make sure we actually generated the expected the number of samples
                # within the timeframe
                n.samples <- tryCatch({
                    nrow(read.csv("mcmc_out/PFRBiomass.csv", header=FALSE))
                }, error=function(e){
                    print("Insufficient samples recorded, retrying.")
                    0
                })

                # If the model converged and we got the right number of samples, dont repeat
                iters <- iters+1
                if((!is.null(convergence.diags) && 
                    convergence.diags$divergences < 0.01 && 
                    convergence.diags$converged && 
                    n.samples == desired.samples) || 
                    iters > 3){
                  break;
                }
                
                print(convergence.diags)
                print(paste0(cr.name, " (sim ", sim.seed, "): ", y, "/", stop.year, " failed, due to lack of convergence."))

                if(n.samples < desired.samples){
                    max.duration <- min(round(max.duration/(n.samples/desired.samples), 0)+1, 30)
                }

            }

            if(iters > 3){
              return(list(success=FALSE, message="Maximum Iterations Reached", last.year=y))
            }
        }

        # Read in the new data files and update the female.spawners param
        dat.files <- read.data.files(model.dir)
        params$female.spawners  <- tail(dat.files$PWS_ASA.dat$perc.female, 1)

        if(!exists("convergence.diags")){
            convergence.diags <- list(total.time=0)
        } 

        print(paste0(cr.name, " (sim ", sim.seed, "): ", y, "/", stop.year, " complete in ", round(as.numeric(convergence.diags$total.time), 2), " seconds."))

    }
    return(list(pop.dyn=pop_dyn, harvest.rate=control.rule, obs=obs_w_err, success=TRUE, last.year=y))
}

# source(file=paste0(here::here("R/operating_model/"),  "fun_operm.R"))
# source(file=paste0(here::here("R/operating_model/"),  "fun_obsm.R"))
  # cr <- list(type="hcr.threshold.linear",      lower.threshold=0, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.20)
  # sim.dir <- paste0(here::here("results"), "/test/")
  # sim.out.f.00 <- run.simulation(cr, nyr.sim=5, sim.seed=1998, write=sim.dir, assessment = FALSE, hindcast=FALSE)

# biomass <- apply(sim.out.f.00$pop.dyn$prefish.spawn.biomass, 1, sum)

# plot(1:25, biomass, type="l")
