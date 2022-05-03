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

files.sources = list.files(here::here("R/operating_model/control_rules"), full.names=TRUE)
sapply(files.sources, source)

#source(file=paste0(here::here("R/utils/"),            "fun_data_reader.R"))

nyr.sim <- 5
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
  prefish.spawn.biomass   <- matrix(0, nyr.sim+1, nage, dimnames=list(c(rownames, nyr.sim+1), colnames))
  seine.catch             <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  gillnet.catch           <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  pound.catch             <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  foodbait.catch          <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  n.spawners              <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  spawn.biomass.age.comp  <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  true.nya                <- matrix(0, nyr.sim+1, nage, dimnames=list(c(rownames, nyr.sim+1), colnames))
  spawn.age.comp          <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  seine.age.comp          <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))

  spawn.biomass           <- rep(0, nyr.sim)
  mdm                     <- rep(0, nyr.sim)
  PWSSC.hydro             <- rep(0, nyr.sim)
  seine.biomass           <- rep(0, nyr.sim)
  egg                     <- rep(0, nyr.sim)
  ADFG.hydro              <- rep(0, nyr.sim)
  annual.age0.devs        <- rep(0, nyr.sim)

  pop_dyn <- listN(survival.summer, survival.winter, maturity, prefish.spawn.biomass, 
                seine.catch, gillnet.catch, pound.catch, foodbait.catch,
                n.spawners, spawn.biomass.age.comp, spawn.biomass, true.nya, 
                mdm, PWSSC.hydro, spawn.age.comp, seine.age.comp, seine.biomass, egg, ADFG.hydro,
                annual.age0.devs)

  return(pop_dyn)
  
}

#set.initial.conditions <- function(pop_dyn, ssb, waa, ssb.nya.conversion=0.03, sim.seed=NA){
set.initial.conditions <- function(dir, pop_dyn){ 
  # true.pop.size <- ssb*ssb.nya.conversion
  # nya.probs <- c(0.206, 0.177, 0.172, 0.136, 0.104, 0.0779, 0.055, 0.030, 0.017, 0.025)
  # if(!is.na(sim.seed)){
  #   set.seed(sim.seed)
  # }
  # true.nya <- as.vector(rmultinom(1, size=true.pop.size, prob=nya.probs))

  year.0.est <- get.assessment.estimates(dir, 0)

  prop.age.structure <- as.numeric(year.0.est$est.nya[3, ]/sum(year.0.est$est.nya[3, ]))

  pop_dyn$true.nya[1, ] <- as.numeric(year.0.est$est.nya[3, ]) # Should these be rounded to integers?
  pop_dyn$prefish.spawn.biomass[1, ] <- as.numeric(year.0.est$est.ssb[3, ])*prop.age.structure
  
  return(pop_dyn)
}

# Computed median weight-at-age from raw WAA data
# over the years 2012-2022. WAA is assumed stable
# on an annual basis.
calculate.waa <- function(dat.files){
  waa.data <- dat.files$PWS_ASA.dat[[4]]
  true.waa <- apply(waa.data[31:41, ], 2, median)
  return(true.waa)
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

run.simulation <- function(hcr.options, nyr.sim, sim.seed=NA, write=NA){
    print(write)
    if(is.na(write)){
      write <- paste0(here::here("results"), "/test")
    }

    # Copy over current stock assessment to use for starting values
    basa.root.dir <- "~/Desktop/Projects/basa/model/"
    model.0.dir <- create.model.dir(write, 0)
    file.remove(model.0.dir) # This is a hack to remove the internal model/ directory
    file.symlink(basa.root.dir, paste0(write, "/year_0/"))

    # Read in data files JUST ONCE, store then write within code
    dat.files <- read.data.files(model.0.dir)

    # These are the observed (not effective) sample sizes 
    seine.comp.obs.ss <- 500
    spawn.comp.obs.ss <- 1500

    # Read these out of the mcmc_out/*.csv files
    # Can be arbitary if non-operational. Can use median naa and waa from past n years
    # to set starting values
    true.waa <- calculate.waa(dat.files)
    pop_dyn <- initialize.popdyn.variables(nyr.sim)

    # Initialize projection estimates (e.g. from forecast model)
    pop_dyn <- set.initial.conditions(write, pop_dyn)

    # This can be back-calculated from CAA = exploitation * selectivity * naa
    # Probably will need to be fine-tuned later.
    # Can assume fully-selected fishery by setting selectivity=1
    # Would need to add some error to the CAA.
    fish.selectivity <- matrix(1, nrow=4, ncol=10)

    params <- read.par.file("~/Desktop/Projects/basa/model/PWS_ASA.par")
    params$female.spawners  <- tail(dat.files$PWS_ASA.dat$perc.female, 1)
    params$pk               <- dat.files$PWS_ASA.dat$pk
    params$waa              <- true.waa
    params$selectivity      <- fish.selectivity
    params$catch.sd         <- 0.1

    age0.error <- 0
    if(!is.na(sim.seed)){
      set.seed(sim.seed)
      age0.error <- rnorm(nyr.sim, 0, 0.5)
    }
    # Generate a new value for annual_age0dev. This will get replaced by a call to BASA
    # eventually, but this works for now.
    pop_dyn$annual.age0.devs <- median(params$annual_age0devs) + age0.error

    # Start loop
    control.rule <- rep(0, nyr.sim)
    ass.biomass <- data.frame(matrix(0, nrow=nyr.sim, ncol=6))
    colnames(ass.biomass) <- c("Year", "Biomass2.5", "Biomass25", "Biomass50", "Biomass75", "Biomass97.5")
    ass.biomass$Year <- 1:nyr.sim

    for(y in 1:nyr.sim){  

      model.dir <- create.model.dir(write, y)
      setwd(model.dir)

      # Perform management procedure
      ass.est <- get.assessment.estimates(write, y-1)
      ass.biomass[y, 2:6] <- t(ass.est$est.ssb)

      hcr.name <- hcr.options$type
      tar_hr <- match.fun(hcr.name)(as.numeric(ass.est$est.ssb[3, ]), as.numeric(ass.est$est.nya[3, ]), hcr.options)
      control.rule[y] <- tar_hr

      true.ssb <- sum(pop_dyn$prefish.spawn.biomass[y, ])
      true.nya <- pop_dyn$true.nya[y, ]

      # Execute fishery following management recommendation
      catch.at.age <- fun_fish(tar_hr, true.ssb, true.nya, params$waa, params$selectivity, params$catch.sd)
      
      # Run operating model
      pop_dyn <- fun_operm(y, pop_dyn, catch.at.age, params, sim.seed)

      # Generate observations with error
      obs_w_err <- fun_obsm(pop_dyn, dat.files$PWS_ASA.dat$waa, dat.files$PWS_ASA.dat$perc.female, y)

      # Write all new data to .dat file (and modify covariate and agecomp_samp_sizes.txt files)
      fun_write_dat(dat.files, catch.at.age, obs_w_err, params, list(seine=seine.comp.obs.ss, spawn=spawn.comp.obs.ss), y)

      # Run BASA
      #  - allow option to skip years of assessment
      convergence.diags <- run.basa.adnuts(model.dir, sim.seed)
      if(convergence.diags$divergences >= 0.005 || convergence.diags$converged == FALSE){
          print("Convergence failed")
          return("Convergence failed")
      }

      dat.files <- read.data.files(model.dir)
      params$female.spawners  <- tail(dat.files$PWS_ASA.dat$perc.female, 1)
      
    }

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
        write.csv(control.rule, "harvest.csv")
        write.csv(ass.biomass, "assessment_biomass.csv")

    }

    return(list(pop.dyn=pop_dyn, harvest.rate=control.rule, obs=obs_w_err))
  }

  cr <- list(type="hcr.hockey.stick", lower.threshold=10000, upper.threshold=30000, min.harvest = 0.0, max.harvest=0.20)
  sim.dir <- paste0(here::here("results"), "/", "lower.b0", "/sim_", "519", "/")
  sim.out.f.00 <- run.simulation(cr, nyr.sim=15, sim.seed=519, write=sim.dir)

