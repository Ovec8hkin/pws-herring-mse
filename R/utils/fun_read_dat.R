# fun_read_dat.r
# Created by John Trochta
# Date created:  07/27/2020
# Summary:
# Project population dynamics to next year with given catches

library(tidyverse)

# Read BASA data files (.dat and .ctl) into list format for easy parameter access.
read.data.files <- function(dat.dir){

    filename <- vector(length=5)
    filename[1] <- paste0(dat.dir, "PWS_ASA.dat")
    filename[2] <- paste0(dat.dir, "PWS_ASA(ESS).ctl")
    filename[3] <- paste0(dat.dir, "PWS_ASA(covariate).ctl")
    filename[4] <- paste0(dat.dir, "agecomp_samp_sizes.txt")
    filename[5] <- paste0(dat.dir, "PWS_ASA_disease.dat")

    PWS_ASA.dat            <- data_reader(filename=filename[1])
    PWS_ASA_ESS.ctl        <- data_reader(filename=filename[2])
    PWS_ASA_covariate.ctl  <- data_reader(filename=filename[3])
    agecomp_samp_sizes.txt <- data_reader(filename=filename[4])
    PWS_ASA_disease.dat    <- data_reader(filename=filename[5])

    par.names <- c("nyr", "nyr_tobefit", "nage", "waa", "fecundity", 
                "pound_catch", "pk", "foodbait_catch", "gillnet_catch", "seine_yield",  "perc.female", 
                "mdm", "egg", "egg_se", "adfg_hydro_year_start", "adfg_hydro", "pwssc_hydro_year_start", "pwssc_hydro", "pwssc_hydro_se",
                "seine_age_comp", "spawn_age_comp", "juvenile_survey")
    names(PWS_ASA.dat) <- par.names

    par.names <- c("ctl_file", "seine_ess", "spawn_ess", "vhsv_ess", "ich_ess")
    names(PWS_ASA_ESS.ctl) <- par.names

    par.names <- c("std_covs", "num_recruit_cov", "recruit_fixed_rand", "cov_on", "regime_shift_89", "r_beta_change", 
                    "num_mort_covs", "mort_fixed_rand", "mort_on", "mort_age_impact", "disease_covs", "winter_mort_devs", "m_btea_change")
    names(PWS_ASA_covariate.ctl) <- par.names

    par.names <- c("seine_sample_size", "spawn_sample_size", "vhsv_sample_size", "ich_smaple_size")
    names(agecomp_samp_sizes.txt) <- par.names

    par.names <- c("vhsv_age_prevalence", "vhsv_obs_start", "vhsv_est_start", "vhsv_recov_prob", 
                   "ich_age_prevalence", "ich_obs_start", "ich_est_start", "ich_recov_prob")
    names(PWS_ASA_disease.dat) <- par.names

    return(listN(PWS_ASA.dat,
                  PWS_ASA_ESS.ctl,
                  PWS_ASA_covariate.ctl,
                  agecomp_samp_sizes.txt,
                  PWS_ASA_disease.dat))
}

# Read and coerce BASA data into appropriate R data type (vector, list, matrix, etc.).
data_reader <- function(filename) {
    #  The user needs to make sure there is a blank line at the end of
    #  each file and that each data type (vector number or matrix) is
    #  separated by a blank line
  
    # This is kind of convoluted
    text <- readLines(filename)
    values <- grep("^\\s{0,2}[0-9]",text)
    signed.values <- grep("^\\s{0,2}[-]",text)
    read.these <- sort(c(values,signed.values))
    nlines <- length(text) 
    indices <- seq(1:nlines)
    indices <- indices[read.these]
    first.differences <- c(diff(indices),5) # This accounts for the last data
    data.types <- length(first.differences[first.differences>1])

    data <- vector("list", data.types)
    j <- 1
    temp <- NA
    for(i in 1:length(indices)){
      temp.1 <- scan(filename,skip=indices[i]-1,nlines=1,quiet=TRUE,flush=FALSE)
      if(first.differences[i]>1){
        if(all(is.na(temp))){
          data[[j]] <- temp.1
          j <- j+1
        } else{
          temp <- rbind(temp,temp.1)
          data[[j]] <- temp
          rownames(data[[j]]) <- NULL
          temp <- NA
          j <- j+1
        }
      } else{
        if(all(is.na(temp))){
          temp <- temp.1
        } else{
          temp <- rbind(temp,temp.1)
        }
      }
    }
    return(data)
}

# Read BASA parameter estimates (from .par file) into name list
read.par.file <- function(filename){
    library(stringr)
    par.vals <- data_reader(filename)
    text <- readLines(filename)
    par.names <- text[grep("^#", text)]
    par.names <- apply(as.array(par.names[2:length(par.names)]), 1, function(x) stringr::str_match(x, "[a-zA-Z0-9_]+"))
    names(par.vals) <- par.names
    return(par.vals)
}

read.biomass.estimates <- function(model.dir, nyr=NA){
  fname <- paste0(here::here(model.dir), "/mcmc_out/PFRBiomass.csv")
  biomass.data <- read.table(fname, header = FALSE, sep = ",", dec=".")

  nyr <- ifelse(is.na(nyr), ncol(biomass.data)-1, nyr)
  years <- 1980:(1980+nyr)

  colnames(biomass.data) <- years

  return(biomass.data[,1:nyr])

}

read.exploit.rates <- function(model.dir, nyr=NA){
  raw.data <- read.data.files(model.dir)

  pfrb <- read.biomass.estimates(model.dir, nyr=nyr)
  nyr <- ifelse(is.na(nyr), ncol(pfrb), nyr)

  total.catch.biomass <- compute.catch.biomass(raw.data$PWS_ASA.dat, nyr)
  exploit.rate <- t(total.catch.biomass/t(pfrb))

  return(exploit.rate)

}

read.true.exploit.rates <- function(cr, sims, nyr){
    catch <- read.catch.data(cr, sims, nyr) %>% na.omit() %>%
        group_by(year, control.rule, sim) %>%
        summarise(catch=sum(catch)) %>%
        mutate(year=2022+year)
    biomass <- read.true.biomass.data(cr, sims, nyr)

    exploit.rates <- catch %>% 
        left_join(biomass, by=c("year", "control.rule", "sim")) %>%
        mutate(
          exploit = catch/biomass
        ) %>%
        print(n=10)

    return(exploit.rates %>% na.omit())
}

compute.catch.biomass <- function(data, nyr){

    weight.at.age <- data$waa[1:nyr,]

    fb.nya <- data$foodbait_catch[1:nyr,]
    pound.nya <- data$pound_catch[1:nyr,]
    gillnet.nya <- data$gillnet_catch[1:nyr,]
    seine.yield  <- data$seine_yield[1:nyr]
    
    fb.nya.biomass <- weight.at.age * fb.nya 
    pound.nya.biomass <- weight.at.age * pound.nya
    gillnet.nya.biomass <- weight.at.age * gillnet.nya 

    fb.biomass.annual       <- rowSums(fb.nya.biomass) # Now sum the biomass over all age classes for each year
    pound.biomass.annual    <- rowSums(pound.nya.biomass)
    gillnet.biomass.annual  <- rowSums(gillnet.nya.biomass)

    fb.biomass.annual       <- replace(fb.biomass.annual, fb.biomass.annual == 0, NA)
    pound.biomass.annual    <- replace(pound.biomass.annual, pound.biomass.annual == 0, NA)
    gillnet.biomass.annual  <- replace(gillnet.biomass.annual, gillnet.biomass.annual == 0, NA)
    seine.yield             <- replace(seine.yield, seine.yield == 0, NA)
    
    # Matrix of catches by gear type in mt
    total.catch <- cbind(fb.biomass.annual, pound.biomass.annual, gillnet.biomass.annual, seine.yield)
    total.catch[is.na(total.catch)] <- 0 
    total.catch.biomass <- rowSums(total.catch) # total catches by year in mt

    return(total.catch.biomass)

}

read.survey.estimates <- function(model.dir){

    nage=10

    fnames <- list(
        mdm             = paste0(model.dir, "/MDM.csv"),
        pwssc.hydro     = paste0(model.dir, "/HYD_PWSSC.csv"),
        adfg.hydro      = paste0(model.dir, "/HYD_ADFG.csv"),
        spawn.age.comp  = paste0(model.dir, "/SpAC.csv"),
        seine.age.comp  = paste0(model.dir, "/SeAC.csv"),
        egg             = paste0(model.dir, "/EGG.csv"),
        juv.schools     = paste0(model.dir, "/juv_schools.csv")
    )

    spac.data <- read.csv(fnames$spawn.age.comp, header=FALSE)
    spac.data <- apply(spac.data, 2, median)
    spac.data <- matrix(spac.data, ncol=nage, byrow=TRUE)

    seac.data <- read.csv(fnames$seine.age.comp, header=FALSE)
    seac.data <- apply(seac.data, 2, median)
    seac.data <- matrix(seac.data, ncol=nage, byrow=TRUE)

    mdm.data <- read.csv(fnames$mdm, header=FALSE)
    mdm.data <- apply(mdm.data, 2, median)
    mdm.data <- as.vector(mdm.data)

    pwssc.hydro.data <- read.csv(fnames$pwssc.hydro, header=FALSE)
    pwssc.hydro.data <- apply(pwssc.hydro.data, 2, median)
    pwssc.hydro.data <- as.vector(pwssc.hydro.data)

    adfg.hydro.data <- read.csv(fnames$adfg.hydro, header=FALSE)
    adfg.hydro.data <- apply(adfg.hydro.data, 2, median)
    adfg.hydro.data <- as.vector(adfg.hydro.data)

    egg.data <- read.csv(fnames$egg, header=FALSE)
    egg.data <- apply(egg.data, 2, median)
    egg.data <- as.vector(egg.data)
    egg.data[egg.data == 0] <- -9

    juv.schools.data <- read.csv(fnames$juv.schools, header=FALSE)
    juv.schools.data <- apply(juv.schools.data, 2, median)
    juv.schools.data <- as.vector(juv.schools.data)

    return(listN(mdm.data, egg.data, pwssc.hydro.data, adfg.hydro.data, juv.schools.data, spac.data, seac.data))
}

read.catch.data <- function(cr, sims, nyr){

    fnames <- c("foodbait_catch.csv", "gillnet_catch.csv", "pound_catch.csv", "seine_catch.csv")

    data <- data.frame(year=NA, catch=NA, fishery=NA, control.rule=NA, sim=NA)
    for(s in sims){
        for(f in fnames){
            fname <- paste0(here::here("results/"), cr, "/sim_", s, "/year_", nyr, "/results/", f)
            dat.fname <- paste0(here::here("results/"), cr, "/sim_", s, "/year_", nyr, "/model/")

            if(!file.exists(dat.fname)) next;
            
            waa <- read.data.files(dat.fname)$PWS_ASA.dat[[4]]
            waa <- waa[(42+1):(42+1+nyr-1),]
            catch.data <- as.matrix(read.csv(fname))[1:nyr, 2:11]
            catch.biomass <- apply(catch.data*waa, 1, sum)

            d <- data.frame(year=1:nyr, catch=catch.biomass, fishery=str_split(f, "_")[[1]][1], control.rule=cr, sim=s)
            data <- data %>% bind_rows(d)
        }
    }
    return(data)

}

read.biomass.data <- function(cr, sims, nyr){
    data <- data.frame(year=NA, biomass=NA, control.rule=NA, sim=NA) 
    for(s in sims){
        for(i in 1:nyr){
            fname <- paste0(here::here("results/"), cr, "/sim_", s, "/year_", i, "/model/mcmc_out/PFRBiomass.csv")
            
            if(!file.exists(fname)) next;

            biomass.data <- read_csv(fname, col_names=FALSE, show_col_types = FALSE, col_select = last_col()) %>%
                            rename(biomass=1) %>%
                            mutate(
                              year=as.character(2022+i),
                              control.rule=cr,
                              sim=s
                            ) %>%
                            relocate(c(year, biomass, control.rule, sim))

            data <- data %>% bind_rows(biomass.data)
        }
    }
    return(data %>% na.omit())
}

read.true.biomass.data <- function(cr, sims, nyr){
    data <- data.frame(year=NA, biomass=NA, control.rule=NA, sim=NA)
    for(s in sims){

        sim.name <- paste0(here::here("results/"), cr, "/sim_", s)
        years <- list.files(sim.name)
        ys <- sort(apply(as.matrix(years), 1, get.year))
        final.year <- ys[length(ys)]

        final.year <- min(final.year, nyr)

        fname <- paste0(sim.name, "/year_", final.year, "/results/prefish_spawn_biomass.csv")
        if(!file.exists(fname)) next;

        biomass.data <- read_csv(fname, col_names=TRUE, show_col_types = FALSE) %>% 
                          slice(1:nyr) %>%
                          select(-c(1)) %>%
                          mutate(
                            year=2022+row_number(),
                            control.rule=as.character(cr),
                            sim=s
                          ) %>%
                          rowwise() %>%
                          mutate(
                            biomass = sum(c(`0`, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`))
                          ) %>%
                          select(year, biomass, control.rule, sim)
        
        data <- data %>% bind_rows(biomass.data)
    }
    return(data %>% na.omit())
}

get.year <- function(f){
    return(as.numeric(str_split(f, "_")[[1]][2]))
}


# Function for simultaneously creating names of variables/elements within a list
listN <- function(...){
  anonList <- list(...)
  names(anonList) <- as.character(substitute(list(...)))[-1]
  anonList
}