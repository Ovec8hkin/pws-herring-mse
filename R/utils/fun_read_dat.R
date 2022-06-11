# fun_read_dat.r
# Created by John Trochta
# Date created:  07/27/2020
# Summary:
# Project population dynamics to next year with given catches

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

accumulate.results.data <- function(nyr.sim, total.sims, seeds, trial, fnames, byage=FALSE, ncols=1){

    if(length(byage) == 1){
      byage <- rep(byage, length(fnames))
    }

    #d.vec <- ifelse(byage, 10, 1)

    data.matrices <- vector("list", length(fnames))
    for(i in 1:length(fnames)){
        d <- ifelse(byage[i], 10, ncols)
        data.matrices[[i]] <- array(NA, dim=c(nyr.sim, d, total.sims))
    }

    for(i in 1:length(fnames)){
        for(s in 1:length(seeds)){
            f <- paste0(here::here("results"), "/", trial, "/sim_", seeds[s], "/year_", nyr.sim, "/results/", fnames[i])
            data.matrices[[i]][,,s] <- as.matrix(read.csv(f)[1:nyr.sim,-1]) 
        }
        
    }

    return(data.matrices)

}

accumulate.assessment.posteriors <- function(nyr.sim, total.sims, seeds, trial, fname="PFRBiomass.csv"){

    f1 <- paste0(here::here("results"), "/", trial, "/sim_", seeds[1], "/year_1/model/mcmc_out/", fname)
    n.samps <- nrow(read.csv(f1, header=FALSE))
    #print(n.samps)
    data.matrices <- vector("list", length(fname))
    for(i in 1:length(fname)){
        data.matrices[[i]] <- array(NA, dim=c(n.samps, nyr.sim, total.sims))
    }

    for(i in 1:nyr.sim){
        for(s in 1:length(seeds)){
            print(i)
            f <- paste0(here::here("results"), "/", trial, "/sim_", seeds[s], "/year_", i-1, "/model/mcmc_out/", fname)
            #print(f)
            data <- read.csv(f, header=FALSE)
            data.matrices[[1]][,i,s] <- as.matrix(data[(nrow(data)-n.samps+1):nrow(data), ncol(data)])
        } 
    }
    
    return(data.matrices)
}


# Function for simultaneously creating names of variables/elements within a list
listN <- function(...){
  anonList <- list(...)
  names(anonList) <- as.character(substitute(list(...)))[-1]
  anonList
}