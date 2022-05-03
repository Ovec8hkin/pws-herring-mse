# fun_write_dat.r
# Created by John Trochta
# Date created:  07/27/2020
# Summary:
# Project population dynamics to next year with given catches
source(file=paste0(here::here("R/utils"), "/fun_read_dat.R"))

fun_write_dat <- function(dat.files, catches, obs, params, ss, year){
  
  options(scipen=10)

  if(-9 %in% obs$seine.age.comp){ss$seine <- 0}
  if(-9 %in% obs$spawn.age.comp){ss$spawn <- 0}

  # Write main data file with actual catches (known w/ error)
  # and observed survey values (from fun_obsm.R).
  PWS.ASA.dat <- dat.files$PWS_ASA.dat
  # Add new data 
  PWS.ASA.dat[[1]]  <- PWS.ASA.dat[[1]]+1                                         # Number of year in data
  PWS.ASA.dat[[2]]  <- PWS.ASA.dat[[2]]+1                                         # Number of years to be fit to
  # 3 is number of age classes (10). Needs no change.
  PWS.ASA.dat[[4]]  <- rbind(PWS.ASA.dat[[4]],  obs$waa)                          # Weight-at-age (assumed constant 03/2022) 
  PWS.ASA.dat[[5]]  <- rbind(PWS.ASA.dat[[5]],  rep(-9, 10))                      # Fecundity (no data post 1997) 
  PWS.ASA.dat[[6]]  <- rbind(PWS.ASA.dat[[6]],  catches[3, ])                     # Pound catch fishery
  PWS.ASA.dat[[7]]  <- params$pk                                                  # Prop of pound-catch killed (Assumed constant)
  PWS.ASA.dat[[8]]  <- rbind(PWS.ASA.dat[[8]],  catches[4, ])                     # Foodbait fishery
  PWS.ASA.dat[[9]]  <- rbind(PWS.ASA.dat[[9]],  catches[2, ])                     # Gillnet fishery
  PWS.ASA.dat[[10]] <- rbind(PWS.ASA.dat[[10]], sum(catches[1, ]*obs$waa))        # Seine net fishery
  PWS.ASA.dat[[11]] <- rbind(PWS.ASA.dat[[11]], obs$female)                       # % female

  PWS.ASA.dat[[12]] <- rbind(PWS.ASA.dat[[12]], obs$mdm)                          # Observed MDM
  PWS.ASA.dat[[13]] <- rbind(PWS.ASA.dat[[13]], obs$egg)                          # Observed egg deposition
  PWS.ASA.dat[[14]] <- rbind(PWS.ASA.dat[[14]], obs$egg.sd)                       # Observed egg deposition sd
  # 15 is first year of ADFG hydro survey. Needs no change.
  PWS.ASA.dat[[16]] <- rbind(PWS.ASA.dat[[16]], obs$adfg.hydro)                   # Observed ADFG hydroacoustic survey biomass
  # 17 is first year of PWSSC hydro survey. Needs no change.
  PWS.ASA.dat[[18]] <- rbind(PWS.ASA.dat[[18]], obs$pwssc.hydro)                  # Observed PWSSC hydroacoustic survey biomass
  PWS.ASA.dat[[19]] <- rbind(PWS.ASA.dat[[19]], obs$pwssc.hydro.sd)               # Observed PWSSC hydroacoustic survey biomass sd
  PWS.ASA.dat[[20]] <- rbind(PWS.ASA.dat[[20]], obs$seine.age.comp)               # Observed seine age composition
  PWS.ASA.dat[[21]] <- rbind(PWS.ASA.dat[[21]], obs$spawn.age.comp)               # Observed spawner age composition
  PWS.ASA.dat[[22]] <- rbind(PWS.ASA.dat[[22]], obs$aerial.juv)                   # Observed juvenile aerial survey
  

  write.data(PWS.ASA.dat, "PWS_ASA.dat")

  # Write ESS data file with actual sample sizes (fixed).
  PWS.ASA.ESS.ctl <- dat.files$PWS_ASA_ESS.ctl 
  PWS.ASA.ESS.ctl[[2]] <- rbind(PWS.ASA.ESS.ctl[[2]], -9)
  PWS.ASA.ESS.ctl[[3]] <- rbind(PWS.ASA.ESS.ctl[[3]], 30)
  PWS.ASA.ESS.ctl[[4]] <- rbind(PWS.ASA.ESS.ctl[[4]], 0)
  PWS.ASA.ESS.ctl[[5]] <- rbind(PWS.ASA.ESS.ctl[[5]], 0)
  # Write file without headers
  write.data(PWS.ASA.ESS.ctl, "PWS_ASA(ESS).ctl")
  
  
  # Write covariate file (not currently used)
  PWS.ASA.covariate.ctl <- dat.files$PWS_ASA_covariate.ctl
  # Skip indices 1-4 for now, as they do not change with year
  PWS.ASA.covariate.ctl[[5]] <- rbind(PWS.ASA.covariate.ctl[[5]], -1) # Covariate for a regime shift in recruitment
  PWS.ASA.covariate.ctl[[6]] <- rbind(PWS.ASA.covariate.ctl[[6]], 0)
  # Skip indices 7-11 for now, as they do not change with year
  PWS.ASA.covariate.ctl[[12]] <- rbind(PWS.ASA.covariate.ctl[[12]], c(-9,0,0)) # These are the disease prevalence indices
  PWS.ASA.covariate.ctl[[14]] <- rbind(PWS.ASA.covariate.ctl[[14]], 0)
  
  # Write .dat file without headers
  write.data(PWS.ASA.covariate.ctl, "PWS_ASA(covariate).ctl")
  
  # Write age comp sample sizes file with actual 
  # sample sizes (currently fixed).
  age.comp.samp.sizes.txt <- dat.files$agecomp_samp_sizes.txt
  age.comp.samp.sizes.txt[[1]] <- rbind(age.comp.samp.sizes.txt[[1]], ss$seine)
  age.comp.samp.sizes.txt[[2]] <- rbind(age.comp.samp.sizes.txt[[2]], ss$spawn)
  age.comp.samp.sizes.txt[[3]] <- rbind(age.comp.samp.sizes.txt[[3]], 0)
  age.comp.samp.sizes.txt[[4]] <- rbind(age.comp.samp.sizes.txt[[4]], 0)
  # Write file without headers
  write.data(age.comp.samp.sizes.txt, "agecomp_samp_sizes.txt")
  
  PWS.ASA.disease.dat <- dat.files$PWS_ASA_disease.dat
  PWS.ASA.disease.dat[[1]] <- rbind(PWS.ASA.disease.dat[[1]], obs$vhsv.antibody)  # Observed VHSV antibodies
  # Skip indices 2-4, as they do not change with year.
  PWS.ASA.disease.dat[[5]] <- rbind(PWS.ASA.disease.dat[[5]], obs$ich.antibody)   # Observed ichthyophonus antibodies
  # Skip indices 6-8, as they do not change with year.
  write.data(PWS.ASA.disease.dat, "PWS_ASA_disease.dat")

  # Function for simultaneously creating names of variables/elements within a list
  listN <- function(...){
    anonList <- list(...)
    names(anonList) <- as.character(substitute(list(...)))[-1]
    anonList
  }

  options(scipen=0)
  
  return(listN(PWS.ASA.dat,
               PWS.ASA.ESS.ctl,
               PWS.ASA.covariate.ctl,
               age.comp.samp.sizes.txt))
}

write.data <- function(data, fname){
    write.table(
      rbind(" ", data[[1]], " "),
      file = fname, 
      append = F, sep = " ",
      row.names=FALSE, col.names=FALSE,
      quote=F
    )
    for(i in 2:length(data)){
        write.table(
          rbind(data[[i]], " "),
          file = fname, 
          append = T, sep = " ",
          row.names=FALSE, col.names=FALSE,
          quote=F
        )
    }
}