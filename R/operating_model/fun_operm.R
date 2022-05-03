# fun_operm.r
# Created by John Trochta
# Date created:  07/27/2020
# Summary:
# Project population dynamics to next year with given catches

nage <- 10

fun_operm <- function(year, pop.dyn, true.catch.at.age, pars, sim.seed=NA){
  
    
    curr.nya <- pop.dyn$true.nya[year, ]

    # Calculate summer and winter survival
    survival <- calc.survival(pop.dyn$survival.summer[year-1, ], pop.dyn$survival.winter[year-1, ], pars$Z_0_8, pars$Z_9, sim.seed)

    # Calculate maturity
    maturity <- calc.maturity(pars$mat_par_1, pars$mat_par_2)
    
    # Calculate sac-roe seine fishery selectivity
    seine.vuln <- 1/(1+exp(-1.0*pars$seine_vul_beta*(0:(nage-1)-pars$seine_vul_alpha)))
    
    # Calculate pre-fishery biomass
    prefish.spawn.biomass <- maturity * pars$waa * curr.nya
    
    # Re-assign true.catch.at.age to match variables as in BASA (EM)
    catches <- calc.catches(true.catch.at.age, pars$pk)
    
    # Calculate post-fishery spawners (for milt survey)
    n.spawners <- maturity*(curr.nya-catches$spring.total)
    spawn.biomass.age.comp <- pars$waa*n.spawners
    spawn.biomass <- sum(spawn.biomass.age.comp)
    
    # Calculate Nya for next year
    new.nya <- calc.nya(pars$log_MeanAge0, pop.dyn$annual.age0.devs[year], pars$sigma_age0devs, pars$pk, curr.nya, catches, survival)
    
    # Calculate Mile-days of milt
    mdm <- (1-pars$female.spawners)*spawn.biomass/exp(pars$logmdm_c)
    
    # Calculate PWSSC hydroacoustic survey
    PWSSC.hydro <- exp(pars$PWSSC_hydro_q)*sum(prefish.spawn.biomass)
    
    # Spawner survey age-comps
    spawn.age.comp <- (maturity * curr.nya)/sum(maturity * curr.nya)
    
    # Seine fishery age-comps
    if(all(catches$seine)==0){
        seine.age.comp <- -9
    }else{
        seine.age.comp <- (catches$seine)/sum(catches$seine)
    }
    
    # Calculate seine fishery biomass
    seine.biomass <- sum(catches$seine.catch*pars$waa)
    
    # Discontinued surveys
    egg <- -9 # egg deposition                  10^-6*pars$female.spawners*sum(curr.nya*pop.dyn$fecundity[y,])
    ADFG.hydro <- -9 # ADF&G hydroacoustic      exp(pars$ADFG_hydro_q)*sum(prefish.spawn.biomass)
    
    ## Reassign all local matrices and vectors to their pop.dyn equivalents here
    pop.dyn$survival.summer[year, ]         <- survival$summer
    pop.dyn$survival.winter[year, ]         <- survival$winter
    pop.dyn$maturity[year, ]                <- maturity
    pop.dyn$prefish.spawn.biomass[year+1, ] <- prefish.spawn.biomass
    pop.dyn$seine.catch[year, ]             <- catches$seine
    pop.dyn$gillnet.catch[year, ]           <- catches$gillnet
    pop.dyn$pound.catch[year, ]             <- catches$pound
    pop.dyn$foodbait.catch[year, ]          <- catches$foodbait
    pop.dyn$n.spawners[year, ]              <- n.spawners
    pop.dyn$spawn.biomass.age.comp[year, ]  <- spawn.biomass.age.comp
    pop.dyn$true.nya[year+1, ]              <- new.nya
    pop.dyn$spawn.age.comp[year, ]          <- spawn.age.comp
    pop.dyn$seine.age.comp[year, ]          <- seine.age.comp

    pop.dyn$spawn.biomass[year]           <- spawn.biomass
    pop.dyn$seine.biomass[year]           <- seine.biomass
    pop.dyn$mdm[year]                     <- mdm
    pop.dyn$egg[year]                     <- egg
    pop.dyn$PWSSC.hydro[year]             <- PWSSC.hydro
    pop.dyn$ADFG.hydro[year]              <- ADFG.hydro

    return(pop.dyn)
}

## Needs to be revised in the future
## 1. Shoudn't include random error not estaimted by the model
## 2. Should include disease information (I think)
calc.survival <- function(prev.survival.summer, prev.survival.winter, Z_0_8, Z_9, sim.seed){

    if(!is.na(sim.seed)){
        set.seed(sim.seed)
    }
    error <- c(rep(rnorm(1, mean=0, sd=0.05), 3), rep(rnorm(1, mean=0, sd=0.05), 6), rnorm(1, 0, 0.05))

    survival.summer <- rep(0, nage)
    survival.winter <- rep(0, nage)

    survival.summer[1:(nage-1)] <- rep(exp(-0.5*Z_0_8), nage-1)+error[1:(nage-1)]
    survival.winter[1:(nage-1)] <- rep(exp(-0.5*Z_0_8), nage-1)+error[1:(nage-1)]

    survival.summer[nage] <- ifelse(length(prev.survival.summer) == 0, 
                                    exp(-0.5*Z_9),
                                    prev.survival.summer[nage]*(survival.summer[nage-1]/prev.survival.summer[nage-1])
    )
                                    #survival.summer[nage-1]*prev.survival.summer[nage]/prev.survival.summer[nage-1])
    survival.winter[nage] <- ifelse(length(prev.survival.winter) == 0, 
                                    exp(-0.5*Z_9),
                                    prev.survival.winter[nage]*(survival.winter[nage-1]/prev.survival.winter[nage-1])
    )

    max.survival.threshhold <- 0.95
    survival.summer[survival.summer > max.survival.threshhold] <- max.survival.threshhold
    survival.winter[survival.winter > max.survival.threshhold] <- max.survival.threshhold

    return(list(summer=survival.summer, winter=survival.winter))
}

calc.maturity <- function(matur_age3_per1, matur_age4_per1){
    
    maturity <- rep(0, nage)
    maturity[1:3]    <- 0
    maturity[4]      <- matur_age3_per1*matur_age4_per1
    maturity[5]      <- matur_age4_per1
    maturity[6:nage] <- 1

    return(maturity)
}

calc.catches <- function(true.caa, pk){

    seine.catch     <- rep(0, nage)
    gillnet.catch   <- rep(0, nage)
    pound.catch     <- rep(0, nage)
    foodbait.catch  <- rep(0, nage)

    seine.catch    <- true.caa[1, ]
    gillnet.catch  <- true.caa[2, ]
    pound.catch    <- true.caa[3, ]
    foodbait.catch <- true.caa[4, ]

    spring.caa <- seine.catch+gillnet.catch+pk*pound.catch

    return(
        list(seine=seine.catch, gillnet=gillnet.catch, pound=pound.catch, foodbait=foodbait.catch, 
             spring.total=spring.caa
        )
    )

}

## Needs to reconsider age-3 recruitment
calc.nya <- function(log.mean.age0, annual.age0.devs, sigma.age0.devs, pk, curr.nya, catches, survival){
    new.nya <- rep(0, nage)

    # Add new age-3 recruits
    new.nya[1] <- calc.recruitment(log.mean.age0, annual.age0.devs, sigma.age0.devs)
    
    # Project new naa (numbers at age) for next year
    #spring.catch <- pop.dyn$seine.catch[year, 1:(nage-2)]+pop.dyn$gillnet.catch[year, 1:(nage-2)]+pk*pop.dyn$pound.catch[year, 1:(nage-2)]

    new.nya[2:(nage-1)] <- (((curr.nya[1:(nage-2)]-(catches$spring.total[1:(nage-2)]))*survival$summer[1:(nage-2)])-catches$foodbait[1:(nage-2)])*survival$winter[1:(nage-2)]
    
    # Calculate plus-group numbers
    spring.catch.8 <- catches$seine[nage-1]+catches$gillnet[nage-1]+pk*catches$pound[nage-1]
    spring.catch.9 <- catches$seine[nage]+catches$gillnet[nage]+pk*catches$pound[nage]
    new.nya[nage] <- ((curr.nya[nage-1]-(spring.catch.8))*survival$summer[nage-1]-catches$foodbait[nage-1])*survival$winter[nage-1] +
                     ((curr.nya[nage]-(spring.catch.9))*survival$summer[nage]-catches$foodbait[nage])*survival$winter[nage]

    new.nya[new.nya < 0] <- 0

    return(new.nya)
}

## Needs to be revised to provide additional variability
## 1. Change recruitment variability (sigmag.age0.devs)?
## 2. Include stock-recruitment curve (beverton-holt w/ steepness)?
## 3. Add temporal autocorrelation?
## 4. Allow for larger recrtuiment events more often?
## 5. Include recruitment regine shift?
calc.recruitment <- function(log.mean.age0, annual.age0.devs, sigma.age0.devs){
    return(exp(log.mean.age0 + annual.age0.devs - 0.5*sigma.age0.devs^2))
}