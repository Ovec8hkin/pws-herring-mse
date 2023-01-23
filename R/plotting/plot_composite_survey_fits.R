library(ggplot2)
library(ggdist)
library(here)
library(dplyr)
library(magrittr)
library(tidyverse)

source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))

# Calculates 95% confidence intervals from survey CV (from Buckland 1992 in References)
calc.buck.cv <- function(cv) {
    return(exp(1.96*sqrt(log(1+(cv^2)))))
}

#Computes posterior predictive intervals for log normal
post.pred.lognorm <- function(model.pred, sd){
    PP <- matrix(NA,nrow=nrow(model.pred),ncol=ncol(model.pred))
    for(i in 1:nrow(model.pred)){
        fits <- as.numeric(model.pred[i,]) # Take the true value of the survey estimate
        errors <- as.numeric(sd[i]*fits) # Take the error based on the true value of the survey estimate
        
        # Need to reparameterize to calculate the log-normal mean and sd for the survey data
        # From: https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
        E <- log(fits^2/sqrt(errors^2+fits^2))
        SD <- sqrt(log(1+(errors^2/fits^2)))
        PP[i,] <- rlnorm(n=length(fits),meanlog=E,sdlog=SD)
        # MDM.PP[i,] <- rnorm(n=length(fits),mean=E,sd=SD) # Generate a random variable based on the expected value (predictive) using the assumed 
    }
    return(PP)
}

post.pred.negbin <- function(model.pred, sd){
    PP <- matrix(NA,nrow=nrow(model.pred),ncol=ncol(model.pred))
    set.seed(100)
    for(i in 1:nrow(model.pred)){
        fits <- as.numeric(model.pred[i,]) # Take the true value of the survey estimate
        dispersion <- as.numeric(sd[i]) # Take the error based on the true value of the survey estimate
        
        E <- fits
        SD <- exp(dispersion)
        PP[i,] <- rnbinom(n=length(fits),size=SD,prob=SD/(SD+E))
        # MDM.PP[i,] <- rnorm(n=length(fits),mean=E,sd=SD) # Generate a random variable based on the expected value (predictive) using the assumed 
    }
    return(PP)
}

generate.post.pred <- function(fname, var, years, nburn=1, dist="lognorm", mask=NA){
    
    if(is.vector(fname)){
        data <- read.table(fname[1], header = FALSE, sep = ",", dec=".")[-c(1:nburn), 1:length(years)]
        for(i in 2:length(fname)){
            d <- read.table(fname[i], header = FALSE, sep = ",", dec=".")[-c(1:nburn), 1:length(years)]
            data <- rbind(data, d)
        }
    }else{
        data <- read.table(fname, header = FALSE, sep = ",", dec=".")[-c(1:nburn), 1:length(years)]
    }
    

    if(is.na(mask)){
        data[data == 0] <- NA
    }else{
        mask <- which(as.numeric(as.vector(mask))==1)
        data[,as.vector(unique(mask))] <- NA
    }

    if(dist == "lognorm"){
        pp <- post.pred.lognorm(data, var)
    }else{
        pp <- post.pred.negbin(data, var)
    }
    
    colnames(pp) <- years

    samples <- as_tibble(pp) %>%
                pivot_longer(everything(), names_to = "year", values_to = "data") %>%
                na.omit() %>%
                group_by(year) %>%
                median_qi(data, .width = c(0.05, 0.5, 0.95)) %>%
                print(n=150)

    return(samples) 
}

plot.survey.fits <- function(fits, survey.data, y.max, title, ylabel="", scale=1, cvs=TRUE){
    
    if(cvs){
        points <- geom_pointinterval(data=survey.data, aes(x=year, y=data/scale, ymin=lower/scale, ymax=upper/scale))
    }else{
        points <- geom_point(data=survey.data, aes(x=year, y=data/scale))
    }
    
    return(
        ggplot(fits) +
            geom_lineribbon(aes(x=year, y=data/scale, ymin=.lower/scale, ymax=.upper/scale, group=1), size=0)+
            scale_fill_brewer(palette = "Blues") +
            points +
            geom_line(data=survey.data, aes(x=year, y=data/scale, group=1))+
            scale_x_discrete("Year", breaks=seq(min(years), max(years), by=5))+
            labs(y=ylabel, title=title)+
            coord_cartesian(xlim=c(1, length(years)), ylim=c(0, y.max/scale))
    )
}

composite.raw.data <- function(model.dirs, nyr){
    #read.data.files(base.dir)$PWS_ASA.dat
    raw.data <- list(
        egg_se          = matrix(NA, nrow=nyr, ncol=length(model.dirs)),
        pwssc_hydro_se  = matrix(NA, nrow=nyr, ncol=length(model.dirs)),
        mdm             = matrix(NA, nrow=nyr, ncol=length(model.dirs)),
        egg             = matrix(NA, nrow=nyr, ncol=length(model.dirs)),
        adfg_hydro      = matrix(NA, nrow=nyr, ncol=length(model.dirs)),
        pwssc_hydro     = matrix(NA, nrow=nyr, ncol=length(model.dirs)),
        juvenile_survey = matrix(NA, nrow=nyr, ncol=length(model.dirs))
    )
    raw.data.names <- names(raw.data)

    for(i in 1:length(model.dirs)){
        d <- model.dirs[i]
        rd <- read.data.files(d)$PWS_ASA.dat
        rd <- rd[raw.data.names]

        for(n in names(rd)){
            rd[[n]][rd[[n]]==-9] <- NA
            raw.data[[n]][,i] <- rd[[n]]
        }
    }

    mean.raw.data <- list()
    for(n in raw.data.names){
        mean.raw.data[[n]] <- apply(raw.data[[n]], 1, median, na.rm=TRUE)
    }

    return(mean.raw.data)

}

start.year <- 1980
curr.year <- 2022
nyr.sim <- 35
years <- seq(start.year, curr.year+nyr.sim-1)
nyr <- length(years)

sims <- c(197, 649, 1017, 1094, 1144, 1787, 1998, 2078, 2214, 2241, 2255, 2386, 2512, 3169, 3709, 4288, 4716, 4775, 7251, 7915, 8004, 8388, 8462, 8634, 8789, 8904, 8935, 9204, 9260, 9716, 9725)
cr <- c("base")

model.dirs <- apply(
                as.matrix(expand.grid(sims, cr)), 
                1, 
                function(x) 
                    paste0(here::here("results"), "/", x[2], "/sim_", as.numeric(x[1]), "/year_", nyr.sim, "/model/")
            )

base.dir <- "/Users/jzahner/Desktop/Projects/basa/model/"
variances <- read.table(paste0(base.dir, "/mcmc_out/VARSReport.csv"), sep=",", header=FALSE)[-c(1:1),]
names(variances) <- c("mdm", "egg", "adfg.hydro", "pwssc.hydro")

add.vars <- apply(variances, 2, median)
raw.data <- composite.raw.data(model.dirs, nyr)
cvs <- list(
    mdm         = as.vector(rep(add.vars[1], nyr)),
    egg         = as.vector(sqrt(raw.data$egg_se[1:nyr]^2 + add.vars[2]^2)),
    adfg.hydro  = as.vector(rep(add.vars[3], nyr)),
    pwssc.hydro = as.vector(sqrt((raw.data$pwssc_hydro_se[1:nyr]^2) + (add.vars[4]^2)))
)

mdm.fname <- apply(as.matrix(model.dirs), 1, function(x) paste0(x, "/mcmc_out/MDM.csv"))
egg.fname <- apply(as.matrix(model.dirs), 1, function(x) paste0(x, "/mcmc_out/EGG.csv"))
adfg.hydro.fname <- apply(as.matrix(model.dirs), 1, function(x) paste0(x, "/mcmc_out/HYD_ADFG.csv"))
pwssc.hydro.fname <- apply(as.matrix(model.dirs), 1, function(x) paste0(x, "/mcmc_out/HYD_PWSSC.csv"))
# juv.schools.fname <- apply(as.matrix(model.dirs), 1, function(x) paste0(x, "/mcmc_out/juv_schools.csv"))
# overdisp.fname <- paste0(base.dir, "/mcmc_out/iterations.csv")

# juv.overdisp <- read.table(overdisp.fname,  header = TRUE, sep = ",")[-c(1:1), 20]
# juv.schools.mask <- c(raw.data$juvenile_survey == NaN)[1:nyr]

mdm.pp <- generate.post.pred(mdm.fname, variances$mdm, years)
egg.pp <- generate.post.pred(egg.fname, variances$egg, years)
adfg.hydro.pp <- generate.post.pred(adfg.hydro.fname, variances$adfg.hydro, years)
pwssc.hydro.pp <- generate.post.pred(pwssc.hydro.fname, variances$pwssc.hydro, years)
#juv.schools.pp <- generate.post.pred(juv.schools.fname, juv.overdisp, years, dist="negbin", mask=juv.schools.mask)


# Read in and format raw survey data for each annual survey
mdm.data            <- data.frame(year=as.character(years), data=raw.data$mdm[1:length(years)], lower=raw.data$mdm[1:length(years)]/calc.buck.cv(cvs$mdm), upper=raw.data$mdm[1:length(years)]*calc.buck.cv(cvs$mdm))
egg.data            <- data.frame(year=as.character(years), data=raw.data$egg[1:length(years)], lower=raw.data$egg[1:length(years)]/calc.buck.cv(cvs$egg), upper=raw.data$egg[1:length(years)]*calc.buck.cv(cvs$egg))
adfg.hydro.data     <- data.frame(year=as.character(years), data=raw.data$adfg_hydro[1:length(years)], lower=raw.data$adfg_hydro[1:length(years)]/calc.buck.cv(cvs$adfg.hydro), upper=raw.data$adfg_hydro[1:length(years)]*calc.buck.cv(cvs$adfg.hydro))
pwssc.hydro.data    <- data.frame(year=as.character(years), data=raw.data$pwssc_hydro[1:length(years)], lower=raw.data$pwssc_hydro[1:length(years)]/calc.buck.cv(cvs$pwssc.hydro), upper=raw.data$pwssc_hydro[1:length(years)]*calc.buck.cv(cvs$pwssc.hydro))
#aer.juvenile.data   <- data.frame(year=as.character(years), data=raw.data$juvenile_survey[1:length(years)])

data <- list(
    mdm = mdm.data,
    egg = egg.data,
    adfg.hydro = adfg.hydro.data,
    pwssc.hydro = pwssc.hydro.data#,
    #juvenile = aer.juvenile.data
)
for(i in 1:length(data)){
    missing <- data[[i]]$data == -9
    data[[i]]$data[missing] <- NA
    data[[i]]$lower[missing] <- NA
    data[[i]]$upper[missing] <- NA
}


mdm.fit.plot <- plot.survey.fits(mdm.pp, data$mdm, y.max=500, title="Mile Days of Milt")
egg.fit.plot <- plot.survey.fits(egg.pp, data$egg, y.max=21, title="Egg Deposition (trillions")
adfg.hydro.fit.plot <- plot.survey.fits(adfg.hydro.pp, data$adfg.hydro, y.max=175000, title="ADF&G Hydroacoustic Biomass (1000s tons)", scale=1000)
pwssc.hydro.fit.plot <- plot.survey.fits(pwssc.hydro.pp, data$pwssc.hydro, y.max=205000, title="PWSSC Hydroacoustic Biomass (1000s tons)", scale=1000)
#aer.juvenile.fit.plot <- plot.survey.fits(juv.schools.pp, data$juvenile, y.max=400000, title="Age 1 Schools", cvs=FALSE)

library(ggpubr)

ggarrange(
    mdm.fit.plot, egg.fit.plot, adfg.hydro.fit.plot, pwssc.hydro.fit.plot, #aer.juvenile.fit.plot,
    nrow=2,
    ncol=2,
    common.legend = TRUE, legend="right"
)


# years.in <- unique(juv.schools.pp$year)
# years.out <- setdiff(years, years.in)

# for(w in c(0.05, 0.5, 0.95)){
#     for(y in years.out){
#         d <- c(as.character(y), NA, NA, NA, w, "median", "qi")
#         juv.schools.pp <- rbind(juv.schools.pp, d)
#     }
# }

# juv.schools.pp[c("data", ".lower", ".upper", ".width")] <- sapply(juv.schools.pp[c("data", ".lower", ".upper", ".width")], as.numeric)
