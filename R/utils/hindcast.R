source(file=paste0(here::here("R/utils/"),            "fun_read_dat.R"))
source(file=paste0(here::here("R/operating_model/"),  "fun_obsm.R"))
source(file=paste0(here::here("R/estimation_model/",  "run_basa.R")))


hindcast <- function(dir){
    # Set age-composition sample sizes
    sample.sizes <- list(
        seac = 500,
        spac = 1500
    )

    # Read in data that is know without error from the most recent assessment
    data <- read.data.files("/Users/jzahner/Desktop/Projects/basa/model/")
    waa <- data$PWS_ASA.dat$waa
    fec <- data$PWS_ASA.dat$fecundity
    perc.female <- data$PWS_ASA.dat$perc.female

    juv.overdisp <- 2.17084637519

    # Read in survey estimates from most recent assessment
    survey.indices <- read.survey.estimates("/Users/jzahner/Desktop/Projects/basa/model/mcmc_out/")
    survey.indices$vhsv <- data$PWS_ASA_disease.dat$vhsv_age_prevalence     # These are assumed known without error for now
    survey.indices$ich <- data$PWS_ASA_disease.dat$ich_age_prevalence
    names(survey.indices) <- c("mdm", "egg", "PWSSC.hydro", "ADFG.hydro", "juv.schools", "spawn.age.comp", "seine.age.comp", "vhsv.antibody", "ich.antibody")

    dat.files <- data

    for(y in 1:data$PWS_ASA.dat$nyr){
        obs <- fun_obsm(survey.indices, waa, fec, perc.female, juv.overdisp, sample.sizes, y, 9716,
                        survey.controls = list(mdm=TRUE, egg=TRUE, pwssc.hydro=TRUE, adfg.hydro=TRUE, spac=TRUE, seac=TRUE, juv=TRUE, vhsv=TRUE, ich=TRUE))
        dat.files$PWS_ASA.dat$mdm[y]                <- ifelse(dat.files$PWS_ASA.dat$mdm[y]              == -9, -9, obs$mdm)
        dat.files$PWS_ASA.dat$egg[y]                <- ifelse(dat.files$PWS_ASA.dat$egg[y]              == -9, -9, obs$egg)
        dat.files$PWS_ASA.dat$egg_se[y]             <- ifelse(dat.files$PWS_ASA.dat$egg_se[y]           == -9, -9, obs$egg.sd)
        dat.files$PWS_ASA.dat$adfg_hydro[y]         <- ifelse(dat.files$PWS_ASA.dat$adfg_hydro[y]       == -9, -9, obs$adfg.hydro)
        dat.files$PWS_ASA.dat$pwssc_hydro[y]        <- ifelse(dat.files$PWS_ASA.dat$pwssc_hydro[y]      == -9, -9, obs$pwssc.hydro)
        dat.files$PWS_ASA.dat$pwssc_hydro_se[y]     <- ifelse(dat.files$PWS_ASA.dat$pwssc_hydro_se[y]   == -9, -9, obs$pwssc.hydro.sd)
        dat.files$PWS_ASA.dat$seine_age_comp[y, ]   <- ifelse(dat.files$PWS_ASA.dat$seine_age_comp[y,]  == rep(-9, 10), rep(-9, 10), obs$seine.age.comp)
        dat.files$PWS_ASA.dat$spawn_age_comp[y, ]   <- ifelse(dat.files$PWS_ASA.dat$spawn_age_comp[y,]  == rep(-9, 10), rep(-9, 10), obs$spawn.age.comp)
        dat.files$PWS_ASA.dat$juvenile_survey[y]    <- ifelse(dat.files$PWS_ASA.dat$juvenile_survey[y]  == -9, -9, obs$aerial.juv)

        dat.files$agecomp_samp_sizes.txt$seine_sample_size[y] <- ifelse(dat.files$agecomp_samp_sizes.txt$seine_sample_size[y] == 0, 0, sample.sizes$seac)
        dat.files$agecomp_samp_sizes.txt$spawn_sample_size[y] <- ifelse(dat.files$agecomp_samp_sizes.txt$spawn_sample_size[y] == 0, 0, sample.sizes$spac)

    }

    setwd(dir)

    fnames <- c("PWS_ASA.dat", "PWS_ASA(ESS).ctl", "PWS_ASA(covariate).ctl", "agecomp_samp_sizes.txt", "PWS_ASA_disease.dat")
    for(d in 1:length(dat.files)){
        write.data(dat.files[[d]], fnames[d])
    }

    convergence.diags <- run.basa.adnuts(dir, 9716, max.duration = 10)
}

# dir <- here::here("results/test/year_0/model/")
# hind <- hindcast(dir)

# # Compare assessment to underlying model truth
# ass.biomass <- read.csv(paste0(dir, "mcmc_out/PFRBiomass.csv"), header=FALSE)
# true.biomass <- read.csv("/Users/jzahner/Desktop/Projects/basa/model/mcmc_out/PFRBiomass.csv", header=FALSE)
# ass.biomass.quants <- as.data.frame(t(apply(ass.biomass, 2, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))))
# true.biomass.quants <- as.data.frame(t(apply(true.biomass, 2, quantile, c(0.025, 0.5, 0.975))))

# ggplot(ass.biomass.quants)+
#     geom_ribbon(aes(x=1:43, ymin=`2.5%`, ymax=`97.5%`), fill=rgb(0.75, 0.75, 0.75))+
#     geom_ribbon(aes(x=1:43, ymin=`25%`, ymax=`75%`), fill=rgb(0.5, 0.5, 0.5))+
#     geom_line(aes(x=1:43, y=`50%`))+
#     geom_point(data=true.biomass.quants, aes(x=1:43, y=`50%`), color="white")+
#     labs(x="Simulation Year", y="Biomass", title="Full Simulation")+
#     ylim(c(0, 200000))
    
# # Compute assessment bias
# ass.biomass.median <- apply(ass.biomass, 2, median)
# true.biomass.median <- apply(true.biomass, 2, median)

# ass.bias <- (ass.biomass.median - true.biomass.median)/true.biomass.median
# mean(ass.bias)


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
