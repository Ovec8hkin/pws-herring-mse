source(file=paste0(here::here("R/plotting/", "compute_plot_products.R")))
source(file=paste0(here::here("R/plotting/", "plot_utils.R")))

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

start.year <- 1980
curr.year <- 2023
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)
nyr <- length(years)

model.dir <- "/Users/jzahner/Desktop/Projects/basa/model/"#here::here("results/base/sim_649/year_0/model/")

sim.rec.devs <- generate.recruitment.deviates(30, 0998)$devs
sim.age0.recs <- exp(5.7276+sim.rec.devs)
sim.age3.recs <- sim.age0.recs*exp(-0.25)^3

sim.rec.devs.df <- data.frame(
  year=as.character(rep(curr.year:(curr.year+30-1), 2)),
  recruits=rep(sim.age3.recs, 2),
  .lower=0,
  .upper=0,
  .width=rep(c(0.50, 0.95), each=length(curr.year:(curr.year+30-1))),
  .point="median",
  .interval="qi"
)

recruit.df <- compute.recruitment(model.dir, nyr, years)
# plot <- plot.recruitment.posterior(recruit.df, years)
# plot +
#      scale_fill_grey(start=0.8, end=0.6)
recruit.df <- recruit.df %>% bind_rows(sim.rec.devs.df)

recruit.df <- recruit.df %>% 
                mutate(
                  regime=as.factor(ifelse(year < 1992 | (year > 2021 & year < 2037), "high", "low"))
                ) %>% 
                print(n=200)

years <- as.numeric(unique(recruit.df$year))

median.lines.df <- data.frame(
  xmin=c(0, 13, 43, 58),
  xmax=c(13, 43, 58, 73),
  y = c(
    recruit.df %>% filter(year > 1980 & year < 1993) %>% pull(recruits) %>% mean,
    recruit.df %>% filter(year > 1992 & year < 2023) %>% pull(recruits) %>% mean,
    recruit.df %>% filter(year > 2022 & year < 2038) %>% pull(recruits) %>% mean,
    recruit.df %>% filter(year > 2037 & year < 2053) %>% pull(recruits) %>% mean
  )
)

ggplot(recruit.df) +
    #geom_lineribbon(aes(x=year, y=recruits, ymin=.lower, ymax=.upper, group=1), size=0.75)+
    geom_line(data=recruit.df %>% filter(.width==0.5), aes(x=year, y=recruits, color=regime, group=1), size=0.75)+
    geom_vline(xintercept=43, linetype="longdash")+
    geom_hline(yintercept=exp(5.7276)*exp(-0.25)^3, linetype="dotted")+
    scale_fill_grey(start=0.8, end=0.6)+
    geom_segment(
      data=median.lines.df, 
      aes(x=xmin, xend=xmax, y=y, yend=y), 
      linetype="solid"
    )+
    scale_color_manual("Regime", values=c("red", "blue", "red", "blue"))+
    scale_x_discrete("Year", breaks=seq(min(years), max(years), by=10))+
    scale_y_continuous("Age-3 Recruits (millions)", breaks=seq(0, 2000, by=500))+
    coord_cartesian(ylim=c(0, 2000), expand=0, clip="off")+
    annotate(
      "segment", x = 43, y = 2000, xend = 77, yend = 2000, size=2,
      arrow = arrow(type = "closed", length = unit(0.02, "npc"))
    )+
    annotate("text", x=45, y=1850, label="Simulated", hjust=0, size=5)+
    ggtitle("Age-3 Recruitment")+
    theme_minimal()+
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(),
      plot.title = element_text(size=16),
      axis.title = element_text(size=12),
      axis.text = element_text(size=10),
      axis.ticks = element_line()
    )+
    guides(color="none")

#ggsave("/Users/jzahner/Desktop/recruitment.jpg", dpi=320)
#ggsave(file.path(here::here(), "figures", "publication", "Fig2_recruitment.jpg"), dpi=300, width=170, height=105, units="mm")

ggsave(file.path(here::here(), "figures", "publication", "Fig2_recruitment.pdf"), dpi=300, width=170, height=105, units="mm") 
