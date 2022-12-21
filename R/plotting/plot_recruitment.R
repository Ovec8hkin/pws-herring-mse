source(file=paste0(here::here("R/plotting/", "compute_plot_products.R")))
source(file=paste0(here::here("R/plotting/", "plot_utils.R")))

start.year <- 1980
curr.year <- 2022
nyr.sim <- 35
years <- seq(start.year, curr.year+nyr.sim-1)
nyr <- length(years)

model.dir <- here::here("results/base/sim_1017/year_35/model/")

recruit.df <- compute.recruitment(model.dir, nyr, years)
# plot <- plot.recruitment.posterior(recruit.df, years)
# plot +
#      scale_fill_grey(start=0.8, end=0.6)

regimes.df <- data.frame(
                xmin=c(0, 13.25, 48.25, 63.25),
                ymin=rep(0, 4),
                ymax=rep(2000, 4),
                xmax=c(13, 48, 63, 77),
                c=as.factor(rep(c("red", "blue"), 2))
              )

regime.text <- data.frame(
    x = c(4, 17, 52, 68),
    y = rep(1900, 4),
    text = rep(c("High", "Low"), 2)
)

ggplot(recruit.df) +
    geom_lineribbon(aes(x=year, y=recruits, ymin=.lower, ymax=.upper, group=1), size=0.75)+
    geom_vline(xintercept=42, linetype="longdash")+
    scale_fill_grey(start=0.8, end=0.6)+
    geom_rect(
        data=regimes.df, 
        aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, color=c), 
        alpha=0.4, fill=NA, size=0.75
    )+
    geom_text(
        data = regime.text,
        aes(x=x, y=y, label=text),
        size=8
    )+
    scale_color_manual(values=c("blue", "red"))+
    scale_x_discrete("Year", breaks=seq(min(years), max(years), by=5))+
    scale_y_continuous("Age-3 Recruits (millions)", breaks=seq(0, 2000, by=500))+
    coord_cartesian(ylim=c(0, 2000), expand=c(0, 0))+
    ggtitle("Age-3 Recruitment")

ggsave("/Users/jzahner/Desktop/recruitment.eps", device="eps", dpi=320)
