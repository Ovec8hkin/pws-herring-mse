source(paste0(here::here("R/operating_model/control_rules/"), "cr_threshold.R"))
source(paste0(here::here("R/operating_model/control_rules/"), "cr_constant_f.R"))
source(paste0(here::here("R/operating_model/control_rules/"), "cr_gradient.R"))
source(paste0(here::here("R/operating_model/control_rules/"), "cr_agestructure.R"))
source(paste0(here::here("R/plotting/plot_util_vals.R")))

b.lim <- 19958
b.tar <- 38555

f.min <- 0
f.max <- 0.2

biomasses <- seq(0, 100000, by=50)
rel.ssb.change <- seq(0.0, 4.0, 0.005)
evenness <- seq(0.0, 1.0, by=0.005)

############# Gradient Control Rule #############
fs.grad <- apply(matrix(biomasses), 1, function(b) apply(matrix(rel.ssb.change), 1, function(r) hcr.gradient(b, r)))
fs.grad <- t(as.matrix(fs.grad, nrow=length(rel.ssb.change)))
catches.grad <- apply(fs.grad, 2, function(x) x*biomasses)

gradient.plot <- plot.hcr.gradient(biomasses, rel.ssb.change, fs.grad, catches.grad, f.only=TRUE)
gradient.plot$p+
     scale_fill_gradient(low="white", high="red", na.value = "transparent", name="Harvest Rate")+
     scale_x_continuous("Pre-Fishery Biomass (mt)", breaks=seq(0, 100000, by=10000), labels=seq(0, 100, 10))+
     coord_cartesian(xlim=c(0, 100000), expand=c(0, 0))+
     theme(
          legend.position="right",
          legend.direction="vertical",
          legend.key.height=unit(1.5, "cm")
     )
############# Age Structured Control Rule #############
fs.as <- apply(matrix(biomasses), 1, function(b) apply(matrix(evenness), 1, function(e) hcr.agestructure(b, c(1), e=e)))
fs.as <- t(as.matrix(fs.as, nrow=length(evenness)))
catches.as <- apply(fs.as, 2, function(x) x*biomasses)

as.plot <- plot.hcr.age.structure(biomasses, evenness, fs.as, catches.as, f.only=TRUE)
as.plot$p+
     scale_fill_gradient(low="white", high="red", na.value = "transparent", name="Harvest Rate")+
     scale_x_continuous("Pre-Fishery Biomass (mt)", breaks=seq(0, 100000, by=10000), labels=seq(0, 100, 10))+
     coord_cartesian(xlim=c(0, 100000), expand=c(0, 0))+
     theme(
          legend.position="right",
          legend.direction="vertical",
          legend.key.height=unit(1.5, "cm")
     )

############# Simple Control Rule #############
fs.default      <- apply(matrix(biomasses), 1, hcr.threshold.linear)
fs.three.step   <- apply(matrix(biomasses), 1, hcr.threshold.multi, options=list(lower.threshold=19958, middle.threshold=38855, upper.threshold=60000, min.harvest=0.0, mid.harvest=0.2, max.harvest=0.6))
fs.low.thresh   <- apply(matrix(biomasses), 1, hcr.threshold.linear, options=list(lower.threshold=10000))
fs.high.thresh  <- apply(matrix(biomasses), 1, hcr.threshold.linear, options=list(lower.threshold=30000))
fs.high.f       <- apply(matrix(biomasses), 1, hcr.threshold.linear, options=list(max.harvest=0.4))
fs.low.f        <- apply(matrix(biomasses), 1, hcr.threshold.linear, options=list(max.harvest=0.1))
fs.constant.00  <- apply(matrix(biomasses), 1, hcr.constant.f, f.rate=0.0)
fs.big.fish     <- apply(matrix(biomasses), 1, hcr.threshold.linear)

fs.df <- data.frame(biomass=biomasses, 
                    base=fs.default, 
                    low.biomass=fs.low.thresh,
                    high.biomass=fs.high.thresh,
                    high.harvest=fs.high.f,
                    low.harvest=fs.low.f,
                    three.step.thresh=fs.three.step,
                    big.fish=fs.big.fish,
                    constant.f.00=fs.constant.00)

catches.default     <- fs.default*biomass
catches.logistic    <- fs.logistic*biomass
catches.exponential <- fs.exponential*biomass
catches.logarithmic <- fs.logarithmic*biomass
catches.low.thresh  <- fs.low.thresh*biomass
catches.high.thresh <- fs.high.thresh*biomass
catches.high.f      <- fs.high.f*biomass
catches.low.f       <- fs.low.f*biomass

par(mfrow=c(1, 2), mar=c(5,6,4,2.5)+0.1)
plot.control.rule(fs.default, catches.default, b.lim=19958, b.tar=38555, name="Default")
mtext("Default", side = 3, line = -3, outer=TRUE, cex=2.25)

par(mfrow=c(1, 2), mar=c(5,6,4,2.5)+0.1)
plot.control.rule(fs.logistic, catches.logistic, b.lim=19958, b.tar=38555, name="Default")
mtext("Logistic", side = 3, line = -3, outer=TRUE, cex=2.25)

par(mfrow=c(1, 2), mar=c(5,6,4,2.5)+0.1)
plot.control.rule(fs.exponential, catches.exponential, b.lim=19958, b.tar=38555, name="Default")
mtext("Exponential", side = 3, line = -3, outer=TRUE, cex=2.25)

par(mfrow=c(1, 2), mar=c(5,6,4,2.5)+0.1)
plot.control.rule(fs.logarithmic, catches.logarithmic, b.lim=19958, b.tar=38555, name="Default")
mtext("Logarithmic", side = 3, line = -3, outer=TRUE, cex=2.25)

par(mfrow=c(1, 2), mar=c(5,6,4,2.5)+0.1)
plot.control.rule(fs.low.thresh, catches.low.thresh, b.lim=10000, b.tar=38555, name="Default")
mtext("Low Threshold", side = 3, line = -3, outer=TRUE, cex=2.25)

par(mfrow=c(1, 2), mar=c(5,6,4,2.5)+0.1)
plot.control.rule(fs.high.thresh, catches.high.thresh, b.lim=30000, b.tar=38555, name="Default")
mtext("High Threshold", side = 3, line = -3, outer=TRUE, cex=2.25)

par(mfrow=c(1, 2), mar=c(5,6,4,2.5)+0.1)
plot.control.rule(fs.high.f, catches.high.f, b.lim=19958, b.tar=38555, name="Default")
mtext("High F", side = 3, line = -3, outer=TRUE, cex=2.25)

par(mfrow=c(1, 2), mar=c(5,6,4,2.5)+0.1)
plot.control.rule(fs.low.f, catches.low.f, b.lim=19958, b.tar=38555, name="Default")
mtext("Low F", side = 3, line = -3, outer=TRUE, cex=2.25)



#e <- compute.f.thresh(evenness, b.lim=0.4, b.tar=0.8, f.max=0.5)
e <- apply(matrix(evenness), 1, compute.f.thresh, b.lim=0.4, b.tar=0.8, f.max=0.5, scaling="linear")+0.5

par(mfrow=c(1, 2), mar=c(5,6,4,2.5)+0.1)
plot(biomass, fs.default, type="l", ylim=c(0, 1.0), col="black", lwd=6,
         ylab="Exploitation Rate (µ)", xlab="Biomass (1000 mt)",
         yaxs="i", xaxs="i", ax=FALSE, cex.lab=1.75)
abline(v=b.lim, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
abline(v=b.tar, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
text(x=b.lim-8500, y=0.5, expression(bold("B"["limit"])), cex=2.00)
text(x=b.tar+8500, y=0.5, expression(bold("B"["target"])), cex=2.00)
#title(main=name, cex.main=2.25)
axis(1, at=seq(0, 60000, 10000), labels=seq(0, 60, 10), cex.axis=1.75)
axis(2, at=seq(0, 1.0, 0.1), labels=seq(0, 1.0, 0.1), cex.axis=1.75)

plot(evenness, e, type="l", ylim=c(0, 1.0), col="black", lwd=6,
         ylab="Evenness Index", xlab="Shannon-Weiner Evenness (J')",
         yaxs="i", xaxs="i", ax=FALSE, cex.lab=1.75)
abline(v=0.4, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
abline(v=0.8, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
text(x=0.4-0.1, y=0.8, expression(bold("J"["limit"])), cex=2.00)
text(x=0.8+0.1, y=0.8, expression(bold("J"["target"])), cex=2.00)
#title(main=name, cex.main=2.25)
axis(1, at=seq(0, 1.0, 0.1), labels=seq(0, 1.0, 0.1), cex.axis=1.75)
axis(2, at=seq(0, 1.0, 0.1), labels=seq(0, 1.0, 0.1), cex.axis=1.75)
mtext("Age Structure CR Components", side = 3, line = -3, outer=TRUE, cex=2.25)


plot.control.rule <- function(exp.rates, catches, b.lim, b.tar, name){
    plot(biomass, exp.rates, type="l", ylim=c(0, 1.0), col="black", lwd=6,
         ylab="Exploitation Rate (µ)", xlab="Biomass (1000 mt)",
         yaxs="i", xaxs="i", ax=FALSE, cex.lab=1.75)
    abline(v=b.lim, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
    abline(v=b.tar, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
    text(x=b.lim-8500, y=0.5, expression(bold("B"["limit"])), cex=2.00)
    text(x=b.tar+8500, y=0.5, expression(bold("B"["target"])), cex=2.00)
    #title(main=name, cex.main=2.25)
    axis(1, at=seq(0, 60000, 10000), labels=seq(0, 60, 10), cex.axis=1.75)
    axis(2, at=seq(0, 1.0, 0.1), labels=seq(0, 1.0, 0.1), cex.axis=1.75)

    plot(biomass, catches, type="l", ylim=c(0, 30000), col="black", lwd=6,
         ylab="Catch (1000 mt)", xlab="Biomass (1000 mt)",
         yaxs="i", xaxs="i", ax=FALSE, cex.lab=1.75)
    abline(v=b.lim, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
    abline(v=b.tar, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
    text(x=b.lim-8500, y=26000, expression(bold("B"["limit"])), cex=2.00)
    text(x=b.tar+8500, y=26000, expression(bold("B"["target"])), cex=2.00)
   # title(main=name, cex.main=2.25)
    axis(1, at=seq(0, 60000, 10000), labels=seq(0, 60, 10), cex.axis=1.75)
    axis(2, at=seq(0, 30000, 5000), labels=seq(0, 30, 5), cex.axis=1.75)

}

library(tidyverse)

fs.df <- fs.df %>% pivot_longer(
     !biomass,
     names_to="control.rule",
     values_to="harvest.rate"
)

fs.df$control.rule <- factor(fs.df$control.rule,
                              levels=c("base", "low.harvest", "high.harvest", "low.biomass", "high.biomass", "three.step.thresh", "big.fish", "constant.f.00"),
                              labels=c("Default", "Low Harvest", "High Harvest", "Low Threshold", "High Threshold", "Three Step Threshold", "Big Fish (Fish > 130g)", "No Fishing"))

simple.cr.plot <- ggplot(fs.df, aes(x=biomass, y=harvest.rate, color=control.rule))+
     geom_line(size=1.5)+
     #geom_vline(aes(xintercept=20000))+
     scale_color_manual(values=as.vector(hcr.colors))+
     scale_y_continuous(limits=c(-0.005, 0.70), expand=c(0, 0), breaks=c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), labels=c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))+
     scale_x_continuous(breaks=seq(0, 100000, by=20000), labels=seq(0, 100, by=20), expand=c(0, 0))+
     facet_wrap(~control.rule)+
     labs(x="Biomass (1000 mt)", y="Harvest Rate", color="Control Rule", title="Simple Threshold Rules")+
     theme(
          legend.position="none",
          panel.grid.minor=element_blank(),
          panel.spacing.x=unit(0.5, "cm"),
          plot.margin = unit(c(0, 30, 0, 0), "pt"),
     )

grad.plot.sing <- gradient.plot$p+
     scale_fill_gradient(low="white", high="red", na.value = "transparent", name="Harvest Rate")+
     scale_x_continuous("Pre-Fishery Biomass (1000 mt)", expand=c(0, 0), breaks=seq(0, 100000, 10000), labels=seq(0, 100, 10))+
     labs(title="Gradient Rule")+
     theme(
          legend.position="right",
          legend.direction="vertical",
          legend.key.width=unit(1.0, "cm"),
          legend.key.height=unit(0.75, "cm")
     )

as.plot.sing <- as.plot$p+
     scale_fill_gradient(low="white", high="red", na.value = "transparent", name="Harvest Rate")+
     scale_x_continuous(expand=c(0, 0), name="", breaks=seq(0, 100000, 10000), labels=seq(0, 100, 10))+
     labs(title="Evenness Rule")+
     theme(
          legend.position="right",
          legend.direction="vertical",
          legend.key.width=unit(1.0, "cm"),
          legend.key.height=unit(0.75, "cm")
     )

library(patchwork)

(tag_facet(simple.cr.plot) + (as.plot.sing / grad.plot.sing)) +
     plot_annotation(
          title="Harvest Control Rules",
          tag_levels="A"
     )

ggsave("/Users/jzahner/Desktop/hcrs.eps", device="eps", dpi=320)


tag_facet <- function(p, open = "(", close = ")", tag_pool = 1:100, x = -Inf, y = Inf, 
    hjust = -0.5, vjust = 2, fontface = 1, family = "", ...) {

    gb <- ggplot_build(p)
    lay <- gb$layout$layout
    tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
    p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
        vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}
