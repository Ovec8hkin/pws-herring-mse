b.lim <- 19958
b.tar <- 38555

f.min <- 0
f.max <- 0.2

biomass <- seq(0, 60000, by=50)
rel.ssb.change <- seq(0.5, 2.0, 0.005)
evenness <- seq(0.0, 1.0, by=0.005)

compute.f.thresh <- function(biomass, b.lim=19956, b.tar=38556, f.max=0.2, scaling="linear"){
    f <- 0
    if(biomass <= b.lim){
        f <- 0
    }else if(biomass > b.lim && biomass <= b.tar){
        if(scaling == "linear"){
            f <- (biomass - b.lim) * f.max / (b.tar - b.lim)
        }else if(scaling == "logistic"){
            f <- f.max/(1+exp(-0.0007*(biomass - (b.lim+((b.tar-b.lim)/2)))))
        }else if(scaling == "exponential"){
            f <- (f.max*2)/(1+exp(-0.0004*(biomass - b.tar)))
        }else if(scaling == "logarithmic"){
            f <- (f.max*2)/(1+exp(-0.0004*(biomass - b.lim)))-f.max
        }else{
            f <- NA
        }
    }else{
        f <- f.max
    }

    return(f)
}

compute.f.gradient <- function(biomass, rel.ssb.change, b.lim=19958){
    f.def <- compute.f.thresh(biomass)
    f.grad <- f.def*rel.ssb.change
    p=1.0
    return(p*f.def + (1-p)*f.grad)
}

compute.f.age.struct <- function(biomass, evenness, b.lim=19956){
    f <- compute.f.thresh(biomass)
    e <- compute.f.thresh(evenness, b.lim=0.4, b.tar=0.8, f.max=0.5)

    return(f*(e + 0.5))
}

# image.plot(x=biomass, y=rel.ssb.change, z=fs.gradient)
# contour(x=biomass, y=rel.ssb.change, z=fs.gradient, levels=c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5), add=TRUE)

fs.gradient <- apply(matrix(biomass), 1, function(b) apply(matrix(rel.ssb.change), 1, function(r) compute.f.gradient(b, r)))
fs.gradient <- t(as.matrix(fs.gradient, nrow=length(rel.ssb.change)))
catch.gradient <- apply(fs.gradient, 2, function(x) x*biomass)

gradient.cr.df <- data.frame(biomass=rep(biomass, each=length(rel.ssb.change)), ssb.change=rel.ssb.change, fs=as.vector(t(fs.gradient)), catch=as.vector(t(catch.gradient)))

fs.plot <- ggplot(gradient.cr.df, aes(x=biomass, y=ssb.change, fill=fs, z=fs))+
    geom_raster(alpha=0.9)+
    geom_contour(breaks=c(0.01, 0.05, 0.1, 0.15, 0.20, 0.25, 0.40, 0.50), color="black", size=1)+
    geom_label_contour(breaks=c(0.01, 0.05, 0.1, 0.15, 0.20, 0.25), skip=0, label.placer=label_placer_fraction(0.5))+
    scale_fill_viridis_c(option="turbo", name="Exploitation Rate")+
    scale_y_continuous(expand=c(0, 0), name="3-year Biomass Gradient")+
    scale_x_continuous(expand=c(0, 0), name="Pre-Fishery Biomass (mt)", breaks=seq(0, 60000, 10000), labels=seq(0, 60, 10))+
    theme(legend.position="bottom", legend.key.width=unit(1.5, "cm"))

catches.plot <- ggplot(gradient.cr.df, aes(x=biomass, y=ssb.change, fill=catch, z=catch))+
    geom_raster(alpha=0.9)+
    geom_contour(breaks=c(1, 2000, 5000, 7500, 10000, 15000), color="black", size=1)+
    geom_label_contour(breaks=c(1, 2000, 5000, 7500, 10000, 15000), skip=0,label.placer=label_placer_fraction(0.5))+
    scale_fill_viridis_c(option="turbo", name="Catch (mt)")+
    scale_y_continuous(expand=c(0, 0), name="3-year Biomass Gradient")+
    scale_x_continuous(expand=c(0, 0), name="Pre-Fishery Biomass (mt)", breaks=seq(0, 60000, 10000), labels=seq(0, 60, 10))+
    theme(legend.position="bottom", legend.key.width=unit(1.5, "cm"))

ggarrange(fs.plot, catches.plot+rremove("ylab"))






par(mfrow=c(1, 2), mar=c(5,6,4,2.5)+0.1)
image.plot(x=biomass, y=rel.ssb.change, z=fs.gradient)
contour(x=biomass, y=rel.ssb.change, z=fs.gradient, levels=c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5), add=TRUE)
image.plot(x=biomass, y=rel.ssb.change, z=catch.gradient)
contour(x=biomass, y=rel.ssb.change, z=catch.gradient, levels=c(2000, 5000, 7500, 10000, 12000), add=TRUE)


fs.age.struct <- apply(matrix(biomass), 1, function(b) apply(matrix(evenness), 1, function(e) compute.f.age.struct(b, e)))
fs.age.struct <- t(as.matrix(fs.age.struct, nrow=length(evenness)))
catch.age.struct <- apply(fs.age.struct, 2, function(x) x*biomass)

age.struct.cr.df <- data.frame(biomass=rep(biomass, each=length(evenness)), evenness=evenness, fs=as.vector(t(fs.age.struct)), catch=as.vector(t(catch.age.struct)))

fs.plot <- ggplot(age.struct.cr.df, aes(x=biomass, y=evenness, fill=fs, z=fs))+
    geom_raster(alpha=0.9)+
    geom_contour(breaks=c(0.01, 0.05, 0.1, 0.15, 0.20), color="black", size=1)+
    geom_label_contour(breaks=c(0.01, 0.05, 0.1, 0.15, 0.20), skip=0, label.placer=label_placer_fraction(0.5))+
    scale_fill_viridis_c(option="turbo", name="Exploitation Rate")+
    scale_y_continuous(expand=c(0, 0), name="Age Structure Evenness Index")+
    scale_x_continuous(expand=c(0, 0), name="Pre-Fishery Biomass (mt)", breaks=seq(0, 60000, 10000), labels=seq(0, 60, 10))+
    theme(legend.position="bottom", legend.key.width=unit(1.5, "cm"))

catches.plot <- ggplot(age.struct.cr.df, aes(x=biomass, y=evenness, fill=catch, z=catch))+
    geom_raster(alpha=0.9)+
    geom_contour(breaks=c(1, 2000, 5000, 7500, 10000), color="black", size=1)+
    geom_label_contour(breaks=c(1, 2000, 5000, 7500, 10000), skip=0,label.placer=label_placer_fraction(0.5))+
    scale_fill_viridis_c(option="turbo", name="Catch (mt)")+
    scale_y_continuous(expand=c(0, 0), name="Age Structure Evenness Index")+
    scale_x_continuous(expand=c(0, 0), name="Pre-Fishery Biomass (mt)", breaks=seq(0, 60000, 10000), labels=seq(0, 60, 10))+
    theme(legend.position="bottom", legend.key.width=unit(1.5, "cm"))

ggarrange(fs.plot, catches.plot+rremove("ylab"))

# image.plot(x=biomass, y=evenness, z=fs.age.struct, horizontal=TRUE)
# contour(x=biomass, y=evenness, z=fs.age.struct, levels=c(-0.1, 0.0, 0.1, 0.2), add=TRUE)

# image.plot(x=biomass, y=evenness, z=catch.age.struct)
# contour(x=biomass, y=evenness, z=catch.age.struct, levels=c(2000, 5000, 7500, 10000, 12000), add=TRUE)



fs.default      <- apply(matrix(biomass), 1, compute.f.thresh, scaling="linear")
fs.logistic     <- apply(matrix(biomass), 1, compute.f.thresh, scaling="logistic")
fs.exponential  <- apply(matrix(biomass), 1, compute.f.thresh, scaling="exponential")
fs.logarithmic  <- apply(matrix(biomass), 1, compute.f.thresh, scaling="logarithmic")
#fs.as.dynamic   <- apply(matrix(evenness), 1, compute.f.thresh, b.lim=0.60, b.tar=0.75, f.max=0.2, scaling="linear")
fs.low.thresh   <- apply(matrix(biomass), 1, compute.f.thresh, b.lim=10000)
fs.high.thresh  <- apply(matrix(biomass), 1, compute.f.thresh, b.lim=30000)
fs.high.f       <- apply(matrix(biomass), 1, compute.f.thresh, f.max=0.4)
fs.low.f        <- apply(matrix(biomass), 1, compute.f.thresh, f.max=0.1)

catches.default     <- fs.default*biomass
catches.logistic    <- fs.logistic*biomass
catches.exponential <- fs.exponential*biomass
catches.logarithmic <- fs.logarithmic*biomass
catches.low.thresh  <- fs.low.thresh*biomass
catches.high.thresh <- fs.high.thresh*biomass
catches.high.f      <- fs.high.f*biomass
catches.low.f       <- fs.low.f*biomass


attach(mtcars)
layout(matrix(c(0,1,0,2,3,4,5,6,7), 3, 3, byrow = TRUE))

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

plot(biomass, fs.default, type="l", ylim=c(0, 0.5), col="black", lwd=6,
     ylab="Exploitation Rate (µ)", xlab="Biomass (1000 mt)",
     yaxs="i", xaxs="i", ax=FALSE, cex.lab=1.75)
abline(v=b.lim, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
abline(v=b.tar, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
text(x=b.lim-8500, y=0.4, expression(bold("B"["limit"])), cex=2.00)
text(x=b.tar+8500, y=0.4, expression(bold("B"["target"])), cex=2.00)
title(main="Default", cex.main=2.25)
axis(1, at=seq(0, 60000, 10000), labels=seq(0, 60, 10), cex.axis=1.75)
axis(2, at=seq(0, 0.5, 0.1), labels=seq(0, 0.5, 0.1), cex.axis=1.75)

plot(biomass, fs.logistic, type="l", ylim=c(0, 0.5), col="black", lwd=6,
     ylab="Exploitation Rate (µ)", xlab="Biomass (1000 mt)",
     yaxs="i", xaxs="i", ax=FALSE, cex.lab=1.75)
abline(v=b.lim, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
abline(v=b.tar, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
text(x=b.lim-8500, y=0.4, expression(bold("B"["limit"])), cex=2.00)
text(x=b.tar+8500, y=0.4, expression(bold("B"["target"])), cex=2.00)
title(main="Logistic", cex.main=2.25)
axis(1, at=seq(0, 60000, 10000), labels=seq(0, 60, 10), cex.axis=1.75)
axis(2, at=seq(0, 0.5, 0.1), labels=seq(0, 0.5, 0.1), cex.axis=1.75)

plot(biomass, fs.exponential, type="l", ylim=c(0, 0.5), col="black", lwd=6,
     ylab="Exploitation Rate (µ)", xlab="Biomass (1000 mt)",
     yaxs="i", xaxs="i", ax=FALSE, cex.lab=1.75)
abline(v=b.lim, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
abline(v=b.tar, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
text(x=b.lim-8500, y=0.4, expression(bold("B"["limit"])), cex=2.00)
text(x=b.tar+8500, y=0.4, expression(bold("B"["target"])), cex=2.00)
title(main="Exponential", cex.main=2.25)
axis(1, at=seq(0, 60000, 10000), labels=seq(0, 60, 10), cex.axis=1.75)
axis(2, at=seq(0, 0.5, 0.1), labels=seq(0, 0.5, 0.1), cex.axis=1.75)

plot(biomass, fs.logarithmic, type="l", ylim=c(0, 0.5), col="black", lwd=6,
     ylab="Exploitation Rate (µ)", xlab="Biomass (1000 mt)",
     yaxs="i", xaxs="i", ax=FALSE, cex.lab=1.75)
abline(v=b.lim, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
abline(v=b.tar, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
text(x=b.lim-8500, y=0.4, expression(bold("B"["limit"])), cex=2.00)
text(x=b.tar+8500, y=0.4, expression(bold("B"["target"])), cex=2.00)
title(main="Logarithmic", cex.main=2.25)
axis(1, at=seq(0, 60000, 10000), labels=seq(0, 60, 10), cex.axis=1.75)
axis(2, at=seq(0, 0.5, 0.1), labels=seq(0, 0.5, 0.1), cex.axis=1.75)

plot(biomass, fs.low.thresh, type="l", ylim=c(0, 0.5), col="black", lwd=6,
     ylab="Exploitation Rate (µ)", xlab="Biomass (1000 mt)",
     yaxs="i", xaxs="i", ax=FALSE, cex.lab=1.75)
abline(v=10000, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
abline(v=b.tar, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
text(x=10000+8500, y=0.4, expression(bold("B"["limit"])), cex=2.00)
text(x=b.tar+8500, y=0.4, expression(bold("B"["target"])), cex=2.00)
title(main="Low Threshold", cex.main=2.25)
axis(1, at=seq(0, 60000, 10000), labels=seq(0, 60, 10), cex.axis=1.75)
axis(2, at=seq(0, 0.5, 0.1), labels=seq(0, 0.5, 0.1), cex.axis=1.75)

plot(biomass, fs.high.thresh, type="l", ylim=c(0, 0.5), col="black", lwd=6,
     ylab="Exploitation Rate (µ)", xlab="Biomass (1000 mt)",
     yaxs="i", xaxs="i", ax=FALSE, cex.lab=1.75)
abline(v=30000, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
abline(v=b.tar, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
text(x=30000-8500, y=0.4, expression(bold("B"["limit"])), cex=2.00)
text(x=b.tar+8500, y=0.4, expression(bold("B"["target"])), cex=2.00)
title(main="High Threshold", cex.main=2.25)
axis(1, at=seq(0, 60000, 10000), labels=seq(0, 60, 10), cex.axis=1.75)
axis(2, at=seq(0, 0.5, 0.1), labels=seq(0, 0.5, 0.1), cex.axis=1.75)

plot(evenness, fs.as.dynamic, type="l", ylim=c(0, 0.5), col="black", lwd=6,
     ylab="Exploitation Rate (µ)", xlab=paste("Shannon Evenness", expression("H^'/ln(S)")),
     yaxs="i", xaxs="i", ax=FALSE, cex.lab=1.75)
abline(v=0.6, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
abline(v=0.75, col=rgb(0.5, 0.5, 0.5), lwd=4, lty=2)
text(x=0.6-0.15, y=0.4, expression(bold("B"["limit"])), cex=2.00)
text(x=0.75+0.15, y=0.4, expression(bold("B"["target"])), cex=2.00)
title(main="Dynamic Age Structure", cex.main=2.25)
axis(1, at=seq(0, 1.0, 0.2), labels=seq(0, 1.0, 0.2), cex.axis=1.75)
axis(2, at=seq(0, 0.5, 0.1), labels=seq(0, 0.5, 0.1), cex.axis=1.75)
