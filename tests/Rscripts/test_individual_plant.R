# R script to test:
dat = read.delim("~/codes/Plant-FATE/assim1.txt")
dat$leaf_area = dat$crown_area * dat$lai
dat$heartwood_fraction = 1-dat$sapwood_fraction

# par(mfrow=c(4,2), mar=c(4,4,1,1), oma=c(1,1,1,1))

# plot(dat$height~dat$diameter, type="l", ylab="Height", xlab="Diameter")
# plot(dat$crown_area~I(dat$height*dat$diameter), type="l", ylab="Crown area", xlab="DH")
# plot(I(dat$total_rep/dat$total_prod)~dat$height, type="l", ylab="Frac alloc to\nreproduction", xlab="Height")
# plot(I(dat$total_rep/dat$total_prod)~dat$diameter, type="l", ylab="Frac alloc to\nreproduction", xlab="Diameter")

par(mfrow=c(4,4), mar=c(4,8,1,1), oma=c(1,1,1,1), mgp=c(4,1,0))

plot(dat$height~dat$i, ylab="Height", xlab="Year", type="l", lwd=2)
plot(dat$diameter~dat$i, ylab="Diameter", xlab="Year", col="brown", type="l", lwd=2)

matplot(y=cbind(dat$total_mass,
                dat$total_rep,
                dat$litter_mass,
                dat$total_mass+dat$total_rep+dat$litter_mass,
                dat$total_prod), 
        x=dat$i, col=c("green4", "purple", "yellow4", "black", "red"), log="", lty=c(1,1,1,1,2), lwd=c(1,1,1,2,1), type="l",
        ylab="Biomass pools", xlab="Year")
abline(h=0, col="grey")

matplot(y=cbind(dat$total_mass,
                dat$leaf_mass,
                dat$root_mass,
                dat$coarse_root_mass,
                dat$stem_mass),
        x=dat$i, col=c("black", "green3", "purple", "purple3", "brown"), log="", lty=c(1,1,1,1,1), lwd=c(1,1,1,1,1), type="l",
        ylab="Biomass pools", xlab="Year")
abline(h=0, col="grey")


matplot(y=cbind(dat$fitness), 
        x=dat$i, col=c("purple2", "magenta", "pink"), log="", lty=c(1,1,1), lwd=c(1,1,2), type="l",
        ylab="Fitness", xlab="Year")
abline(h=0, col="grey")

# plot(dat$total_mass~dat$i, ylab="Total biomass", xlab="Year")
# plot(dat$total_mass~dat$i, ylab="Total biomass", xlab="Year")
# 
# points(y=dat$total_prod-dat$litter_mass, x=dat$i, type="l", col="red")

# plot(y=1 - (dat$total_prod-dat$litter_mass)/dat$total_mass, x=dat$i, ylab="Total biomass", xlab="Year")

# plot(dat$ppfd[01:1000]~dat$i[01:1000], type="l")
# plot(dat$assim_gross[901:1000]~dat$ppfd[901:1000], type="l")


matplot(y=cbind(dat$assim_gross/dat$crown_area,
                dat$assim_net/dat$crown_area,
                dat$assim_net/dat$assim_gross * max(dat$assim_gross/dat$crown_area)), 
                x=dat$i, col=c("green3", "green4", scales::alpha("yellow3", 0.5)), log="", lty=1, type="l", lwd=c(1,1,3),
                ylab=expression(atop("GPP, NPP", "(kg m"^"-2"*"Yr"^"-1"*")")), xlab="Year")
abline(h=0, col="grey")

matplot(y=cbind(dat$rr/dat$crown_area,
                dat$rs/dat$crown_area,
                dat$rl/dat$crown_area), 
                x=dat$i, col=c("pink2", "pink3", "pink4"), log="", lty=1, type="l",
                ylab="Respiration", xlab="Year")

matplot(y=cbind(dat$tr/dat$crown_area,
                dat$tl/dat$crown_area), 
                x=dat$i, col=c("orange3", "orange4"), log="", lty=1, type="l",
                ylab="Turnover", xlab="Year",
                add=F)

plot(I(dat$transpiration/dat$crown_area/1000*1000)~dat$i, type="l", col="blue", ylab="Transpitation\n(mm/yr)")

plot(dat$dpsi~dat$i, type="l", col="cyan", ylab=expression(atop(Delta*psi, "(MPa)")), xlab="Year")

plot(dat$vcmax~dat$i, type="l", col="limegreen", ylab=expression(atop(V[cmax], "("*mu*"mol m"^"-2"*"s"^"-1"*")")), xlab="Year")

plot(I(dat$leaf_area/dat$crown_area)~dat$i, ylab="LAI", xlab="Year", type="l")

plot(I(dat$mortality)~dat$i, ylab="Cumulative\nMortality", xlab="Year", type="l")

plot(I(dat$mortality_inst[dat$diameter<0.5])~dat$diameter[dat$diameter<0.5], ylab="Instantaneous\nmortality rate", xlab="Diameter", type="l")

matplot(y=cbind(dat$sapwood_fraction, 
              dat$heartwood_fraction),
        x=dat$i, 
        ylab="Sap/Heart wood\nfraction", xlab="Year", type="l", col=c("yellowgreen", "brown4"), lty=1, lwd=1)

# plot(I(dat$ppfd)~dat$i, ylab="PPFD", xlab="Year", type="l", col="yellowgreen")

matplot(y=(cbind(dat$leaf_lifespan, 
                dat$fineroot_lifespan)),
        x=dat$height, 
        ylab="Leaf/fine root\nlifespan", xlab="Height", type="l", col=c("yellowgreen", "brown4"), lty=1, lwd=1)


# plot(y=2+(dat$rl+dat$rr+dat$rs)/dat$tl,x=dat$i)

# plot(dat$leaf_lifespan~dat$height)

# plot(y=(dat$height[2:2000]-dat$height[1:1999]), x=dat$i[1:1999], ylab="Height growth rate", type="l")

# plot(y=(dat$diameter[2:2000]-dat$diameter[1:1999])/dat$diameter[1:1999], x=dat$i[1:1999], ylab="Diameter RGR", type="l")

# detrend = function(x,y){
#   mod = lm(y~x)
#   y-fitted(mod)
# }
# plot(detrend(dat$i, dat$assim_gross)~dat$i, type="l")

# par(mfrow=c(2,1), mar=c(4,4,1,1))
# plot(scale((dat$assim_net/dat$crown_area))~dat$i, type="l")
# points(scale(dat$lai)~dat$i, type="l", col="red", lwd=3)
# 
# ccf(dat$lai, (dat$assim_net/dat$crown_area), lag.max=100)
# abline(v=0, col="red")
# abline(v=5, col="blue")

