# R script to test:
dat = read.delim("~/codes/tmodel_cpp/assim.txt")
dat$leaf_area = dat$crown_area * dat$lai
dat$heartwood_fraction = 1-dat$sapwood_fraction

# par(mfrow=c(4,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
# 
# plot(dat$height~dat$diameter)
# plot(dat$crown_area~I(dat$height*dat$diameter))
# plot(dat$crown_area~dat$height)
# plot(dat$crown_area~dat$diameter)
# plot(dat$leaf_area~dat$crown_area)
# plot(dat$sapwood_fraction~dat$height)
# # plot(sqrt(4*dat$crown_area/pi)~dat$height)
# plot(I(dat$total_rep/dat$total_prod)~dat$height)
# plot(I(dat$total_rep/dat$total_prod)~dat$diameter)

par(mfrow=c(4,3), mar=c(4,4,1,1), oma=c(1,1,1,1))

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

matplot(y=cbind(dat$seed_pool,
                dat$germinated), 
        x=dat$i, col=c("purple2", "magenta"), log="", lty=c(1,1), lwd=c(1,1), type="l",
        ylab="Seed pools", xlab="Year")
abline(h=0, col="grey")

# plot(dat$total_mass~dat$i, ylab="Total biomass", xlab="Year")
# plot(dat$total_mass~dat$i, ylab="Total biomass", xlab="Year")
# 
# points(y=dat$total_prod-dat$litter_mass, x=dat$i, type="l", col="red")

# plot(y=1 - (dat$total_prod-dat$litter_mass)/dat$total_mass, x=dat$i, ylab="Total biomass", xlab="Year")

# plot(dat$ppfd[01:1000]~dat$i[01:1000], type="l")
# plot(dat$assim_gross[901:1000]~dat$ppfd[901:1000], type="l")


matplot(y=cbind(dat$assim_gross/dat$crown_area,
                dat$assim_net/dat$crown_area), 
                x=dat$i, col=c("green3", "green4"), log="", lty=1, type="l",
                ylab="GPP, NPP", xlab="Year")
abline(h=0, col="grey")

matplot(y=cbind(dat$rr,
                dat$rs,
                dat$rl), 
                x=dat$i, col=c("pink2", "pink3", "pink4"), log="", lty=1, type="l",
                ylab="Respiration", xlab="Year")

matplot(y=cbind(dat$tr,
                dat$tl), 
                x=dat$i, col=c("orange3", "orange4"), log="", lty=1, type="l",
                ylab="Turnover", xlab="Year",
                add=F)

plot(I(dat$transpiration/dat$crown_area/1000*1000)~dat$i, type="l", col="blue", ylab="Transpitation (mm/yr)")
plot(dat$dpsi~dat$i, type="l", col="cyan")
plot(dat$vcmax~dat$i, type="l", col="limegreen")

plot(I(dat$leaf_area/dat$crown_area)~dat$i, ylab="LAI", xlab="Year", type="l")
# plot(y=(dat$height[2:2000]-dat$height[1:1999]), x=dat$i[1:1999], ylab="Height growth rate", type="l")
plot(y=(dat$diameter[2:2000]-dat$diameter[1:1999])/dat$diameter[1:1999], x=dat$i[1:1999], ylab="Diameter RGR", type="l")

# detrend = function(x,y){
#   mod = lm(y~x)
#   y-fitted(mod)
# } 
# plot(detrend(dat$i, dat$assim_gross)~dat$i, type="l") 

# par(mfrow=c(2,1), mar=c(4,4,1,1))
# plot(scale((dat$assim_net/dat$crown_area)[100:200])~dat$i[100:200], type="l")
# points(scale(dat$lai[100:200])~dat$i[100:200], type="l", col="red", lwd=3)
# 
# ccf(dat$lai[201:2000], (dat$assim_net/dat$crown_area)[201:2000], lag.max=100)
# abline(v=0, col="red")
# abline(v=5, col="blue")

