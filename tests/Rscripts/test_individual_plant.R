# R script to test:
dat = read.delim("~/codes/tmodel_cpp/assim.txt")
dat$leaf_area = dat$crown_area * dat$lai
dat$heartwood_fraction = 1-dat$sapwood_fraction

par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))

plot(dat$height~dat$diameter)
plot(dat$crown_area~I(dat$height*dat$diameter))
plot(dat$crown_area~dat$height)
plot(dat$crown_area~dat$diameter)
plot(dat$leaf_area~dat$crown_area)
plot(dat$sapwood_fraction~dat$height)
# plot(sqrt(4*dat$crown_area/pi)~dat$height)

par(mfrow=c(3,3), mar=c(4,4,1,1), oma=c(1,1,1,1))

plot(dat$height~dat$i, ylab="Height", xlab="Year")
plot(dat$diameter~dat$i, ylab="Diameter", xlab="Year", col="brown")
plot(dat$total_mass~dat$i, ylab="Total biomass", xlab="Year")
points(dat$total_prod~dat$i, type="l", col="red")

plot(I(dat$leaf_area/dat$crown_area)~dat$i, ylab="LAI", xlab="Year", type="l")

# plot(dat$ppfd[01:1000]~dat$i[01:1000], type="l")
# plot(dat$assim_gross[901:1000]~dat$ppfd[901:1000], type="l")


matplot(y=cbind(dat$assim_gross,
                dat$assim_net), 
                x=dat$i, col=c("green3", "green4"), log="", lty=1, type="l",
                ylab="GPP, NPP", xlab="Year")

matplot(y=cbind(dat$rr,
                dat$rs,
                dat$rl), 
                x=dat$i, col=c("pink2", "pink3", "pink4"), log="", lty=1, type="l",
                ylab="Respiration, Turnover", xlab="Year")

matplot(y=cbind(dat$tr,
                dat$tl), 
                x=dat$i, col=c("orange3", "orange4"), log="", pch=20,
                ylab="Turnover", xlab="Year",
                add=T)


plot(I(dat$transpiration/dat$crown_area/1000*1000)~dat$i, type="l", col="blue", ylab="Transpitation (mm/yr)")
plot(dat$dpsi~dat$i, type="l", col="cyan")
plot(dat$vcmax~dat$i, type="l", col="limegreen")

# detrend = function(x,y){
#   mod = lm(y~x)
#   y-fitted(mod)
# }
# plot(detrend(dat$i, dat$assim_gross)~dat$i, type="l")
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot(scale(dat$ppfd)~dat$i, type="l")
points(scale(dat$lai)~dat$i, type="l", col="red", lwd=3)

plot(scale(dat$ppfd[1:100])~dat$i[1:100], type="l")
points(scale(dat$lai[1:100])~dat$i[1:100], type="l", col="red", lwd=3)

       