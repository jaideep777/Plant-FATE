# R script to test:
dat = read.delim("~/codes/tmodel_cpp/assim.txt")
dat$heartwood_fraction = 1-dat$sapwood_fraction

par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))

plot(dat$height~dat$diameter)
plot(dat$crown_area~I(dat$height*dat$diameter))
plot(dat$crown_area~dat$height)
plot(dat$crown_area~dat$diameter)
plot(dat$leaf_area~dat$crown_area)
plot(dat$sapwood_fraction~dat$height)
# plot(sqrt(4*dat$crown_area/pi)~dat$height)
  
plot(dat$height~dat$i)
plot(dat$diameter~dat$i)
plot(dat$total_mass~dat$i)
points(dat$total_prod~dat$i, type="l", col="red")

# plot(dat$ppfd[01:1000]~dat$i[01:1000], type="l")
# plot(dat$assim_gross[901:1000]~dat$ppfd[901:1000], type="l")


matplot(y=cbind(dat$assim_gross,
                dat$assim_net), 
                x=dat$i, col=c("green3", "green4"), log="", lty=1, type="l")

matplot(y=cbind(dat$rr,
                dat$rs,
                dat$rl), 
                x=dat$i, col=c("pink2", "pink3", "pink4"), log="", lty=1, type="l")

matplot(y=cbind(dat$tr,
                dat$tl), 
                x=dat$i, col=c("orange3", "orange4"), log="", pch=20)
