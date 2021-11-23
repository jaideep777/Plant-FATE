# R script to test:
dat = read.delim("~/codes/tmodel_cpp/geometric_growth_2.txt")
dat$heartwood_fraction = 1-dat$sapwood_fraction
dat$hv = (pi*0.25*dat$diameter^2*dat$sapwood_fraction)/dat$leaf_area

par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
plot(dat$height~dat$diameter)
plot(dat$crown_area~I(dat$height*dat$diameter))
plot(dat$crown_area~dat$height)
plot(dat$crown_area~dat$diameter)
plot(dat$leaf_area~dat$crown_area)
plot(dat$sapwood_fraction~dat$height)
plot(dat$sapwood_turnover_rate~dat$height)
# plot(sqrt(4*dat$crown_area/pi)~dat$height)
  
par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
plot(dat$sapwood_fraction_ode~dat$height)
points(dat$sapwood_fraction~dat$height, type="l", col="red")
plot(dat$total_mass~dat$height, log="")
points(dat$sapwood_mass~dat$height, col="brown")
points(dat$heartwood_mass~dat$height, col="grey")
points(dat$heartwood_mass_ode~dat$height, col="blue", type="l")
plot(I(dat$sapwood_mass/dat$leaf_area)~dat$height)
plot(dat$hv~dat$height)


par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
plot(dat$height~dat$i)
plot(dat$diameter~dat$i)
plot(dat$total_mass~dat$i)
points(dat$total_prod~dat$i, type="l", col="red")



