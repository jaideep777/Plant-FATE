dat1 = read.delim("D:/github_stuff/Plant-FATE_patch/muh.txt")
i = 0 
j = 1000
par(mfrow=c(3,3), mar=c(4,4,1,1), oma=c(1,1,1,1))
plot(dat1$swp[i:j], ylab="SWP", xlab="Year", type="l", lwd=2)
plot(dat1$h[i:j], ylab="hydraulics", xlab="Year", type="l", lwd=2)
plot(dat1$mu[i:j], ylab="Mortality", xlab="Year", col="brown", type="l", lwd=2)

dat = read.delim("D:/github_stuff/Plant-FATE_patch/assim.txt")
plot(dat$fecundity_mort[i:j], ylab="fec_mort", xlab="Year", type="l", lwd=2)
plot(dat$fecundity[i:j], ylab="fec", xlab="Year", type="l", lwd=2)
plot(exp(-dat$mortality)[i:j], ylab="mort_cumul", xlab="Year", type="l", lwd=2)
plot(dat$mortality_inst[i:j], ylab="mort_inst", xlab="Year", type="l", lwd=2)
plot(dat$germinated[i:j], ylab="germ", xlab="Year", type="l", lwd=2)
#plot(dat$mortality, ylab="height", xlab="Year", type="l", lwd=2)
#plot(diff((dat$fecundity_mort[50:150])))
