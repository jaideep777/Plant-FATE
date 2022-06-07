dat = read.delim("D:/github_stuff/Plant-FATE_patch/muh.txt")

par(mfrow=c(4,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
plot(dat$swp, ylab="SWP", xlab="Year", type="l", lwd=2)
plot(dat$h, ylab="hydraulics", xlab="Year", type="l", lwd=2)
plot(dat$mu, ylab="Mortality", xlab="Year", col="brown", type="l", lwd=2)

dat1 = read.delim("D:/github_stuff/Plant-FATE_patch/mort.txt")
plot(dat1$X1.31275e.318, ylab="mort_cumul", xlab="Year", type="l", lwd=2)

dat = read.delim("D:/github_stuff/Plant-FATE_patch/assim.txt")
plot(dat$fecundity, type="l", col="black")
