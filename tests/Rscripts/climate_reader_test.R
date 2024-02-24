dat = read.delim("~/codes/Plant-FATE/climate_inst.txt", header=F)
colnames(dat) = c("t", "Temp", "VPD", "PAR", "SWP", "co2")
data = read.delim("~/codes/Plant-FATE/climate_acclim.txt", header=F)
colnames(data) = c("t", "Temp", "VPD", "PAR", "SWP", "co2")

dato = read.csv("~/codes/Plant-FATE/tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE_new.csv", header=T)
# datco2 = read.csv("~/codes/tmodel_cpp/tests/data/CO2_AMB_AmzFACE2000_2100.csv", header=T)
dato$t = dato$Year + (dato$Month-1)/12
dato$VPD = dato$VPD*100
dato$SWP = -dato$SWP
# dato_rep = rbind(dato, dato, dato, dato, dato, dato, dato, dato)
# dato_rep$Year = c(rep(1996,12), rep(1997,12),
#               rep(1998,12), rep(1999,12),
#               rep(2000,12), rep(2001,12),
#               rep(2002,12), rep(2003,12),
#               rep(2004,12), rep(2005,12),
#               rep(2006,12), rep(2007,12),
#               rep(2008,12), rep(2009,12),
#               rep(2010,12), rep(2011,12))
# dato_rep$t = dato_rep$Year + (dato_rep$Month-1)/12

par(mfrow=c(5,1), mar=c(2,4,1,1), oma=c(1,1,1,1), xlab="")
plot(dat$Temp~dat$t, type="l", col="black", xlim=c(2000,2005))
points(dato$Temp~dato$t, pch=20, type="p", col="red", cex=0.6)
plot(dat$Temp~dat$t, type="l", col="black")
points(dato$Temp~dato$t, pch=20, type="l", col="red", cex=0.6)
plot(dat$VPD~dat$t, type="l", col="black")
points(dato$VPD~dato$t, pch=20, type="l", col="red", cex=0.6)
plot(dat$PAR~dat$t, type="l", col="black")
points(dato$PAR~dato$t, pch=20, type="l", col="red", cex=0.6)
plot(dat$SWP~dat$t, type="l", col="black")
points(dato$SWP~dato$t, pch=20, type="l", col="red", cex=0.6)


par(mfrow=c(5,1), mar=c(2,4,1,1), oma=c(1,1,1,1), xlab="")
matplot(y=cbind(dat$Temp,data$Temp), x=dat$t, type="l", lty=1, col=c("black", "green3"), xlim=c(2000,2005))
matplot(y=cbind(dat$VPD,data$VPD), x=dat$t, type="l", lty=1, col=c("black", "green3"), xlim=c(2000,2005))
matplot(y=cbind(dat$co2,data$co2), x=dat$t, type="l", lty=1, col=c("black", "green3"), xlim=c(2000,2005), ylim=c(360,400))
matplot(y=cbind(dat$PAR,data$PAR), x=dat$t, type="l", lty=1, col=c("black", "green3"), xlim=c(2000,2005))
matplot(y=cbind(dat$SWP,data$SWP), x=dat$t, type="l", lty=1, col=c("black", "green3"), xlim=c(2000,2005))








# 
# dat = read.delim("~/codes/tmodel_cpp/climate.txt", header=F)
# dato = read.csv("~/codes/tmodel_cpp/tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv", header=T)
# datco2 = read.csv("~/codes/tmodel_cpp/tests/data/CO2_AMB_AmzFACE2000_2100.csv", header=T)
# dato$t = dato$Year + (dato$Month-1)/12
# dato$VPD = dato$VPD*100
# dato$SWP = -dato$SWP
# # dato = rbind(dato, dato, dato, dato, dato, dato, dato, dato)
# # dato$Year = c(rep(1996,12), rep(1997,12), 
# #               rep(1998,12), rep(1999,12), 
# #               rep(2000,12), rep(2001,12), 
# #               rep(2002,12), rep(2003,12),
# #               rep(2004,12), rep(2005,12),
# #               rep(2006,12), rep(2007,12),
# #               rep(2008,12), rep(2009,12),
# #               rep(2010,12), rep(2011,12))
# dato$t = dato$Year + (dato$Month-1)/12
# 
# par(mfrow=c(5,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
# plot(dat$V2~dat$V1, type="l", col="black")
# points(dato$Temp~dato$t, pch=20, type="p", col="red")
# plot(dat$V3~dat$V1, type="l", col="black")
# points(dato$VPD~dato$t, pch=20, type="p", col="red")
# plot(dat$V4~dat$V1, type="l", col="black")
# points(dato$PAR~dato$t, pch=20, type="p", col="red")
# plot(dat$V5~dat$V1, type="l", col="black")
# points(dato$SWP~dato$t, pch=20, type="p", col="red")
# plot(dat$V6~dat$V1, type="l", col="black")
# points(datco2$CO2~datco2$Year, pch=20, type="p", col="red")
