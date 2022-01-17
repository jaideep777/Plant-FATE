dat = read.delim("~/codes/tmodel_cpp/climate.txt", header=F)
dato = read.csv("~/codes/tmodel_cpp/tests/data/MetData_AmzFACE_Monthly_2000_2001_PlantFATE.csv", header=T)
datco2 = read.csv("~/codes/tmodel_cpp/tests/data/CO2_AMB_AmzFACE2000_2100.csv", header=T)
dato$t = dato$Year + (dato$Month-1)/12
dato$VPD = dato$VPD*100
dato = rbind(dato, dato, dato, dato, dato, dato)
dato$Year = c(rep(2000,12), rep(2001,12), 
              rep(2002,12), rep(2003,12),
              rep(2004,12), rep(2005,12),
              rep(2006,12), rep(2007,12),
              rep(2008,12), rep(2009,12),
              rep(2010,12), rep(2011,12))
dato$t = dato$Year + (dato$Month-1)/12

par(mfrow=c(5,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
plot(dat$V2~dat$V1, type="l", col="black")
points(dato$Temp~dato$t, pch=20, type="p", col="red")
plot(dat$V3~dat$V1, type="l", col="black")
points(dato$VPD~dato$t, pch=20, type="p", col="red")
plot(dat$V4~dat$V1, type="l", col="black")
points(dato$PAR~dato$t, pch=20, type="p", col="red")
plot(dat$V5~dat$V1, type="l", col="black")
points(dato$SWP~dato$t, pch=20, type="p", col="red")
plot(dat$V6~dat$V1, type="l", col="black")
points(datco2$CO2~datco2$Year, pch=20, type="p", col="red")
