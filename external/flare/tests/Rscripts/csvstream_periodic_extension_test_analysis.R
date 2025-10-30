dat = read.delim("~/codes/Flare.v2/csvstream_met.txt", header=F)
colnames(dat) = c("julian", "t", "Year", "Month", "Year.dec", "Temp", "VPD", "PAR", "ppfd_max", "SWP", "X")
dato = read.csv("~/codes/Flare.v2/tests/data/gf-guy_drivers_plantfate.csv", header=T)
# datco2 = read.csv("~/codes/tmodel_cpp/tests/data/CO2_AMB_AmzFACE2000_2100.csv", header=T)
dato$t = dato$decimal_year # dato$year + (dato$month-1)/12
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

par(mfrow=c(5,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
plot(dat$Temp~dat$t, type="l", col="black")
points(dato$temp~dato$t, pch=20, type="l", col="red", cex=0.6)
plot(dat$VPD~dat$t, type="l", col="black")
points(dato$vpd~dato$t, pch=20, type="l", col="red", cex=0.6)
plot(dat$PAR~dat$t, type="l", col="black")
points(dato$par~dato$t, pch=20, type="l", col="red", cex=0.6)
plot(dat$SWP~dat$t, type="l", col="black")
points(dato$swp~dato$t, pch=20, type="l", col="red", cex=0.6)
with(dat[as.integer(dat$t) == 2008,], plot(Temp~t, type="l", col="black"))
with(dato[dato$year == 2008,], points(temp~t, pch=20, type="l", col="red", cex=0.6))


