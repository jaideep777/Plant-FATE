library(tidyverse)

setwd("~/codes/libpspm/demo/")

### RED ####
# Analytical calculation of equilibrium distribution, if available
Ueq = function(x){
  a0<- 0.396
  phiA<- 0.749
  g0<-0.0838
  phiG<-0.7134
  m0<-1
  mort<-0.035
  mu0<-mort*m0/g0
  alpha<-0.10
  temp<- mu0/(1-phiG)
  coverage<-1-(1-alpha)/alpha*mu0/((mu0/(1-phiG))^(phiG/(phiG-1))*exp(mu0/(1-phiG))*as.numeric(pracma::gammainc(mu0/(1-phiG),phiG/(1-phiG)+1)[2]))
  Neq<-coverage/a0/(temp^(phiA/(phiG-1))*exp(temp)*as.numeric(pracma::gammainc(temp,phiA/(1-phiG)+1))[2])
  n0<-Neq*mort/g0
  return(n0*(x/m0)^(-phiG)*exp(mu0/(1-phiG)*(1-(x/m0)^(1-phiG)))*10000)
}

xeq_red = exp(seq(log(1), log(1e6), length.out=1000))
ueq_red = Ueq(xeq_red)

N_red = integrate(Ueq, 1, 1e6, abs.tol = 1e-6, rel.tol = 1e-6)
B_red = integrate(function(x){x*Ueq(x)}, 1, 1e6, abs.tol = 1e-8, rel.tol = 1e-8)

plot1 = function(file, N, title){
  dat = read.delim(file, header=F)
  dat = dat[,-ncol(dat)]
  x = exp(seq(log(1), log(1e6), length.out=N))
  
  cols = rainbow(nrow(dat)/1, start=0, end=0.9, alpha=10/nrow(dat))
  
  matplot(x=x, y = t(dat[seq(1,nrow(dat),by=1),-c(1,2)]), col=cols, 
          type="l", lty=1, log="xy", ylim=c(1e-20, 1e4), xlab="Size", ylab="Density", main=title)
  
  lines(x=xeq_red, y=ueq_red, col=scales::alpha("orange", 0.7), lwd=4)
  points(x=x, y = t(dat[nrow(dat),-c(1,2)]), col="black", 
         type="l", lty=1, lwd=1, cex = 0.7)
}


### DAPHNIA ####

a = 0.75
mu = 0.1
r = 0.5
K = 3

xstar = (mu*(1+mu)*(2+mu)/2/a)^(1/3)
sstar = xstar/(1-xstar)
xeq = seq(0,1,length.out = 1000)
u_equil = function(x){
  a*r*sstar*(1-sstar/K)*(xstar-x)^(mu-1)/xstar^mu
} 
ueq = u_equil(xeq)

N = integrate(u_equil, 0, xstar*0.99999999999, abs.tol = 1e-6, rel.tol = 1e-6)
B = integrate(function(x){x*u_equil(x)}, 0, xstar*0.99999999999, abs.tol = 1e-6, rel.tol = 1e-6)
Feq = integrate(function(x){(a*x^2*sstar/(1+sstar))*u_equil(x)}, 0, xstar*0.99999999999, abs.tol = 1e-6, rel.tol = 1e-6)

plot2 = function(file, N, title){
  dat = read.delim(file, header=F)
  dat = dat[,-ncol(dat)]
  x=seq(0,1,length.out = N)
  # plot(x=x, y=exp(-8*x^3), type="l")
  
  cols = rainbow(nrow(dat)/1, start=0, end=0.9, alpha=10/nrow(dat))
  # cols = scales::viridis_pal(direction = -1, alpha = 10/nrow(dat))(nrow(dat))
  matplot(x=x, y = t(dat[seq(1,nrow(dat),by=1),-c(1,2,3)]), col=cols, 
          type="l", lty=1, log="y", ylim=c(0.1, 1e4), xlab="Size", ylab="Density", main=title)
  abline(v=xstar, col="grey")
  
  lines(x=xeq, y=ueq, col=scales::alpha("orange", 0.7), lwd=4)
  points(x=x, y = t(dat[nrow(dat),-c(1,2,3)]), col="black", 
         type="l", lty=1, lwd=1, cex = 0.7)
  # plot(dat$V2~dat$V1, type="l")
}


# plot1("RED_model_nD/ebtn_Redmodel.txt", 150, "EBTn")


cairo_pdf("RED_Daphnia_explicit.pdf", width = 6/4*1.75*4, height=6)

par(mfrow = c(2,4), mar=c(4,4,4,1), oma=c(1,1,1,1), cex.lab=1.2, cex.axis=1.2)

plot1("RED_model/fmu_Redmodel.txt", 150, "FMU")
plot1("RED_model/ebt_Redmodel.txt", 150, "EBT")
plot1("RED_model/cm_Redmodel.txt", 150, "CM")
plot1("RED_model/abm_Redmodel.txt", 150, "ABM")

plot2("Daphnia_model/fmu_Daphnia.txt", 300, "FMU")
plot2("Daphnia_model/ebt_Daphnia.txt", 300, "EBT")
plot2("Daphnia_model/cm_Daphnia.txt", 300, "CM")
plot2("Daphnia_model/abm_Daphnia.txt", 300, "ABM")

dev.off()

cairo_pdf("RED_Daphnia_implicit.pdf", width = 6/4*1.75*3, height=6)

par(mfrow = c(2,3), mar=c(4,4,4,1), oma=c(1,1,1,1), cex.lab=1.2, cex.axis=1.2)


plot1("RED_model/ifmu_Redmodel.txt", 150, "IFMU")
# plot1("RED_model/ifmu2_Redmodel.txt", 150, "ILUD")
plot1("RED_model/iebt_Redmodel.txt", 150, "IEBT")
plot1("RED_model/icm_Redmodel.txt", 150, "ICM")

plot2("Daphnia_model/ifmu_Daphnia.txt", 300, "IFMU")
# plot2("Daphnia_model/ifmu2_Daphnia.txt", 300, "ILUD")
plot2("Daphnia_model/iebt_Daphnia.txt", 300, "IEBT")
plot2("Daphnia_model/icm_Daphnia.txt", 300, "ICM")

dev.off()


cols = scales::alpha(c("purple", "green3", "mediumspringgreen", "darkgoldenrod2", "red3", "pink", "#2b8cbe"), alpha=0.7)
cols = scales::alpha(c("darkgreen", 
                       "yellowgreen", 
                       # "green3", 
                       "magenta", 
                       "purple", 
                       "darkgoldenrod2", 
                       "darkgoldenrod3", 
                       "turquoise2"), alpha=0.7)
names = c("FMU", 
          "IFMU", 
          # "ILUD", 
          "EBT", 
          "IEBT", 
          "CM", 
          "ICM", 
          "ABM")

plotS = function(file, N, title){
  dat = read.delim(file, header=F)
  dat = dat[,-ncol(dat)]
  # plot(x=x, y=exp(-8*x^3), type="l")

  lines(x=dat$V1, y=dat$V3,col=cols[idx], type="l", lwd=2)
  idx <<- idx+1
}


par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
idx = 1
plot(x=1,y=NA, xlim=c(0,100), ylim=c(0,3))
plotS("Daphnia_model/fmu_Daphnia.txt", 300, "FMU")
plotS("Daphnia_model/ifmu_Daphnia.txt", 300, "IFMU")
# plotS("Daphnia_model/ifmu2_Daphnia.txt", 300, "ILUD")
plotS("Daphnia_model/ebt_Daphnia.txt", 300, "EBT")
plotS("Daphnia_model/iebt_Daphnia.txt", 300, "EBT")
plotS("Daphnia_model/cm_Daphnia.txt", 300, "CM")
plotS("Daphnia_model/icm_Daphnia.txt", 300, "ICM")
plotS("Daphnia_model/abm_Daphnia.txt", 300, "ABM")



plot_seeds = function(file){
  dat = read.delim(file, header=F)
  dat = dat[,-ncol(dat)]
  # plot(x=x, y=exp(-8*x^3), type="l")
  
  lines(x=dat$V1, y=dat$V2,col=cols[idx], type="l", lwd=2)
  idx <<- idx+1
}

cairo_pdf("RED_Daphnia_reproduction_rate_timeseries.pdf", width = 8.7, height=5)

par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
idx = 1
plot(x=1,y=NA, xlim=c(0,100), ylim=c(0,0.5), ylab="Reproduction rate", xlab="Time")
plot_seeds("Daphnia_model/fmu_Daphnia.txt")
mtext(side=3, line=0.5, text = "Daphnia model")
plot_seeds("Daphnia_model/ifmu_Daphnia.txt")
# plot_seeds("Daphnia_model/ifmu2_Daphnia.txt")
plot_seeds("Daphnia_model/ebt_Daphnia.txt")
plot_seeds("Daphnia_model/iebt_Daphnia.txt")
plot_seeds("Daphnia_model/cm_Daphnia.txt")
plot_seeds("Daphnia_model/icm_Daphnia.txt")
plot_seeds("Daphnia_model/abm_Daphnia.txt")
abline(h=Feq$value, col="black")

idx = 1
plot(x=1,y=NA, xlim=c(0,5000), ylim=c(0,50), ylab="Reproduction rate", xlab="Time")
mtext(side=3, line=0.5, text = "RED model")
plot_seeds("RED_model/fmu_Redmodel.txt")
plot_seeds("RED_model/ifmu_Redmodel.txt")
# plot_seeds("RED_model/ifmu2_Redmodel.txt")
plot_seeds("RED_model/ebt_Redmodel.txt")
plot_seeds("RED_model/iebt_Redmodel.txt")
plot_seeds("RED_model/cm_Redmodel.txt")
plot_seeds("RED_model/icm_Redmodel.txt")
plot_seeds("RED_model/abm_Redmodel.txt")
abline(h=37.5845, col="black")
legend(x = 2500, y=26, legend = names, col=cols, lwd=2, bty = "n", cex=1.)

dev.off()
