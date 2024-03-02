Ueq = function(x){
  (1-x)^2/(1+x)^4
}

plot1 = function(file, title){
  x = seq(0,1,length.out=26)
  dat = read.delim(file, header=F)
  dat = dat[,-ncol(dat)]
  Y = dat
  matplot(x=x, y=t(Y[,-c(1,2)]), type = "l", lty=1, col=rainbow(nrow(Y)/1, start=0, end=0.9, alpha=10/nrow(Y)), xlab="Size", ylab="Density", ylim=c(0,2), main=title, log="")
  points(x=x, y=Ueq(x))
}

# png("~/codes/libpspm/testmodel.png", width=950*3, height = 740*3, res=300)
par(mfrow = c(3,3), mar=c(4,4,4,1), oma=c(1,1,1,1), cex.lab=1.2, cex.axis=1.2)
plot1("~/codes/libpspm/cm_testmodel_equil.txt", "CM")
plot1("~/codes/libpspm/ebt_testmodel_equil.txt", "EBT")
plot1("~/codes/libpspm/fmu_testmodel_equil.txt", "FMU")
plot1("~/codes/libpspm/ifmu_testmodel_equil.txt", "IFMU")
plot1("~/codes/libpspm/ifmu2_testmodel_equil.txt", "IFMU(O2)")
plot1("~/codes/libpspm/iebt_testmodel_equil.txt", "IEBT")
plot1("~/codes/libpspm/abm_testmodel_equil.txt", "ABM")
plot1("~/codes/libpspm/icm_testmodel_equil.txt", "ICM")



setwd("~/codes/libpspm")
dat2 = read.delim("fmu_testmodel_equil.txt", header=F)
dat3 = read.delim("ifmu_testmodel_equil.txt", header=F)
dat4 = read.delim("ifmu2_testmodel_equil.txt", header=F)
dat5 = read.delim("ebt_testmodel_equil.txt", header=F)
dat6 = read.delim("iebt_testmodel_equil.txt", header=F)
dat7 = read.delim("cm_testmodel_equil.txt", header=F)
dat8 = read.delim("abm_testmodel_equil.txt", header=F)

cols_m = c("purple", "green3", "mediumspringgreen", "darkgoldenrod2", "red3", "pink", "#2b8cbe")

par(mfrow=c(1,1), mar=c(4,4,3,1), oma=c(1,1,1,1))
matplot(y=cbind(dat2$V2, dat3$V2, dat4$V2, dat5$V2, dat6$V2, dat7$V2, dat8$V2),
              x=dat2$V1, type="l", ylab = "u(0)", xlab = "Time", main="Long-term performance",
            col=scales::alpha(cols_m, 0.75), lty=1, lwd=c(2,2,2,2,2,0.75))
# dev.off()


