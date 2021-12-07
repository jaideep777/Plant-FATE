setwd("/home/jaideep/codes/tmodel_cpp/")

# dat_crown=read.delim("crown_area.txt", header=F)
# matplot(x= dat_crown[,1], y=dat_crown[,3:4], type="l", lty=1)
# abline(v = dat_crown[1,2], col="grey")
  
n = 201
Dp   = read.delim("species_0_X.txt", header=F, col.names = paste0("V", 1:n))
hp   = read.delim("species_0_height.txt", header=F, col.names = paste0("V", 1:n))
Up   = read.delim("species_0_u.txt", header=F, col.names = paste0("V", 1:n))
Lp   = read.delim("species_0_lai.txt", header=F, col.names = paste0("V", 1:n))
Mp   = read.delim("species_0_mort.txt", header=F, col.names = paste0("V", 1:n))
Sp   = read.delim("species_0_seeds.txt", header=F, col.names = paste0("V", 1:n))
Zp   = read.delim("z_star.txt", header=F, col.names = paste0("V", 1:50))
ids_x = 2:51
ids_h = which(diff(as.numeric(hp[1,ids_x])) > 0)[-1]
t = Dp[,1]

par(mfrow=c(3,3), mar=c(4,5,1,.5))

cols = rainbow(n = length(ids_x)+5, start = 0, end = 0.85)
cols_t = scales::alpha(cols, 0.5)

# matplot(hp$V1, hp[,-1], lty=1, col=cols, type="l",
#         las=1, xlab="Time (years)", ylab="Height (m)")

# plot(x=hp$V1, y=hp$V1*NA, ylim=c(0,20), xlab="Time (years)", ylab="Height (m)")
# for (i in 2:n){
#   alph = as.numeric(Up[i,-1]/max(Up[,-1], na.rm=T))
#   alph[is.na(alph)]=0
#   segments(hp$V1[i-1], as.numeric(hp[i-1,-1]), hp$V1[i], as.numeric(hp[i,-1]), col=rgb(t(col2rgb(cols)/255), alpha=alph) )
# }

  
image(x=Dp[,1], y=as.numeric(Dp[1,ids_x]), z=log(1e-20+as.matrix(Up[,ids_x])), log="y", xlab="Time", ylab = "Diameter distribution", col=scales::viridis_pal()(100))

image(x=hp[,1], y=as.numeric(hp[1,ids_h]), z=log(1e-20+as.matrix(Up[,ids_h])), log="", xlab="Time", ylab = "Height distribution, Z*", col=scales::viridis_pal()(100))
matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*", add=T)

matplot(Up$V1, (as.matrix(Up[,-1])), lty=1, col=cols_t, type="l",
        las=1, xlab="Time (years)", ylab="Density (m-2)")

# matplot(Dp$V1, Dp[,-1], lty=1, col=cols_t, type="l",
#         las=1, xlab="Time (years)", ylab="Diameter (m)")


matplot(Mp$V1, Mp[,-1], lty=1, col=cols_t, type="l",
        las=1, xlab="Time (years)", ylab="Mortality (yr-1)")

matplot(Sp$V1, Sp[,-1]/1e6, lty=1, col=cols_t, type="l",
        las=1, xlab="Time (years)", ylab="Seeds in pool (Millions)")

matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")

matplot(x=(t(hp[,-1])), y=(t(Up[,-1])), lty=1, col=scales::alpha(rainbow(n = nrow((hp)), start = 0, end = 0.75),0.05), type="l",
        las=1, xlab = "Height", ylab="Density", log="xy", ylim=c(exp(-20),exp(10)))
Up_mean = colMeans(Up[t>50, ids_x])
# points(as.numeric(Up_mean)~as.numeric(hp[1,ids_x]))

matplot(x=(t(Dp[,-1])), y=(t(Up[,-1])), lty=1, col=scales::alpha(rainbow(n = nrow((Dp)), start = 0, end = 0.75),0.05), type="l",
        las=1, xlab = "Diameter", ylab="Density", log="xy", ylim=c(exp(-20), exp(10)))
Up_mean = colMeans(Up[t>50, ids_x])
# points(as.numeric(Up_mean)~as.numeric(Dp[1,ids_x]))


matplot(Lp$V1, Lp[,-1], lty=1, col=cols_t, type="l",
        las=1, xlab="Time (years)", ylab="LAI")
# image(x=Dp[,1], y=as.numeric(Dp[1,ids_x]), z=as.matrix(Lp[,ids_x]), xlab="Time", ylab = "Diameter", col=scales::viridis_pal()(100))


# matplot(x=(t(hp[nrow(hp),-1])), y=(t(Up[nrow(Up),-1])), log="xy", lty=1, col=scales::alpha(rainbow(n = nrow(hp), start = 0, end = 0.85),0.4), type="l",
#         las=1, xlab = "Height", ylab="Density")
# matplot(x=(t(Dp[nrow(hp),-1])), y=(t(Up[nrow(Up),-1])), log="xy", lty=1, col=scales::alpha(rainbow(n = nrow(hp), start = 0, end = 0.85),0.4), type="l",
#         las=1, xlab = "Diameter", ylab="Density")
