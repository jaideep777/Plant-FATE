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
  
par(mfrow=c(3,2), mar=c(4,5,3,.5))

matplot(hp$V1, hp[,-1], lty=1, col=scales::alpha("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", 0))

matplot(Dp$V1, Dp[,-1], lty=1, col=scales::alpha("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Diameter (m)", main = paste0("With p1 | seed rain = ", 0))

matplot(Lp$V1, Lp[,-1], lty=1, col=scales::alpha("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="LAI", main = paste0("With p1 | seed rain = ", 0))

matplot(Up$V1, (as.matrix(Up[,-1])), lty=1, col=scales::alpha("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Density (m-2)", main = paste0("With p1 | seed rain = ", 0))

matplot(Mp$V1, Mp[,-1], lty=1, col=scales::alpha("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Mortality (yr-1)", main = paste0("With p1 | seed rain = ", 0))

matplot(Sp$V1, Sp[,-1]/1e6, lty=1, col=scales::alpha("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Seeds in pool (Millions)", main = paste0("With p1 | seed rain = ", 0))
