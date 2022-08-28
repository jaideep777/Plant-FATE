output_dir = "pspm_output_4"
prefix = "timestep_test"

solver = "IEBT0.1_succ2"
setwd(paste0("~/codes/Plant-FATE/",output_dir,"/",prefix,"_",solver))

plot_to_file = F
n_species = 5
n = 101

seeds = read.delim("seeds.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
Zp = read.delim("z_star.txt", header=F, col.names = paste0("V", 1:50))
BA = read.delim("basal_area.txt", header=F, col.names = paste0("V", 1:(n_species+2)))

par(mfcol=c(3,3), mar=c(6,6,1,1), oma=c(1,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(4,1,0))
matplot(seeds$V1, seeds[,-1], lty=1, col=rainbow(n = n_species+1, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Species Seed output", log="")
mtext(line=0.5, side=3, text=solver)
matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")

matplot(BA$V1, cbind(BA[,-1], rowSums(BA[,-1], na.rm=T))*1e4, lty=1, col=c(rainbow(n = n_species+1, start = 0, end = 0.85), "black"), type="l",
        las=1, xlab="Time (years)", ylab="Basal area", log="")
# 
# maxcol = 1200
# dp = read.delim("species_0_height.txt", header=F, col.names = paste0("V", 1:maxcol))
# up = read.delim("species_0_u.txt", header=F, col.names = paste0("V", 1:maxcol))
# mp = read.delim("species_0_mort.txt", header=F, col.names = paste0("V", 1:maxcol))
# rgr = read.delim("species_0_rgr.txt", header=F, col.names = paste0("V", 1:maxcol))
# lai = read.delim("species_0_lai.txt", header=F, col.names = paste0("V", 1:maxcol))
# 
# # id = as.integer(seq(1,min(nrow(up), nrow(dp)),length.out = 100))
# # matplot(x=t(dp[id,-1]), y=t(up[id,-1]), type="l", lty=1, col=rainbow(n = 101, start = 0, end = 0.85, alpha = 10/100), log="xy", ylim=c(1e-6,1e2))
# # points(x=colMeans(dp[-(1:1000),-1]), y=colMeans(up[-(1:1000),-1]))
# 
# nrows = min(nrow(dp), nrow(mp), nrow(rgr))-1
# ids = max((nrows-1000),1):nrows
# smoothScatter(as.numeric(as.matrix((mp[ids,-1])))~as.numeric(as.matrix((dp[ids,-1]))), xlim=c(0.01,.5))
# smoothScatter(as.numeric(as.matrix(mp[ids,-1]))~as.numeric(as.matrix(log(1000*rgr[ids,-1]*dp[ids,-1]))), xlim=c(0,0.2))
# smoothScatter(as.numeric(as.matrix((lai[ids,-1])))~as.numeric(as.matrix((dp[ids,-1]))), xlim=c(0.01, 2))

dat = read.delim("AmzFACE_D_PFATE_AMB_LD.txt")
plot(dat$LAI~dat$YEAR, type="l")
matplot(y=cbind(dat$GPP, dat$NPP)*1e-3*365, x=dat$YEAR, type="l", lty=1, col=c("green4", "green3"), ylab="GPP (kgC/m2/yr")
matplot(y=cbind(dat$GS), x=dat$YEAR, type="l", lty=1, col=c("cyan3"))
matplot(y=cbind(dat$ET), x=dat$YEAR, type="l", lty=1, col=c("blue"))
matplot(y=cbind(dat$CL+dat$CR+dat$CW)/1e3, x=dat$YEAR, type="l", lty=1, col=c("blue"))
matplot(y=cbind(dat$CFR)/1e3, x=dat$YEAR, type="l", lty=1, col=c("blue"))
