output_dir = "pspm_output_4"
prefix = "lma_test"

solver = "IEBT0.1_succ3_nodist"
setwd(paste0("~/codes/Plant-FATE/",output_dir,"/",prefix,"_",solver))

plot_to_file = F
n_species = 3
n = 101

seeds = read.delim("seeds.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
Zp = read.delim("z_star.txt", header=F, col.names = paste0("V", 1:50))
BA = read.delim("basal_area.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
co = read.delim("canopy_openness.txt", header=F, col.names = paste0("V", 1:50))
lai_v = read.delim("lai_profile.txt", header=F, col.names = paste0("V", 1:27))

par(mfcol=c(3,4), mar=c(6,6,1,1), oma=c(1,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(4,1,0))
matplot(seeds$V1, seeds[,-1], lty=1, col=rainbow(n = n_species+1, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Species Seed output", log="")
mtext(line=0.5, side=3, text=solver)

matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")
matplot(co$V1, co[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Io")
matplot(y=1:25, x=t(lai_v[,2:26]), lty=1, col=rainbow(n = 1000, start = 0, end = 0.85, alpha=0.05), type="l",
        las=1, xlab="Cumulative LAI", ylab="Height")

matplot(BA$V1, cbind(BA[,-1], rowSums(BA[,-1], na.rm=T))*1e4, lty=1, col=c(rainbow(n = n_species+1, start = 0, end = 0.85), "black"), type="l",
        las=1, xlab="Time (years)", ylab="Basal area", log="")
abline(h=31.29, col="grey40")
# 
maxcol = 1200
dp = read.delim("species_0_height.txt", header=F, col.names = paste0("V", 1:maxcol))
up = read.delim("species_0_u.txt", header=F, col.names = paste0("V", 1:maxcol))
# mp = read.delim("species_0_mort.txt", header=F, col.names = paste0("V", 1:maxcol))
rgr = read.delim("species_0_rgr.txt", header=F, col.names = paste0("V", 1:maxcol))
lai = read.delim("species_0_lai.txt", header=F, col.names = paste0("V", 1:maxcol))
# 
# # id = as.integer(seq(1,min(nrow(up), nrow(dp)),length.out = 100))
# # matplot(x=t(dp[id,-1]), y=t(up[id,-1]), type="l", lty=1, col=rainbow(n = 101, start = 0, end = 0.85, alpha = 10/100), log="xy", ylim=c(1e-6,1e2))
# # points(x=colMeans(dp[-(1:1000),-1]), y=colMeans(up[-(1:1000),-1]))
# 
nrows = min(nrow(dp), nrow(lai), nrow(rgr))-1
ids = max((nrows-50),1):nrows
# smoothScatter(as.numeric(as.matrix((mp[ids,-1])))~as.numeric(as.matrix((dp[ids,-1]))), xlim=c(0.01,.5))
# smoothScatter(as.numeric(as.matrix(mp[ids,-1]))~as.numeric(as.matrix(log(1000*rgr[ids,-1]*dp[ids,-1]))), xlim=c(0,0.2))

dat = read.delim("AmzFACE_D_PFATE_AMB_LD.txt")

plot(dat$LAI~dat$YEAR, type="l", col="red", ylim=c(0,max(dat$LAI,6.5)))
abline(h=c(5.3, 6.2), col=scales::muted("red"))
abline(h=c(3.5), col=scales::muted("grey100"))

matplot(y=cbind(dat$GPP, dat$NPP)*1e-3*365, x=dat$YEAR, type="l", lty=1, col=c("green4", "green3"), ylab="GPP (kgC/m2/yr")
abline(h=c(3,3.5), col=scales::muted("green4"))
abline(h=c(1.31), col=scales::muted("green3"))

matplot(y=cbind(dat$GS), x=dat$YEAR, type="l", lty=1, col=c("cyan3"))
abline(h=c(0.16), col=scales::muted("cyan3"))

matplot(y=cbind(dat$CL+dat$CW)/1e3, x=dat$YEAR, type="l", lty=1, col=c("yellow3"), ylab="AGB (kg/m2)")
abline(h=c(16.9, 20.7), col=scales::muted("yellow3"))

matplot(y=cbind(dat$CFR)/1e3, x=dat$YEAR, type="l", lty=1, col=c("brown"), ylab="C-FR (kg/m2)", ylim=c(0, max(dat$CFR/1e3,0.7)))
abline(h=c(0.48, 0.66), col=scales::muted("brown"))

matplot(y=cbind(dat$VCMAX), x=dat$YEAR, type="l", lty=1, col=c("green3"), ylab="Vcmax (umol/m2/s)")
abline(h=c(40), col=scales::muted("green3"))

smoothScatter(as.numeric(as.matrix((lai[ids,-1])))~as.numeric(as.matrix((dp[ids,-1]))))


# matplot(y=cbind(dat$ET), x=dat$YEAR, type="l", lty=1, col=c("blue"))
