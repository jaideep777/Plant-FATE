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
Gp   = read.delim("species_0_g.txt", header=F, col.names = paste0("V", 1:n))
Ap   = read.delim("species_0_gpp.txt", header=F, col.names = paste0("V", 1:n))

Zp   = read.delim("z_star.txt", header=F, col.names = paste0("V", 1:50))
COp   = read.delim("canopy_openness.txt", header=F, col.names = paste0("V", 1:50))
BAp   = read.delim("basal_area.txt", header=F, col.names = paste0("V", 1:50))
LAIp   = read.delim("LAI.txt", header=F, col.names = paste0("V", 1:2))

ids_x = 2:51
ids_h = which(diff(as.numeric(hp[1,ids_x])) > 0)[-1]
t = Dp[,1]

seeds = read.delim("seeds.txt", header=F, col.names = paste0("V", 1:50))

par(mfrow=c(5,3), mar=c(4,5,1,.5))

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


image(x=Dp[,1], y=as.numeric(Dp[1,ids_x]), z=log(1+1e8*pmax(as.matrix(Up[,ids_x])/max(as.matrix(Up[,ids_x]),na.rm=T),0)), log="y", xlab="Time", ylab = "Diameter distribution", col=scales::viridis_pal()(100))

image(x=hp[,1], y=as.numeric(hp[1,ids_h]), z=log(1+1e8*pmax(as.matrix(Up[,ids_h])/max(as.matrix(Up[,ids_h]),na.rm=T),0)), log="", xlab="Time", ylab = "Height distribution, Z*", col=scales::viridis_pal()(100))
matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*", add=T)

image(x=hp[,1], y=as.numeric(hp[1,ids_h]), z=log(1+1e4*pmax(as.matrix(Gp[,ids_h])/max(as.matrix(Gp[,ids_h]),na.rm=T),0)), log="", xlab="Time", ylab = "R. Growth rate dist, Z*", col=scales::viridis_pal()(100))
matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*", add=T)

image(x=hp[,1], y=as.numeric(hp[1,ids_h]), z=log(1+1e4*pmax(as.matrix(Ap[,ids_h])/max(as.matrix(Ap[,ids_h]),na.rm=T),0)), log="", xlab="Time", ylab = "GPP dist, Z*", col=scales::viridis_pal()(100))
matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*", add=T)

image(x=hp[,1], y=as.numeric(hp[1,ids_h]), z=as.matrix(Lp[,ids_h]), zlim=c(0,4), log="", xlab="Time", ylab = "LAI distribution, Z*", col=scales::viridis_pal()(100))
matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*", add=T)

# matplot(Dp$V1, Dp[,-1], lty=1, col=cols_t, type="l",
#         las=1, xlab="Time (years)", ylab="Diameter (m)")


image(x=Dp[,1], y=(as.numeric(Dp[1,ids_x])), z=log(as.matrix(Mp[,ids_x])), log="y", xlab="Time", ylab = "Mortality rate dist", col=scales::viridis_pal()(100))


matplot(Up$V1, (as.matrix(Up[,-1])*1e4/100), lty=1, col=cols_t, type="l",
        las=1, xlab="Time (years)", ylab="Density (Ind/cm/Ha)")

# matplot(Sp$V1, Sp[,-1]/1e6, lty=1, col=cols_t, type="l",
#         las=1, xlab="Time (years)", ylab="Seeds in pool (Millions)")

matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")

matplot(COp$V1, COp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Canopy openness")

matplot(BAp$V1, BAp[,-1]*1e4, lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Basal area (m2 / Ha)\n")

matplot(LAIp$V1, LAIp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Community LAI")

matplot(seeds$V1, seeds[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Species Seed output")


matplot(x=(t(hp[,-1])), y=(t(Up[,-1]))*1e4/100, lty=1, col=scales::alpha(rainbow(n = nrow((hp)), start = 0, end = 0.75),0.05), type="l",
        las=1, xlab = "Height", ylab="Density (Ind/cm/Ha)\n", log="xy", ylim=c(exp(-5),exp(10)))
# Up_mean = colMeans(Up[t>150, ids_h]*1e4/100)
# points(as.numeric(Up_mean)~as.numeric(hp[1,ids_h]))

matplot(x=(t(Dp[,-1])), y=(t(Up[,-1]))*1e4/100, lty=1, col=scales::alpha(rainbow(n = nrow((Dp)), start = 0, end = 0.75),0.05), type="l",
        las=1, xlab = "Diameter", ylab="Density (Ind/cm/Ha)\n", log="xy", ylim=c(exp(-15), exp(10)))
# Up_mean = colMeans(Up[t>150, ids_x]*1e4/100)
# points(as.numeric(Up_mean)~as.numeric(Dp[1,ids_x]))

matplot(y=t(as.matrix(Mp[c(1,nrow(Mp)),2:50])), x=as.numeric(Dp[1,2:50]), type="l", col=c("red", "black"), log="x", ylab="Mortality rate @(t0, tf)", xlab="Diameter")

# matplot(Lp$V1, Lp[,-1], lty=1, col=cols_t, type="l",
#         las=1, xlab="Time (years)", ylab="LAI")
# image(x=Dp[,1], y=as.numeric(Dp[1,ids_x]), z=as.matrix(Lp[,ids_x]), xlab="Time", ylab = "Diameter", col=scales::viridis_pal()(100))

# library(tidyverse)
# dat = list(rgr = seq(0.001,0.8, length.out = 100), D = seq(0.01,4, length.out = 100)) %>% cross_df() %>% 
#   mutate(logit_rgr = 1.51 + 0.64*log(1+100*rgr)) %>% 
#   mutate(mu_rgr = log(1 +exp(-logit_rgr)) ) %>% 
#   mutate(mu_d = exp(-5-0.3*log(D) + .3*(D)) ) %>% 
#   mutate(mu = mu_d + mu_rgr)
# 
# plot(dat$mu_rgr~dat$rgr)
# plot(dat$mu_d~dat$D)
# ii+# dat %>% ggplot(aes(y=mu, x=rgr, group=D, col=D))+ geom_line()

# matplot(x=(t(hp[nrow(hp),-1])), y=(t(Up[nrow(Up),-1])), log="xy", lty=1, col=scales::alpha(rainbow(n = nrow(hp), start = 0, end = 0.85),0.4), type="l",
#         las=1, xlab = "Height", ylab="Density")
# matplot(x=(t(Dp[nrow(hp),-1])), y=(t(Up[nrow(Up),-1])), log="xy", lty=1, col=scales::alpha(rainbow(n = nrow(hp), start = 0, end = 0.85),0.4), type="l",
#         las=1, xlab = "Diameter", ylab="Density")
