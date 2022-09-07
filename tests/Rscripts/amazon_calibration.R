library(tidyverse)
rm(list=ls())

output_dir = "pspm_output_4"
prefix = "lma_test"

solver = "IEBT0.1_succ3_nodist_wmutants_3spp_wd"
setwd(paste0("~/codes/Plant-FATE/",output_dir,"/",prefix,"_",solver))

plot_to_file = F
n_species = 9
n = 101

# seeds1 = read.delim("seeds.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
Zp = read.delim("z_star.txt", header=F, col.names = paste0("V", 1:50))
# BA1 = read.delim("basal_area.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
co = read.delim("canopy_openness.txt", header=F, col.names = paste0("V", 1:50))
lai_v = read.delim("lai_profile.txt", header=F, col.names = paste0("V", 1:27))
traits = read.delim("traits.txt")

dat = read.delim("AmzFACE_D_PFATE_AMB_LD.txt")
dat2 = read.delim("AmzFACE_Y_PFATE_AMB_LD.txt")


par(mfcol=c(3,4), mar=c(6,6,1,1), oma=c(1,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(4,1,0))
seeds = dat2 %>% select(YEAR, PID, SEEDS) %>% spread(value = "SEEDS", key = "PID")
matplot(seeds$YEAR, seeds[,-1], lty=1, col=rainbow(n = n_species+1, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Species Seed output", log="")
mtext(line=0.5, side=3, text=solver)

# matplot(seeds1$V1, seeds1[,-1], lty=1, col=rainbow(n = n_species+1, start = 0, end = 0.85), type="l",
#         las=1, xlab="Time (years)", ylab="Species Seed output", log="")
# mtext(line=0.5, side=3, text=solver)

matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")
matplot(co$V1, co[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Io")
matplot(y=1:25, x=t(lai_v[,2:26]), lty=1, col=rainbow(n = 1000, start = 0, end = 0.85, alpha=0.05), type="l",
        las=1, xlab="Cumulative LAI", ylab="Height")

BA = dat2 %>% select(YEAR, PID, BA) %>% spread(value = "BA", key = "PID")
matplot(BA$YEAR, cbind(BA[,-1], rowSums(BA[,-1], na.rm=T))*1e4, lty=1, col=c(rainbow(n = n_species, start = 0, end = 0.85), "black"), type="l",
        las=1, xlab="Time (years)", ylab="Basal area", log="")
abline(h=31.29, col="grey40")

# matplot(BA1$V1, cbind(BA1[,-1], rowSums(BA1[,-1], na.rm=T))*1e4, lty=1, col=c(rainbow(n = n_species+1, start = 0, end = 0.85), "black"), type="l",
#         las=1, xlab="Time (years)", ylab="Basal area", log="")
# abline(h=31.29, col="grey40")

# 
# maxcol = 1200
# dp = read.delim("species_0_height.txt", header=F, col.names = paste0("V", 1:maxcol))
# up = read.delim("species_0_u.txt", header=F, col.names = paste0("V", 1:maxcol))
# # mp = read.delim("species_0_mort.txt", header=F, col.names = paste0("V", 1:maxcol))
# rgr = read.delim("species_0_rgr.txt", header=F, col.names = paste0("V", 1:maxcol))
# lai = read.delim("species_0_lai.txt", header=F, col.names = paste0("V", 1:maxcol))
# 
# # id = as.integer(seq(1,min(nrow(up), nrow(dp)),length.out = 100))
# # matplot(x=t(dp[id,-1]), y=t(up[id,-1]), type="l", lty=1, col=rainbow(n = 101, start = 0, end = 0.85, alpha = 10/100), log="xy", ylim=c(1e-6,1e2))
# # points(x=colMeans(dp[-(1:1000),-1]), y=colMeans(up[-(1:1000),-1]))
# 


# nrows = min(nrow(dp), nrow(lai), nrow(rgr))-1
# ids = max((nrows-50),1):nrows
# smoothScatter(as.numeric(as.matrix((mp[ids,-1])))~as.numeric(as.matrix((dp[ids,-1]))), xlim=c(0.01,.5))
# smoothScatter(as.numeric(as.matrix(mp[ids,-1]))~as.numeric(as.matrix(log(1000*rgr[ids,-1]*dp[ids,-1]))), xlim=c(0,0.2))


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

cprops = read.delim("cohort_props.txt")

cprops %>% filter(speciesID==0) %>% filter(t > max(t)-50) %>% with(smoothScatter(lai~height))
cprops %>% filter(speciesID==0) %>% filter(t > max(t)-50) %>% with(smoothScatter(gpp~height))
cprops %>% filter(speciesID==0) %>% filter(t > max(t)-50) %>% with(smoothScatter(mort~height))

dists = read.delim("size_distributions.txt", header=F)
dists %>% filter(V2 == 0) %>% select(-V1, -V2, -V103) %>% as.matrix() %>% t() %>% matplot(x=exp(seq(log(0.01), log(10), length.out=100)), y=., type="l", lty=1, col=rainbow(n = ncol(.), start = 0, end = 0.85, alpha = 10/ncol(.)), log="xy", ylim=c(1e-6,1e2))
dists %>% filter(V2 == 1) %>% select(-V1, -V2, -V103) %>% as.matrix() %>% t() %>% matplot(x=exp(seq(log(0.01), log(10), length.out=100)), y=., type="l", lty=1, col=rainbow(n = ncol(.), start = 0, end = 0.85, alpha = 10/ncol(.)), log="xy", ylim=c(1e-6,1e2))
dists %>% filter(V2 == 2) %>% select(-V1, -V2, -V103) %>% as.matrix() %>% t() %>% matplot(x=exp(seq(log(0.01), log(10), length.out=100)), y=., type="l", lty=1, col=rainbow(n = ncol(.), start = 0, end = 0.85, alpha = 10/ncol(.)), log="xy", ylim=c(1e-6,1e2))

dists %>% filter(V2 == 2) %>% filter(V1 == 1500) %>% select(-V1, -V2, -V103) %>% as.numeric() %>% plot(x=exp(seq(log(0.01), log(10), length.out=100)), y=., type="l", lty=1, col="black", log="xy", ylim=c(1e-6,1e2))



# 
# par(mfcol=c(3,4), mar=c(6,6,1,1), oma=c(1,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(4,1,0))
# for (isp in 1:n_species){
#   traits_sp = traits %>% filter(SPP == isp) %>% filter(YEAR>1010)
#   traits_sp %>% with(matplot(y=log(cbind(r0_last, r0_avg, r0_exp)), x=YEAR, type="l", lty=1, col=rainbow(4, end = 0.75)))
# }
# 
# 
# par(mfcol=c(3,4), mar=c(6,6,1,1), oma=c(1,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(4,1,0))
# for (isp in 1:n_species){
#   traits_sp = traits %>% filter(SPP == isp) %>% filter(YEAR>1010)
#   traits_sp %>% with(matplot(y=WD, x=YEAR, type="l", lty=1, col=rainbow(4, end = 0.75)))
# }
  

traits %>% 
  filter(YEAR > 1120) %>% 
  # filter(resident==T) %>% 
  ggplot(aes(y=LMA, x=WD))+
  theme_classic(base_size = 18)+
  geom_point(aes(col=YEAR, size=RES), alpha=0.7)+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5))+
  ylim(c(0.115, 0.125))
#   xlim(c(450, 950))




# matplot(y=cbind(dat$ET), x=dat$YEAR, type="l", lty=1, col=c("blue"))
