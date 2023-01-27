library(tidyverse)
rm(list=ls())

ref_dir = "pspm_output_cont_test3/cont_test_ref"
spinup_dir = "pspm_output_cont_test3/cont_test_spinup"
main_dir   = "pspm_output_cont_test3/cont_test_main"

setwd(paste0("~/codes/Plant-FATE/",spinup_dir))
dat2_spin = read.delim("AmzFACE_Y_PFATE_ELE_HD.txt")
dat_spin = read.delim("AmzFACE_D_PFATE_ELE_HD.txt")
Yend = max(dat2_spin$YEAR)

setwd(paste0("~/codes/Plant-FATE/",main_dir))
dat2_main = read.delim("AmzFACE_Y_PFATE_ELE_HD.txt")
dat_main = read.delim("AmzFACE_D_PFATE_ELE_HD.txt")

setwd(paste0("~/codes/Plant-FATE/",ref_dir))
dat2_ref = read.delim("AmzFACE_Y_PFATE_ELE_HD.txt")
dat_ref = read.delim("AmzFACE_D_PFATE_ELE_HD.txt")


dat2_main = dat2_main %>% filter(YEAR > Yend)
dat_main = dat_main %>% filter(YEAR > Yend)

dat2 = rbind(dat2_spin, dat2_main)
dat = rbind(dat_spin, dat_main)

plot_data = function(dat, dat2, main){
  n_species = length(unique(dat2$PID))
  n_year = length(unique(dat2$YEAR))
  
  par(mfrow=c(1,3), mar=c(6,6,1,1), oma=c(1,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(3.2,1,0), las=1)
  
  seeds = dat2 %>% select(YEAR, PID, SEEDS) %>% spread(value = "SEEDS", key = "PID")
  matplot(seeds$YEAR, seeds[,-1], lty=1, col=rainbow(n = n_species+1, start = 0, end = 0.85, alpha = 0.5), type="l",
          las=1, xlab="Time (years)", ylab="Species seed\noutput", log="", main=main)
  abline(v=1200, col="pink")
  
  BA = dat2 %>% select(YEAR, PID, BA) %>% spread(value = "BA", key = "PID")
  matplot(BA$YEAR, cbind(BA[,-1], rowSums(BA[,-1], na.rm=T))*1e4, lty=1, col=c(rainbow(n = n_species, start = 0, end = 0.85), "black"), type="l",
          las=1, xlab="Time (years)", ylab="Basal area", log="")
  abline(v=1200, col="pink")
  
  plot(dat$LAI~dat$YEAR, type="l", col="red3", ylim=c(0,max(dat$LAI,6.5)), xlab="Time (years)", ylab="Total LAI")
  abline(v=1200, col="pink")
}

plot_data(dat_ref, dat2_ref, "Ref")
plot_data(dat, dat2, "Main")
