library(tidyverse)

get_wd_cwm = function(dat2, traits, yr){
  dat2 %>% select(YEAR, PID, BA) %>% 
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
    filter(YEAR == yr) %>% summarize(cwm= sum(BA*WD)/sum(BA)) %>% as.numeric()
}

get_wd_cwm_t = function(dat2, traits){
  dat2 %>% select(YEAR, PID, BA) %>% 
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
    group_by(YEAR) %>% summarize(cwm= sum(BA*WD)/sum(BA))
}

scale_minmax = function(x, xmin, xmax){
  xmin + (x-min(x))*(xmax-xmin)/(max(x)-min(x))
}


N=2
zeta = c(0.08, 0.2)
output_dir = "pspm_output_calib_amb2"

wd_mat_t = matrix(nrow = 2001, ncol = N)
wd_amb_ele_mat = matrix(nrow=2, ncol=N)
ba_amb_ele_mat = matrix(nrow=2, ncol=N)
npp_amb_ele_mat = matrix(nrow=2, ncol=N)
agb_amb_ele_mat = matrix(nrow=2, ncol=N)
cfr_amb_ele_mat = matrix(nrow=2, ncol=N)
for(i in 1:N){
  expt = sprintf("HIST_ELE7_zeta_%.6f", zeta[i])
  cat(expt, "\n")
  dat = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",expt,"/AmzFACE_Y_PFATE_ELE_HD.txt"))
  dat1 = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",expt,"/AmzFACE_D_PFATE_ELE_HD.txt"))
  traits = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",expt,"/traits_ELE_HD.txt"))
  wd_amb_ele_mat[,i] = c(get_wd_cwm(dat, traits, 2000), get_wd_cwm(dat, traits, 3000))
  npp_amb_ele_mat[1,i] = dat1 %>% filter(YEAR > 1950 & YEAR < 2000) %>% select(NPP) %>% colMeans()
  npp_amb_ele_mat[2,i] = dat1 %>% filter(YEAR > 2950 & YEAR < 3000) %>% select(NPP) %>% colMeans()
  agb_amb_ele_mat[1,i] = dat1 %>% filter(YEAR > 1950 & YEAR < 2000) %>% mutate(AGB=CL+CW) %>% select(AGB) %>% colMeans()
  agb_amb_ele_mat[2,i] = dat1 %>% filter(YEAR > 2950 & YEAR < 3000) %>% mutate(AGB=CL+CW) %>% select(AGB) %>% colMeans()
  cfr_amb_ele_mat[1,i] = dat1 %>% filter(YEAR > 1950 & YEAR < 2000) %>% select(CFR) %>% colMeans()
  cfr_amb_ele_mat[2,i] = dat1 %>% filter(YEAR > 2950 & YEAR < 3000) %>% select(CFR) %>% colMeans()
  wd_mat_t[,i] = get_wd_cwm_t(dat, traits)$cwm
}


matplot(y=wd_mat_t, x=1000:3000, type="l", lty=1, lwd=2, col=scales::viridis_pal(end = 0.8)(N), ylab="Wood density", xlab="Year")

cairo_pdf(filename = "~/codes/Plant-FATE/paper_figs2/zeta_x_co2_barplots.pdf", height = 19/3, width = 8/3)
par(mfrow=c(4,1), mar=c(4,4,.5,1), oma=c(1,1,1,1))

barplot(wd_amb_ele_mat*NA, beside = T, ylim=c(0, 1), xlim=c(0,1), xpd=F, 
        legend.text = c("Ambient CO2 (368.9 ppm)", "Elevated CO2 (614 ppm)"), 
        axis.lty = 0, axes = F,
        #col=scales::viridis_pal(end = 0.8)(N),
        ylab="",
        args.legend = list(x=0.9,y=0.35, bty = "n")
)

barplot(wd_amb_ele_mat, beside = T, ylim=c(700, 850), xpd=F, 
#        legend.text = c("Ambient CO2 (368.9 ppm)", "Elevated CO2 (614 ppm)"), 
        axis.lty = 1, 
        #col=scales::viridis_pal(end = 0.8)(N),
        ylab="Wood density",
        names.arg=paste0("zeta = ", zeta))

barplot(npp_amb_ele_mat*1e-3*365, beside = T, ylim=c(0,3), xpd=F, 
        #legend.text = c("Ambient CO2 (368.9 ppm)", "Elevated CO2 (614 ppm)"), 
        axis.lty = 1, 
        ylab="NPP",
        names.arg=paste0("zeta = ", zeta))

barplot(agb_amb_ele_mat*1e-3, beside = T, ylim=c(20,45), xpd=F, 
        #legend.text = c("Ambient CO2 (368.9 ppm)", "Elevated CO2 (614 ppm)"), 
        axis.lty = 1, 
        ylab="AGB",
        names.arg=paste0("zeta = ", zeta))

dev.off()



#### CO2 scan ####

N=7
wd_vec = numeric(N)
npp_vec = numeric(N)
agb_vec = numeric(N)
wd_mat = matrix(nrow = 501, ncol = N)
zstar_vec = matrix(nrow = 20, ncol = N)
co2 = seq(360, 600, length.out=12)
for (i in 1:N){
  output_dir = "pspm_output_co2_scan"
  expt = sprintf("scan_co2_%0.6f", co2[i])
  cat(expt, "\n")
  dat1 = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",expt,"/AmzFACE_D_PFATE_ELE_HD.txt"))
  dat2 = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",expt,"/AmzFACE_Y_PFATE_ELE_HD.txt"))
  zstar = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",expt,"/z_star.txt"), header=F, col.names = paste0("V", 1:21))
  traits = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",expt,"/traits_ELE_HD.txt"))
  wd_vec[i] = get_wd_cwm(dat2, traits, 1500)
  npp_vec[i] = dat1 %>% filter(YEAR > 1450) %>% select(NPP) %>% colMeans()
  agb_vec[i] = dat1 %>% filter(YEAR == 1500) %>% mutate(AGB=(CL+CW)/1e3) %>% select(AGB) %>% as.numeric()
  wd_mat[,i] = get_wd_cwm_t(dat2, traits)$cwm
  zstar_vec[,i] = zstar %>% filter(V1 > 1450) %>% select(-V1) %>% colMeans() %>% as.numeric()
}

# par(mfrow=c(1,1), mar=c(4,6,1,1), cex.lab=1.3, cex.axis = 1.3, mgp=c(3,1,0), oma=c(1,1,1,1))
# 
# matplot(y=wd_mat, x=1000:1500, type="l", lty=1, lwd=2, col=scales::viridis_pal()(N), ylab="Wood density", xlab="Year")
# a=get_wd_cwm_t(dat, traits)$cwm
# matlines(y=a, x=1000:(1000+length(a)-1), type="l", lty=1, lwd=2, col=scales::viridis_pal()(N)[N])

cairo_pdf(filename = "paper_figs/co2_zeta_scan.pdf", height = 19/3, width = 18/3)
par(mfcol=c(4,2), mar=c(2,6,1,1), cex.lab=1.3, cex.axis = 1.3, mgp=c(3,1,0), oma=c(5,1,1,1) )

matplot(y=t(zstar_vec[,1:N]), x=co2[1:N], lty=1, lwd=1,
        col=scales::viridis_pal(direction = -1, end = 0.9)(4), type="l",
        las=1, xlab="", ylab="Canopy layer\nheights (m)")
abline(v=400, col="grey")

plot(wd_vec[1:N]~co2[1:N], type="o", xlab="", ylab="Wood density\n(kg m-3)")
abline(v=400, col="grey")

plot(y=npp_vec[1:N]*1e-3*365, x=co2[1:N], type="o", xlab="", ylab="NPP\n(kgC m-2 yr-1)")
abline(v=400, col="grey")

plot(y=agb_vec[1:N], x=co2[1:N], type="o", xlab="", ylab="AGB\n(kgC m-2)")
abline(v=400, col="grey")
mtext(side=1, line=3, text = "CO2 (ppm)")

# bio_conv = agb_vec[1:N]/npp_vec[1:N]
# points(y=bio_conv %>% scale_minmax(min(agb_vec[1:N]), max(agb_vec[1:N])), x=co2[1:N], type="l", col="cyan3")
# abline(v=400, col="grey")

# dev.off()

#### Zeta scan

N=11
wd_vec = numeric(N)
npp_vec = numeric(N)
agb_vec = numeric(N)
wd_mat = matrix(nrow = 501, ncol = N)
zstar_vec = matrix(nrow = 20, ncol = N)
zeta = seq(0.10-3*0.2/9, 0.30, length.out=13)
for (i in 1:N){
  output_dir = "pspm_output_zeta_scan"
  expt = sprintf("HIST_zeta_%0.6f", zeta[i])
  cat(expt, "\n")
  dat1 = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",expt,"/AmzFACE_D_PFATE_ELE_HD.txt"))
  dat2 = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",expt,"/AmzFACE_Y_PFATE_ELE_HD.txt"))
  zstar = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",expt,"/z_star.txt"), header=F, col.names = paste0("V", 1:21))
  traits = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",expt,"/traits_ELE_HD.txt"))
  wd_vec[i] = get_wd_cwm(dat2, traits, 1500)
  npp_vec[i] = dat1 %>% filter(YEAR > 1450) %>% select(NPP) %>% colMeans()
  agb_vec[i] = dat1 %>% filter(YEAR == 1500) %>% mutate(AGB=(CL+CW)/1e3) %>% select(AGB) %>% as.numeric()
  wd_mat[,i] = get_wd_cwm_t(dat2, traits)$cwm
  zstar_vec[,i] = zstar %>% filter(V1 > 1450) %>% select(-V1) %>% colMeans() %>% as.numeric()
}

# par(mfrow=c(1,1), mar=c(4,6,1,1), cex.lab=1.3, cex.axis = 1.3, mgp=c(3,1,0), oma=c(1,1,1,1))
# 
# matplot(y=wd_mat, x=1000:1500, type="l", lty=1, lwd=2, col=scales::viridis_pal()(N), ylab="Wood density", xlab="Year")
# # a=get_wd_cwm_t(dat, traits)$cwm
# # matlines(y=a, x=1000:(1000+length(a)-1), type="l", lty=1, lwd=2, col=scales::viridis_pal()(N)[N])
# 
# cairo_pdf(filename = "paper_figs/zeta_scan.pdf", height = 19/3, width = 10/3)
# par(mfrow=c(4,1), mar=c(4,6,1,1), cex.lab=1.3, cex.axis = 1.3, mgp=c(3,1,0), oma=c(1,1,1,1))

matplot(y=t(zstar_vec[,1:N]), x=zeta[1:N], lty=1, lwd=1,
        col=scales::viridis_pal(direction = -1, end = 0.9)(4), type="l",
        las=1, xlab="", ylab="Canopy layer\nheights (m)")
abline(v=0.2, col="grey")

plot(wd_vec[1:N]~zeta[1:N], type="o", xlab="", ylab="Wood density\n(kg m-3)")
abline(v=0.2, col="grey")

plot(y=npp_vec[1:N]*1e-3*365, x=zeta[1:N], type="o", xlab="", ylab="NPP\n(kgC m-2 yr-1)")
abline(v=0.2, col="grey")

plot(y=agb_vec[1:N], x=zeta[1:N], type="o", xlab="", ylab="AGB\n(kgC m-2)")
abline(v=0.2, col="grey")
mtext(side=1, line=4, text="Fine root mass per unit\nleaf area (kg m-2)")
# bio_conv = agb_vec[1:N]/npp_vec[1:N]
# points(y=bio_conv %>% scale_minmax(min(agb_vec[1:N]), max(agb_vec[1:N])), x=zeta[1:N], type="l", col="cyan3")
# abline(v=0.2, col="grey")

dev.off()



#### Optimization ####

library(PlantFATE)

# Ambient Env
zp_amb = c(19.22669, 15.71661, 3.15710, 0.00000)
co_amb = c(1.0000000, 0.4712540, 0.2218492, 0.0711068)

# eCO2 env
zp_ele = c(20.695389, 18.106550, 14.087510, 2.206985, 0.000000) 
co_ele = c(1.000000, 4.712543e-01, 2.220795e-01, 1.032055e-01, 2.763073e-02)

setwd("~/codes/Plant-FATE/")

fitness_wd = function(wd, zp, co, co2, zeta=0.2, lma, p50, pfile="tests/params/p.ini"){
  lho = new(LifeHistoryOptimizer, pfile)
  lho$set_metFile("tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv")
  lho$traits0$lma = lma
  lho$traits0$p50_xylem = p50 
  lho$traits0$zeta = zeta
  lho$env$clim$co2 = co2
  lho$env$z_star = zp
  lho$env$canopy_openness = co
  lho$init()
  lho$set_traits(c(0.12, wd))
  # lho$printPlant()
  lho$calcFitness()
}

wd = seq(350, 900, length.out=20)

dat_amb = matrix(nrow=50, ncol=length(wd))
for (i in 1:50){
  cat(i, "\n")
  fitness = wd %>% purrr::map_dbl(fitness_wd, zp=zp_amb, co=co_amb, co2=414, lma=runif(1,0.05,0.25), p50 = runif(1,-4.5,-0.5))
  dat_amb[i,] = fitness
}
matplot(y=t(dat_amb), x=wd, lty=1, type="l", col=scales::alpha("black", 0.2))
dat_amb_opt = wd[apply(dat_amb, MARGIN=1, FUN = function(x){which(x==max(x))})]

dat_ele_base = wd %>% purrr::map_dbl(fitness_wd, zp=zp_amb, co=co_amb, co2=614)
dat_ele = wd %>% purrr::map_dbl(fitness_wd, zp=zp_ele, co=co_ele, co2=614)
# dat_amb_nfert = wd %>% purrr::map_dbl(fitness_wd, zp=zp_amb, co=co_amb, co2=368, pfile="tests/params/p_zeta_0.1.ini")

wd_amb = wd[which(dat_amb==max(dat_amb))]
wd_ele_base = wd[which(dat_ele_base==max(dat_ele_base))]
wd_ele = wd[which(dat_ele==max(dat_ele))]
# wd_amb_nfert = wd[which(dat_amb_nfert==max(dat_amb_nfert))]


cairo_pdf(filename = "paper_figs/optimal_wd.pdf", height = 7, width = 5.5)
par(mfrow=c(2,1), mar=c(3,6,1,1), oma=c(4,1,1,1))
matplot(y=cbind(dat_amb/max(dat_amb), dat_ele_base/max(dat_ele_base), dat_ele/max(dat_ele)), 
        x = wd, type="l", lty=1, col=c("black", "grey", "yellow3", "brown"), 
        ylab="Fitness", xlab="", cex.lab=1.3, lwd=2,
        xlim=c(350, 900))

barplot(c(wd_amb, wd_ele_base, wd_ele)[c(3,2,1)], xlim=c(330,950), xpd=F, 
        names.arg=c("aC+aE", "eC+aE", "eC+eE")[c(3,2,1)],
        col=c("black", "grey", "yellow3", "brown")[c(3,2,1)], horiz = T, las=1)

mtext(text="Wood density (kg m-3)", side=1, line=3)
dev.off()

