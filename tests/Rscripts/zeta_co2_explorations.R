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


N=4
zeta = c(0.0667, 0.1, 0.2, 0.3)
output_dir = "pspm_output_amazon_nut2"

wd_mat_t = matrix(nrow = 2001, ncol = N)
wd_amb_ele_mat = matrix(nrow=2, ncol=N)
for(i in 1:4){
  expt = paste0("ELE_HD_zeta_", zeta[i])
  cat(expt, "\n")
  dat = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",expt,"/AmzFACE_Y_PFATE_ELE_HD.txt"))
  traits = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",expt,"/traits_ELE_HD.txt"))
  wd_amb_ele_mat[,i] = c(get_wd_cwm(dat, traits, 2000), get_wd_cwm(dat, traits, 3000))
  wd_mat_t[,i] = get_wd_cwm_t(dat, traits)$cwm
}


barplot(wd_amb_ele_mat, beside = T, ylim=c(550, 800), xpd=F, 
        legend.text = c("Ambient CO2 (368.9 ppm)", "Elevated CO2 (614 ppm)"), 
        axis.lty = 1, 
        main="Basal-area-weighted mean wood density",
        names.arg=paste0("zeta = ", zeta))

matplot(y=wd_mat_t, x=1000:3000, type="l", lty=1, lwd=2, col=scales::viridis_pal()(N))


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

par(mfrow=c(1,1), mar=c(4,6,1,1), cex.lab=1.3, cex.axis = 1.3, mgp=c(3,1,0), oma=c(1,1,1,1))

matplot(y=wd_mat, x=1000:1500, type="l", lty=1, lwd=2, col=scales::viridis_pal()(N), ylab="Wood density", xlab="Year")
# a=get_wd_cwm_t(dat, traits)$cwm
# matlines(y=a, x=1000:(1000+length(a)-1), type="l", lty=1, lwd=2, col=scales::viridis_pal()(N)[N])

par(mfrow=c(4,1), mar=c(4,6,1,1), cex.lab=1.3, cex.axis = 1.3, mgp=c(3,1,0), oma=c(1,1,1,1))

plot(y=npp_vec[1:N]*1e-3*365, x=zeta[1:N], type="o", xlab="Fine roor mass per unit leaf area (kg m-2)", ylab="NPP\n(kgC m-2 yr-1)")
abline(v=0.2, col="grey")

matplot(y=t(zstar_vec[,1:N]), x=zeta[1:N], lty=1, lwd=1,
        col=scales::viridis_pal(direction = -1, end = 0.9)(4), type="l",
        las=1, xlab="Fine roor mass per unit leaf area (kg m-2)", ylab="Canopy layer\nheights (m)")
abline(v=0.2, col="grey")

plot(wd_vec[1:N]~zeta[1:N], type="o", xlab="Fine roor mass per unit leaf area (kg m-2)", ylab="Wood density\n(kg m-3)")
abline(v=0.2, col="grey")

plot(y=agb_vec[1:N], x=zeta[1:N], type="o", xlab="Fine roor mass per unit leaf area (kg m-2)", ylab="AGB\n(kgC m-2)")
abline(v=0.2, col="grey")

bio_conv = agb_vec[1:N]/npp_vec[1:N]
points(y=bio_conv %>% scale_minmax(min(agb_vec[1:N]), max(agb_vec[1:N])), x=zeta[1:N], type="l", col="cyan3")
abline(v=0.2, col="grey")



