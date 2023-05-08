get_wd_cwm = function(dat2){
  c(dat2 %>% select(YEAR, PID, BA) %>%
      left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
      filter(YEAR == 2000) %>% summarize(cwm= sum(BA*WD)/sum(BA)) %>% as.numeric(),
    dat2 %>% select(YEAR, PID, BA) %>%
      left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
      filter(YEAR == 3000) %>% summarize(cwm= sum(BA*WD)/sum(BA)) %>% as.numeric()
  )
}

output_dir = "pspm_output_amazon_nut3"
prefix = "ELE"

solver = "HD_hi_zeta_0.25"#_old_params"
dat2_hi = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",prefix,"_",solver,"/AmzFACE_Y_PFATE_ELE_HD.txt"))
wd_hi = get_wd_cwm(dat2_hi)

solver = "HD_low_zeta_0.15"#_old_params"
dat2_low = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",prefix,"_",solver,"/AmzFACE_Y_PFATE_ELE_HD.txt"))
wd_low = get_wd_cwm(dat2_low)

solver = "HD_base_zeta_0.2"#_old_params"
dat2_base = read.delim(paste0("~/codes/Plant-FATE/",output_dir,"/",prefix,"_",solver,"/AmzFACE_Y_PFATE_ELE_HD.txt"))
wd_base = get_wd_cwm(dat2_base)

# dev.off()
par(mar=c(5,5,5,1), mgp=c(3,2,0))
barplot(cbind(wd_low, wd_base, wd_hi), beside = T, ylim=c(550, 800), xpd=F, 
        legend.text = c("Ambient CO2 (368.9 ppm)", "Elevated CO2 (614 ppm)"), 
        axis.lty=1, names.arg = c("High N\nzeta=0.0667", "Baseline\nzeta=0.2", "Low N\nzeta=0.3"),
        main="Basal-area-weighted mean wood density")
