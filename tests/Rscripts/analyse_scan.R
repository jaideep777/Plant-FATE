library(tidyverse)
rm(list=ls())

df = list(eg = c(-2, -1.5, -1.0),
     ea = c(-2, -1.5, -1.0),
     ewd = c(-2, -1.5, -1.0, -0.5)) %>% cross_df()

df[,c("GPP", "NPP", "RAU", "GS", "LAI", "VCMAX", "BA", "CL", "AGB", "BGB")] = NA
df$WD = NA

# k = 1
for (k in 1:nrow(df)){

output_dir = "/home/jjoshi/output_data/jaideep/pspm_output_36sims"
expt_dir = sprintf("par_ea%g_eg%g_ewd%g_614ppm", df$ea[k], df$eg[k], df$ewd[k]) #_old_params"

setwd(paste0(output_dir,"/",expt_dir))

plot_to_file = T
plot_trait_space = F

add_band = function(){
  # polygon(x=c(2000,5000,5000,2000), y=c(-1e20,-1e20,1e20,1e20), border = NA, col=scales::alpha("yellow2",0.2))
}

add_hband = function(ylim, col=scales::alpha("grey30",0.2), xlim=c(-1e20,5000)){
  polygon(y=c(ylim[1],ylim[2],ylim[2],ylim[1]), x=c(xlim[1],xlim[1],xlim[2],xlim[2]), border = NA, col=col)
}
# 
# seeds1 = read.delim("seeds.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
Zp = read.delim("z_star.txt", header=F, col.names = paste0("V", 1:50))
# BA1 = read.delim("basal_area.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
co = read.delim("canopy_openness.txt", header=F, col.names = paste0("V", 1:50))
lai_v = read.delim("lai_profile.txt", header=F, col.names = paste0("V", 1:27))
traits = read.delim("traits_ELE_HD.txt")
dat = read.delim("AmzFACE_D_PFATE_ELE_HD.txt")
dat2 = read.delim("AmzFACE_Y_PFATE_ELE_HD.txt")
dist = read.delim("size_distributions.txt", header=F)
dist = dist[,-ncol(dist)]
x = exp(seq(log(0.01), log(10), length.out=100))
traits_obs = read.csv(file = "~/codes/Plant-FATE/tests/data/Amz_trait_filled_HD.csv")
traits_used = read.csv(file = "~/codes/Plant-FATE/tests/data/Traits_random_HD2.csv")

# To get avg size distribution, sum over species and average over years
names(dist)[1:2] = c("YEAR", "SPP")
dist_amb = dist %>% filter(YEAR>1100 & YEAR<2000) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)
dist_ele = dist %>% filter(YEAR>4100 & YEAR<5000) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)

# dist_amb = dist %>% filter(YEAR == 1100) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)
# dist_ele = dist %>% filter(YEAR == 1101) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)

n_species = length(unique(dat2$PID))
n_year = length(unique(dat2$YEAR))

if (plot_to_file) png(paste0(output_dir, "emg_props_", expt_dir,".png"), width=2412*1.5*4/5, height = 1472*1.5, res=300)

par(mfcol=c(4,4), mar=c(4.5,6,.5,1), oma=c(1,1,2,1), cex.lab=1.1, cex.axis=1.1, mgp=c(3.2,1,0), las=1)
seeds = dat2 %>% select(YEAR, PID, SEEDS) %>% spread(value = "SEEDS", key = "PID")
matplot(seeds$YEAR, seeds[,-1], lty=1, col=rainbow(n = n_species+1, start = 0, end = 0.85, alpha = min(10/n_species, 1)), type="l",
        las=1, xlab="Time (years)", ylab="Species seed\noutput", log="")
mtext(line=0.5, side=3, text=expt_dir)
add_band()

# matplot(seeds1$V1, seeds1[,-1], lty=1, col=rainbow(n = n_species+1, start = 0, end = 0.85), type="l",
#         las=1, xlab="Time (years)", ylab="Species Seed output", log="")
# mtext(line=0.5, side=3, text=solver)

BA = dat2 %>% select(YEAR, PID, BA) %>% spread(value = "BA", key = "PID")
matplot(BA$YEAR, cbind(BA[,-1], rowSums(BA[,-1,drop=FALSE], na.rm=T))*1e4, lty=1, col=c(rainbow(n = n_species, start = 0, end = 0.85), "black"), type="l",
        las=1, xlab="Time (years)", ylab="Basal area", log="")
add_hband(c(31.29, 31.29*1.02))
add_band()

matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")
add_band()
# matplot(co$V1, co[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
#         las=1, xlab="Time (years)", ylab="Io")
matplot(y=1:24, x=t(-lai_v[,3:26]+lai_v[,2:25]), lty=1, col=rainbow(n = n_year, start = 0, end = 0.85, alpha=0.05), type="l",
        las=1, xlab="Leaf area density", ylab="Height")
matplot(y=1:25, x=t(lai_v[,2:26]), lty=1, col=rainbow(n = n_year, start = 0, end = 0.85, alpha=0.05), type="l",
        las=1, xlab="Cumulative LAI", ylab="Height")


plot(dat$LAI~dat$YEAR, type="l", col="red3", ylim=c(0,max(dat$LAI,6.5)), xlab="Time (years)", ylab="Total LAI")
# abline(h=c(5.3, 6.2), col=scales::muted("red"))
add_hband(c(5.3, 6.2))#, col=scales::alpha("red3", 0.2))
add_band()
# abline(h=c(3.5), col=scales::muted("grey100"))


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



matplot(y=cbind(dat$GPP, dat$NPP)*1e-3*365, x=dat$YEAR, type="l", lty=1, col=c("green4", "green3"), ylab="GPP, NPP\n(kgC/m2/yr)", xlab="Time (years)")
# points(y=dat$NPP/dat$GPP*4, x=dat$YEAR, type="l", lty=1, col=c("yellow1"))
# abline(h=c(3,3.5), col="grey")
add_hband(c(3,3.5))#, col=scales::alpha("black",0.3))
add_hband(c(1.31,1.3555))#, col=scales::alpha("black",0.3))
# abline(h=c(1.31), col=scales::muted("green3"))
add_band()

matplot(y=cbind(dat$GS), x=dat$YEAR, type="l", lty=1, col=c("cyan3"), ylab="Stomatal conductance\n(mol/m2/s)", xlab="Time (years)")
add_hband(c(0.16, 0.16555))#, col=scales::alpha("cyan4", 0.6))
# abline(h=c(0.16), col=scales::muted("cyan3"))
add_band()

agb = cbind(dat$CL+dat$CW)/1e3
matplot(y=agb, x=dat$YEAR, type="l", lty=1, col=c("yellow4"), ylim=c(0,max(agb)), ylab="AGB\n(kgC/m2)", xlab = "Time (years)")
add_hband(c(16.9, 20.7))#, col=scales::alpha("yellow3", 0.3))
add_band()

matplot(y=cbind(dat$CFR)/1e3, x=dat$YEAR, type="l", lty=1, col=c("brown"), ylab="C-FR\n(kgC/m2)", xlab = "Time (years)", ylim=c(0, max(dat$CFR/1e3,0.7)))
add_hband(c(0.48, 0.66))#, col=scales::alpha("brown", 0.3))
add_band()

matplot(y=cbind(dat$VCMAX), x=dat$YEAR, type="l", lty=1, col=c("green3"), ylab="Vcmax\n(umol/m2/s)", xlab="Time (years)")
add_hband(c(20,45)) #, col=scales::muted("green4"))
add_band()

plot_size_dist = function(){
  # matplot(y=cbind(as.numeric(colMeans(filter(dist, V1>1100 & V1<2000)[, -c(1,2)], na.rm = T)),
  #                 as.numeric(colMeans(filter(dist, V1>2100 & V1<3000)[, -c(1,2)], na.rm = T))
  #                 )*1e-2*1e4, 
  #         x=x, type="l", log="y", lty=1, col=c("black", "brown"), 
  #         xlim=c(0.01, 2), ylim=c(1e-4, 200), ylab="Density (stems/cm/ha)", xlab="Diameter (m)")
  matplot(y=cbind(as.numeric(dist_amb[gtools::mixedsort(names(dist_amb))][-1]),
                  as.numeric(dist_ele[gtools::mixedsort(names(dist_ele))][-1])
  )*1e-2*1e4, # Convert stems m-1 m-2 --> stems cm-1 ha-1
  x=x, type="l", log="y", lty=1, col=c("black", "yellow3"), 
  xlim=c(0.01, 1.2), ylim=c(1e-4, 1000), ylab="Density\n(stems/cm/ha)", xlab="Diameter (m)", las=0)
  
  abline(v=1, col=scales::alpha("red", 0.2))
  
  xobs = c(15,25,35,45,55,65,75,85,95,105)/100
  # Data for Manaus from https://link.springer.com/article/10.1007/s00442-004-1598-z 
  yobs=c(350.5221340921042,
         132.41927918860426,
         62.62503296462008,
         29.61724892214378,
         15.095996574802413,
         5.702923697662178,
         2.3219542502889836,
         1.5968055466971947,
         0.7006940913385968,
         0.5597156879584093)/10
  points(yobs~xobs, pch=20, col=scales::alpha("grey30", 0.4), cex=1.7)
}
try(plot_size_dist())

traits_obs %>% select(Leaf.LMA..g.m2., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =Leaf.LMA..g.m2., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.02), las=0, main="", xlab="LMA", col=NA, lwd=2)
traits_obs %>% select(Leaf.LMA..g.m2., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =Leaf.LMA..g.m2., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
  with(density(x =LMA*1000, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
try(
  dat2 %>% select(YEAR, PID, BA) %>%
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(YEAR == min(3000, max(dat2$YEAR)-1)) %>%
    with(density(x =LMA*1000, weights=BA/sum(BA))) %>% points(col="yellow3", type="l", lwd=1.5)
)

traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.005), las=0, main="", xlab="Wood density", col=NA, lwd=2)
traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
  with(density(x =WD, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
try(
  dat2 %>% select(YEAR, PID, BA) %>%
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(YEAR == min(3000, max(dat2$YEAR)-1)) %>%
    with(density(x =WD, weights=BA/sum(BA))) %>% points(col="yellow3", type="l", lwd=1.5)
)

traits_obs %>% select(Height_Max.m., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =Height_Max.m., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.18), las=0, main="", xlab="Max. height", col=NA, lwd=2)
traits_obs %>% select(Height_Max.m., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =Height_Max.m., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
try(
  dat2 %>% select(YEAR, PID, BA) %>% 
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
    filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
    with(density(x =HMAT, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
)
try(
  dat2 %>% select(YEAR, PID, BA) %>%
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(YEAR == min(3000, max(dat2$YEAR)-1)) %>%
    with(density(x =HMAT, weights=BA/sum(BA))) %>% points(col="yellow3", type="l", lwd=1.5)
)

# traits_obs %>% select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =P50..Mpa., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.5), las=0, main="", xlab="Max. height", col=NA, lwd=2)
# traits_obs %>% select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =P50..Mpa., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
# dat2 %>% select(YEAR, PID, BA) %>% 
#   filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
#   left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>% 
#   drop_na %>% 
#   with(density(x =P50..Mpa., weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
# try(
#   dat2 %>% select(YEAR, PID, BA) %>% 
#     filter(YEAR == 3000) %>% 
#     left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>% 
#     drop_na %>% 
#     with(density(x =P50..Mpa., weights=BA/sum(BA))) %>% points(col="yellow3", type="l", lwd=1.5)
# )


# 
# traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.005), las=0, main="", xlab="Wood density", col=NA, lwd=2)
# traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
# dat2 %>% select(YEAR, PID, BA) %>% 
#   filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
#   left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>% 
#   left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
#   # drop_na %>% 
#   with(density(x =meanWoodDensity..g.cm3.*1000, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
# try(
#   dat2 %>% select(YEAR, PID, BA) %>% 
#     filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
#     left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>% 
#     left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
#     #drop_na %>% 
#     with(density(x =WD, weights=BA/sum(BA))) %>% points(col="yellow3", type="l", lwd=1.5)
# )
# 
# 

cwm_wd = traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% summarise(cwm=sum(meanWoodDensity..g.cm3.*Total.BasalArea_2017.cm2.)/sum(Total.BasalArea_2017.cm2.))*1000
cwm_wd_pred = dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  drop_na %>% group_by(YEAR) %>% 
  summarize(cwm_wd = sum(WD*BA)/sum(BA))

cwm_wd_pred %>% with(plot(cwm_wd~YEAR, type="l")) #, ylim=c(200,900)))
add_hband(c(cwm_wd,cwm_wd+5))
add_band()

# cwm_p50 = traits_obs %>% select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% drop_na %>% summarise(cwm=sum(P50..Mpa.*Total.BasalArea_2017.cm2.)/sum(Total.BasalArea_2017.cm2.))
# cwm_p50_pred = dat2 %>% select(YEAR, PID, BA) %>% 
#   left_join(traits_obs, by = c("PID"="Species")) %>% 
#   drop_na %>% group_by(YEAR) %>% 
#   summarize(cwm_p50 = sum(P50..Mpa.*BA)/sum(BA))
# 
# cwm_p50_pred %>% with(plot(cwm_p50~YEAR, type="l")) #, ylim=c(200,900)))
# add_hband(c(cwm_wd,cwm_wd+5))

if (plot_to_file) dev.off()

mat = cbind(dat$GPP*1e-3*365, 
            dat$NPP*1e-3*365, 
            dat$RAU*1e-3*365,
            dat$GS,
            dat$LAI,
            dat$VCMAX,
            rowSums(BA[,-1], na.rm=T)*1e4,
            (dat$CL)*1e-3,
            (dat$CW+dat$CL)*1e-3,
            (dat$CCR+dat$CFR)*1e-3) %>% as.data.frame()
colnames(mat) = c("GPP", "NPP", "RAU", "GS", "LAI", "VCMAX", "BA", "CL", "AGB", "BGB")

df[k,colnames(mat)] = as.list(colMeans(mat))
df$WD[k] = cwm_wd_pred %>% filter(YEAR > 3000) %>% select(cwm_wd) %>% colMeans()

}

df %>% write.csv(paste0(output_dir,"/output_vars_614ppm.csv"))


#### Fittest species ####

df_414 = read.csv(paste0(output_dir,"/output_vars_414ppm.csv"), header=T)
df_614 = read.csv(paste0(output_dir,"/output_vars_614ppm.csv"), header=T)

p1 = df_414 %>% ggplot(aes(x=ea, y=eg)) + 
  geom_raster(aes(fill=WD)) + 
  facet_wrap(~ewd, nrow = 1)+
  scale_fill_viridis_c(limits = c(400, 1000))

p2 = df_414 %>% ggplot(aes(x=ea, y=eg)) + 
  geom_raster(aes(fill=df_614$WD-df_414$WD)) + 
  facet_wrap(~ewd, nrow = 1)+
  scale_fill_viridis_c(direction = -1, limits=c(-220,0))


df_scan = read.csv("~/codes/Plant-FATE/fitness_scan.csv")
traits_used = read.csv(file = "~/codes/Plant-FATE/tests/data/Traits_random_HD2.csv") %>% arrange(desc(Total.BasalArea_2017.cm2.))
N = nrow(traits_used)

wd_opt = function(wd, fitness){
  sum(wd*fitness)/sum(fitness)
}

wd_diff = traits_used %>% 
  select(Species, meanWoodDensity..g.cm3.) %>% 
  mutate(i = 1:N) %>% 
  right_join(df_scan) %>% 
  group_by(e_alpha, e_gamma, ewd0) %>% summarize(wd_aCaE = wd_opt(wd = meanWoodDensity..g.cm3.*1e3, fitness = fitness_aCaE^10),
                                                 wd_eCaE = wd_opt(wd = meanWoodDensity..g.cm3.*1e3, fitness = fitness_eCaE^10),
                                                 wd_eCeE = wd_opt(wd = meanWoodDensity..g.cm3.*1e3, fitness = fitness_eCeE^10)) %>% 
  mutate(delta_wd_eCaE = wd_eCaE - wd_aCaE,
         delta_wd_eCeE = wd_eCeE - wd_aCaE)


q1 = wd_diff %>% ggplot(aes(x=e_alpha, y=e_gamma)) + 
  geom_raster(aes(fill=wd_aCaE)) + 
  facet_wrap(~ewd0, nrow = 1)+
  scale_fill_viridis_c(limits = c(400, 1000))

q2 = wd_diff %>% ggplot(aes(x=e_alpha, y=e_gamma)) + 
  geom_raster(aes(fill=delta_wd_eCeE)) + 
  facet_wrap(~ewd0, nrow = 1)+
  scale_fill_viridis_c(direction = -1, limits=c(-220,0))


cairo_pdf(paste0(output_dir, "/fittest_wd_414ppm.pdf"), width = 8, height = 4)
cowplot::plot_grid(p1, q1, nrow = 2, align = "hv")
dev.off()

cairo_pdf(paste0(output_dir, "/delta_wd_614-414ppm.pdf"), width = 8, height = 4)
cowplot::plot_grid(p2, q2, nrow = 2, align = "hv")
dev.off()


# rrdat = dat %>% filter(YEAR %in% (2001+16*(-1:0)))
# rr_short = log(rrdat[2,]/rrdat[1,]) %>% select(-YEAR, -DOY) %>% t() 
# rr_short
