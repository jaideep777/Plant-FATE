library(tidyverse)
rm(list=ls())

output_dir = "pspm_output"
prefix = "ELE"
  
solver = "HD"#_old_params"
setwd(paste0("~/codes/Plant-FATE/",output_dir,"/",prefix,"_",solver))

plot_to_file = T

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

n_species = length(unique(dat2$PID))
n_year = length(unique(dat2$YEAR))
  
par(mfcol=c(3,5), mar=c(6,6,1,1), oma=c(1,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(3.2,1,0), las=1)
seeds = dat2 %>% select(YEAR, PID, SEEDS) %>% spread(value = "SEEDS", key = "PID")
matplot(seeds$YEAR, seeds[,-1], lty=1, col=rainbow(n = n_species+1, start = 0, end = 0.85, alpha = 0.5), type="l",
        las=1, xlab="Time (years)", ylab="Species Seed\noutput", log="")
mtext(line=0.5, side=3, text=solver)

# matplot(seeds1$V1, seeds1[,-1], lty=1, col=rainbow(n = n_species+1, start = 0, end = 0.85), type="l",
#         las=1, xlab="Time (years)", ylab="Species Seed output", log="")
# mtext(line=0.5, side=3, text=solver)

BA = dat2 %>% select(YEAR, PID, BA) %>% spread(value = "BA", key = "PID")
matplot(BA$YEAR, cbind(BA[,-1], rowSums(BA[,-1], na.rm=T))*1e4, lty=1, col=c(rainbow(n = n_species, start = 0, end = 0.85), "black"), type="l",
        las=1, xlab="Time (years)", ylab="Basal area", log="")
abline(h=31.29, col="grey40")


matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")
# matplot(co$V1, co[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
#         las=1, xlab="Time (years)", ylab="Io")
matplot(y=1:24, x=t(-lai_v[,3:26]+lai_v[,2:25]), lty=1, col=rainbow(n = n_year, start = 0, end = 0.85, alpha=0.05), type="l",
        las=1, xlab="Leaf area density", ylab="Height")
matplot(y=1:25, x=t(lai_v[,2:26]), lty=1, col=rainbow(n = n_year, start = 0, end = 0.85, alpha=0.05), type="l",
        las=1, xlab="Cumulative LAI", ylab="Height")


plot(dat$LAI~dat$YEAR, type="l", col="red", ylim=c(0,max(dat$LAI,6.5)), xlab="Time (years)", ylab="Total LAI")
abline(h=c(5.3, 6.2), col=scales::muted("red"))
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
abline(h=c(3,3.5), col=scales::muted("green4"))
abline(h=c(1.31), col=scales::muted("green3"))

matplot(y=cbind(dat$GS), x=dat$YEAR, type="l", lty=1, col=c("cyan3"), ylab="Stomatal conductance\n(mol/m2/s)", xlab="Time (years)")
abline(h=c(0.16), col=scales::muted("cyan3"))

agb = cbind(dat$CL+dat$CW)/1e3
matplot(y=agb, x=dat$YEAR, type="l", lty=1, col=c("yellow3"), ylim=c(0,max(agb)), ylab="AGB\n(kgC/m2)", xlab = "Time (years)")
abline(h=c(16.9, 20.7), col=scales::muted("yellow3"))

matplot(y=cbind(dat$CFR)/1e3, x=dat$YEAR, type="l", lty=1, col=c("brown"), ylab="C-FR\n(kgC/m2)", xlab = "Time (years)", ylim=c(0, max(dat$CFR/1e3,0.7)))
abline(h=c(0.48, 0.66), col=scales::muted("brown"))

matplot(y=cbind(dat$VCMAX), x=dat$YEAR, type="l", lty=1, col=c("green3"), ylab="Vcmax\n(umol/m2/s)", xlab="Time (years)")
abline(h=c(40), col=scales::muted("green3"))

# To get avg size distribution, sum over species and average over years
names(dist)[1:2] = c("YEAR", "SPP")
dist_amb = dist %>% filter(YEAR>1100 & YEAR<2000) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)
dist_ele = dist %>% filter(YEAR>2100 & YEAR<3000) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)
# matplot(y=cbind(as.numeric(colMeans(filter(dist, V1>1100 & V1<2000)[, -c(1,2)], na.rm = T)),
#                 as.numeric(colMeans(filter(dist, V1>2100 & V1<3000)[, -c(1,2)], na.rm = T))
#                 )*1e-2*1e4, 
#         x=x, type="l", log="y", lty=1, col=c("black", "brown"), 
#         xlim=c(0.01, 2), ylim=c(1e-4, 200), ylab="Density (stems/cm/ha)", xlab="Diameter (m)")
matplot(y=cbind(as.numeric(dist_amb[gtools::mixedsort(names(dist_amb))][-1]),
                as.numeric(dist_ele[gtools::mixedsort(names(dist_ele))][-1])
                )*1e-2*1e4, # Convert stems m-1 m-2 --> stems cm-1 ha-1
        x=x, type="l", log="y", lty=1, col=c("black", "brown"), 
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
points(yobs~xobs, col=scales::alpha("blue", 0.5))

traits_obs = read.csv(file = "../../tests/data/Amz_trait_filled_HD.csv")

traits_obs %>% select(Leaf.LMA..g.m2., Total.BasalArea_2017.cm2.) %>% drop_na %>% slice_head(n = 100) %>% with(density(x =Leaf.LMA..g.m2., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.02), las=0, main="", xlab="LMA")
dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR == 2000) %>% 
  with(density(x =LMA*1000, weights=BA/sum(BA))) %>% points(col="red", type="l")

traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% slice_head(n = 100) %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.005), las=0, main="", xlab="Wood density")
dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR == 2000) %>% 
  with(density(x =WD, weights=BA/sum(BA))) %>% points(col="red", type="l")


  # cprops = read.delim("cohort_props.txt")
  # 
  # cprops %>% filter(speciesID==0) %>% filter(t > max(t)-50) %>% with(smoothScatter(lai~height))
  # cprops %>% filter(speciesID==0) %>% filter(t > max(t)-50) %>% with(smoothScatter(gpp~height))
  # cprops %>% filter(speciesID==0) %>% filter(t > max(t)-50) %>% with(smoothScatter(mort~height))
  # 
  # dists = read.delim("size_distributions.txt", header=F)
  # dists %>% filter(V2 == 0) %>% select(-V1, -V2, -V103) %>% as.matrix() %>% t() %>% matplot(x=exp(seq(log(0.01), log(10), length.out=100)), y=., type="l", lty=1, col=rainbow(n = ncol(.), start = 0, end = 0.85, alpha = 10/ncol(.)), log="xy", ylim=c(1e-6,1e2))
  # dists %>% filter(V2 == 1) %>% select(-V1, -V2, -V103) %>% as.matrix() %>% t() %>% matplot(x=exp(seq(log(0.01), log(10), length.out=100)), y=., type="l", lty=1, col=rainbow(n = ncol(.), start = 0, end = 0.85, alpha = 10/ncol(.)), log="xy", ylim=c(1e-6,1e2))
  # dists %>% filter(V2 == 2) %>% select(-V1, -V2, -V103) %>% as.matrix() %>% t() %>% matplot(x=exp(seq(log(0.01), log(10), length.out=100)), y=., type="l", lty=1, col=rainbow(n = ncol(.), start = 0, end = 0.85, alpha = 10/ncol(.)), log="xy", ylim=c(1e-6,1e2))
  # 
  # dists %>% filter(V2 == 2) %>% filter(V1 == 1500) %>% select(-V1, -V2, -V103) %>% as.numeric() %>% plot(x=exp(seq(log(0.01), log(10), length.out=100)), y=., type="l", lty=1, col="black", log="xy", ylim=c(1e-6,1e2))


  
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
  
p1 = traits %>% 
  filter(YEAR > 1050) %>% 
  # filter(RES==T) %>% 
  ggplot(aes(y=LMA, x=WD))+
  theme_classic(base_size = 12)+
  geom_point(aes(col=YEAR, size=RES), alpha=0.4)+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5))


  
# matplot(y=cbind(dat$ET), x=dat$YEAR, type="l", lty=1, col=c("blue"))
  

p2 = dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR > 1120 & YEAR < 2000) %>% 
  # filter(RES==T) %>% 
  ggplot(aes(y=LMA, x=WD))+
  theme_classic(base_size = 12)+
  geom_point(aes(col=BA*1e4, size=RES), alpha=0.7)+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5))+
  labs(col="BA", size=NULL)

p3 = dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR > 2120 & YEAR < 3000) %>% 
  # filter(RES==T) %>% 
  ggplot(aes(y=LMA, x=WD))+
  theme_classic(base_size = 12)+
  geom_point(aes(col=BA*1e4, size=RES), alpha=0.7)+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5))+
  labs(col="BA", size=NULL)

p4 = traits_obs[1:100,] %>%  
  ggplot(aes(y=Leaf.LMA..g.m2., x=meanWoodDensity..g.cm3.))+
  theme_classic(base_size = 12)+
  geom_point(aes(col=Total.BasalArea_2017.cm2./1e4), alpha=0.7)+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5))+
  labs(col="BA", size=NULL)


print(
cowplot::plot_grid(p2,p3,p4, align="hv")
)


