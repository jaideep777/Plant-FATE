library(tidyverse)
rm(list=ls())

source("~/codes/Plant-FATE/tests/Rscripts/utils.R")

#' @param xfrac The fraction over from the left side.
#' @param yfrac The fraction down from the top.
#' @param label The text to label with.
#' @param pos Position to pass to text()
#' @param ... Anything extra to pass to text(), e.g. cex, col.
add_label <- function(xfrac = 0.0, yfrac = 0.1, label, pos = 4, ...) {
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, col="grey40", cex=1.2, ...)
}


output_dir = "~/codes/Plant-FATE/pspm_output_alloc_change"
# output_dir = "~/output_data/pspm_output_36sims"
expt_dir = "zeta_0.2_to_0.400000" #_old_params"

setwd(paste0(output_dir,"/",expt_dir))

plot_to_file = T
plot_trait_space = F

add_band = function(){
  polygon(x=c(2000,5000,5000,2000), y=c(-1e20,-1e20,1e20,1e20), border = NA, col=scales::alpha("yellow2",0.2))
}

add_hband = function(ylim, col=scales::alpha("grey30",0.3), xlim=c(-1e20,2000)){
  polygon(y=c(ylim[1],ylim[2],ylim[2],ylim[1]), x=c(xlim[1],xlim[1],xlim[2],xlim[2]), border = NA, col=col)
}

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

# To get avg size distribution, sum over species and average over years
names(dist)[1:2] = c("YEAR", "SPP")
dist_amb = dist %>% filter(YEAR>1100 & YEAR<2000) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)
dist_ele = dist %>% filter(YEAR>4100 & YEAR<5000) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)

# dist_amb = dist %>% filter(YEAR == 1100) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)
# dist_ele = dist %>% filter(YEAR == 1101) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)

n_species = length(unique(dat2$PID))
n_year = length(unique(dat2$YEAR))

eco2_col = rgb(190/255,190/255,0)


if (plot_to_file) cairo_pdf("../../paper_figs3/emg_props.pdf", width=7*1.1, height = 4.5*1.1)


# par(mfrow=c(3,3), mar=c(4.5,5,.5,1), oma=c(3,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(2.6,1,0), las=0)

par(mfrow=c(2,3), mar=c(4.5,6,.5,1), oma=c(3,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(3.5,1,0), las=1)

matplot(y=cbind(dat$GPP, dat$NPP)*1e-3*365, x=dat$YEAR, type="l", lty=1, col=scales::alpha("red4", 0.2), ylab="GPP, NPP\n(kgC m-2 yr-1)", xlab="Year")
matlines(y=cbind(fitted(loess(dat$GPP~dat$YEAR, span=0.006)), 
                 fitted(loess(dat$NPP~dat$YEAR, span=0.006)))*1e-3*365, x=dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
# points(y=dat$NPP/dat$GPP*4, x=dat$YEAR, type="l", lty=1, col=c("yellow1"))
# abline(h=c(3,3.5), col="grey")
add_hband(c(3,3.5))#, col=scales::alpha("black",0.3))
add_hband(c(1.31,1.3555))#, col=scales::alpha("black",0.3))
# abline(h=c(1.31), col=scales::muted("green3"))
add_band()
add_label(label = "A")

matplot(y=cbind(dat$GS), x=dat$YEAR, type="l", lty=1, col=scales::alpha(c("red4"),0.2), ylab="Stomatal conductance\n(mol m-2 s-1)", xlab="Year")
matlines(y=cbind(fitted(loess(dat$GS~dat$YEAR, span=0.006))), x=dat$YEAR, type="l", lty=1, col="black", ylab="GPP, NPP\n(kgC m-2 yr-1)", xlab="Year")
add_hband(c(0.16, 0.16555))#, col=scales::alpha("cyan4", 0.6))
# abline(h=c(0.16), col=scales::muted("cyan3"))
add_band()
add_label(label = "B")


matplot(y=cbind(dat$VCMAX), x=dat$YEAR, type="l", lty=1, col=scales::alpha(c("red4"), 0.2), ylab="Vcmax\n(umol m-2 s-1)", xlab="Year")
matlines(y=cbind(fitted(loess(dat$VCMAX~dat$YEAR, span=0.006))), x=dat$YEAR, type="l", lty=1, col="black", ylab="GPP, NPP\n(kgC m-2 yr-1)", xlab="Year")
add_hband(c(20,45)) #, col=scales::muted("green4"))
add_band()
add_label(label = "C")



plot(dat$LAI~dat$YEAR, type="l", col="black", ylim=c(0,max(dat$LAI,6.5)), xlab="Year", ylab="Total LAI")
# abline(h=c(5.3, 6.2), col=scales::muted("red"))
add_hband(c(5.3, 6.2))#, col=scales::alpha("red3", 0.2))
add_band()
# abline(h=c(3.5), col=scales::muted("grey100"))
add_label(label = "D")



BA = dat2 %>% select(YEAR, PID, BA) %>% spread(value = "BA", key = "PID")
col_spp = sample(rainbow(n = n_species, start = 0, end = 0.9, alpha = 0.4), replace = F)
# col_spp = rep(scales::alpha("blue", 0.2), n_species)
matplot(BA$YEAR, cbind(BA[,-1], rowSums(BA[,-1], na.rm=T))*1e4, lty=1, col=c(col_spp, "black"), type="l",
        las=1, xlab="Year", ylab="Basal area (m2)", log="")
add_hband(c(31.29, 31.29*1.02))
add_band()
add_label(label = "E")


matplot(Zp[,1], Zp[,-1], lty=1, col=scales::viridis_pal(direction = -1, end=0.95)(4), type="l",
        las=1, xlab="Year", ylab="Canopy layer\nheights (m)", log="")
add_band()
add_label(label = "F")

dev.off()


##### Structural Calib ####


cols_amb_ele = scales::viridis_pal(begin = 0.3, end = 0.8)(3)[c(1,3)]

if (plot_to_file) cairo_pdf("../../paper_figs3/emg_props_structural.pdf", width=6*1.1, height = 5.5*1.1)

# par(mfrow=c(3,3), mar=c(4.5,5,.5,1), oma=c(3,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(2.6,1,0), las=0)
par(mfrow=c(2,2), mar=c(4.5,5,.5,1), oma=c(3,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(2.6,1,0), las=1)

matplot(y=cbind(as.numeric(dist_amb[gtools::mixedsort(names(dist_amb))][-1]),
                as.numeric(dist_ele[gtools::mixedsort(names(dist_ele))][-1])
)*1e-2*1e4, # Convert stems m-1 m-2 --> stems cm-1 ha-1
x=x, type="l", log="y", lty=1, col=cols_amb_ele, lwd=2,
xlim=c(0.01, 1.2), ylim=c(1e-4, 1000), ylab="Density\n(stems cm-1 ha-1)", xlab="Diameter (m)", las=0)
# abline(v=1, col=scales::alpha("red", 0.2))
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
add_label(label = "A", xfrac = 0.85, yfrac=-50)



traits_obs = read.csv(file = "../../tests/data/Amz_trait_filled_HD.csv")

traits_obs %>% select(Leaf.LMA..g.m2., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =Leaf.LMA..g.m2., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.03), las=0, main="", xlab="LMA (kg m-2)", col=NA, lwd=2)
traits_obs %>% select(Leaf.LMA..g.m2., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =Leaf.LMA..g.m2., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
  with(density(x =LMA*1000, weights=BA/sum(BA))) %>% points(col=cols_amb_ele[1], type="l", lwd=2)
try(
  dat2 %>% select(YEAR, PID, BA) %>% 
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
    filter(YEAR == 5000) %>% 
    with(density(x =LMA*1000, weights=BA/sum(BA))) %>% points(col=cols_amb_ele[2], type="l", lwd=2)
)
add_label(label = "B")


traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.006), las=0, main="", xlab="Wood density (kg m-3)", col=NA, lwd=2)
traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
  with(density(x =WD, weights=BA/sum(BA))) %>% points(col=cols_amb_ele[1], type="l", lwd=2)
try(
  dat2 %>% select(YEAR, PID, BA) %>% 
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
    filter(YEAR == 5000) %>% 
    with(density(x =WD, weights=BA/sum(BA))) %>% points(col=cols_amb_ele[2], type="l", lwd=2)
)
add_label(label = "C")

traits_obs %>% select(Height_Max.m., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =Height_Max.m., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.25), las=0, main="", xlab="Max. height", col=NA, lwd=2)
traits_obs %>% select(Height_Max.m., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =Height_Max.m., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
try(
  dat2 %>% select(YEAR, PID, BA) %>% 
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
    filter(YEAR == min(2000, max(dat2$YEAR)-1)) %>% 
    with(density(x =HMAT, weights=BA/sum(BA))) %>% points(col=cols_amb_ele[1], type="l", lwd=2)
)
try(
  dat2 %>% select(YEAR, PID, BA) %>%
    left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(YEAR == min(5000, max(dat2$YEAR)-1)) %>%
    with(density(x =HMAT, weights=BA/sum(BA))) %>% points(col=cols_amb_ele[2], type="l", lwd=2)
)
add_label(label = "D")



if (plot_to_file) dev.off()




##### 

if (plot_to_file) cairo_pdf("../../paper_figs3/emg_props_si.pdf", width=7, height = 5)

par(mfrow=c(2,3), mar=c(6,6,1,1), oma=c(1,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(3.2,1,0), las=1)


# matplot(y=1:24, x=t(-lai_v[,3:26]+lai_v[,2:25]), lty=1, col=rainbow(n = n_year, start = 0, end = 0.85, alpha=0.05), type="l",
#         las=1, xlab="Leaf area density", ylab="Height")
matplot(y=1:25, x=t(lai_v[,2:26]), lty=1, col=rainbow(n = n_year, start = 0, end = 0.85, alpha=0.05), type="l",
        las=1, xlab="Cumulative LAI", ylab="Height (m)")



matplot(Zp$V1, Zp[,-1], lty=1, col=scales::viridis_pal(end = 0.95, direction = -1)(5), type="l",
        las=1, xlab="Year", ylab="Canopy layer\nheights (m)")
add_band()

matplot(co$V1, co[,-1]*100, lty=1, col=scales::viridis_pal(end = 0.95, direction = -1)(5), type="l",
        las=1, xlab="Year", ylab="PAR in layer (%)", ylim=c(0,100))

# matplot(BA1$V1, cbind(BA1[,-1], rowSums(BA1[,-1], na.rm=T))*1e4, lty=1, col=c(rainbow(n = n_species+1, start = 0, end = 0.85), "black"), type="l",
#         las=1, xlab="Year", ylab="Basal area", log="")
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



agb = cbind(dat$CL+dat$CW)/1e3
matplot(y=agb, x=dat$YEAR, type="l", lty=1, col=c("yellow4"), ylim=c(0,max(agb)), ylab="AGB\n(kgC m-2)", xlab = "Year")
add_hband(c(16.9, 20.7))#, col=scales::alpha("yellow3", 0.3))
add_band()


matplot(y=cbind(dat$CFR)/1e3, x=dat$YEAR, type="l", lty=1, col=c("brown"), ylab="C-FR\n(kgC m-2)", xlab = "Year", ylim=c(0, max(dat$CFR/1e3,0.7)))
add_hband(c(0.48, 0.66))#, col=scales::alpha("brown", 0.3))
add_band()


if (plot_to_file) dev.off()


BA_change = rowSums(BA[c(1000,2000),-1], na.rm=T)*1e4
BA_change[3] = (BA_change[2]-BA_change[1])/BA_change[1]*100

P_change = t(cbind(dat$GPP*1e-3*365, 
                   dat$NPP*1e-3*365, 
                   dat$RAU*1e-3*365,
                   rowSums(BA[,-1], na.rm=T)*1e4,
                   (dat$CL)*1e-3,
                   (dat$CW+dat$CL)*1e-3,
                   (dat$CCR+dat$CFR)*1e-3)[c(1000,2000),]) %>% as.data.frame()
rownames(P_change) = c("GPP", "NPP", "RAU", "BA", "CL", "AGB", "BGB")
colnames(P_change) = c("AMB", "ELE")
P_change$Delta = (P_change[,2]-P_change[,1])/P_change[,1]*100
P_change

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


#### TRAIT SPACE ####
if (plot_trait_space){

  
p1 = traits %>% 
  filter(YEAR > 1050) %>% 
  # filter(RES==T) %>% 
  ggplot(aes(y=LMA, x=WD))+
  theme_classic(base_size = 12)+
  geom_point(aes(col=YEAR, size=RES), alpha=0.4)+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5))
# print(p1)

  
# matplot(y=cbind(dat$ET), x=dat$YEAR, type="l", lty=1, col=c("blue"))
  

p2 = dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR > 1120 & YEAR < 2000) %>% 
  # filter(RES==T) %>% 
  ggplot(aes(y=LMA, x=WD))+
  theme_classic(base_size = 12)+
  geom_point(aes(col=BA*1e4, size=RES), alpha=0.7)+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5), guide = F)+
  labs(col="BA")+
  ggtitle("Ambient")

p3 = dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR > 2120 & YEAR < 3000) %>% 
  # filter(RES==T) %>% 
  ggplot(aes(y=LMA, x=WD))+
  theme_classic(base_size = 12)+
  geom_point(aes(col=BA*1e4, size=RES), alpha=0.7)+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5), guide = F)+
  labs(col="BA")+
  ggtitle("Elevated")

p4 = traits_obs[1:100,] %>%  
  ggplot(aes(y=Leaf.LMA..g.m2., x=meanWoodDensity..g.cm3.))+
  theme_classic(base_size = 12)+
  geom_point(aes(col=Total.BasalArea_2017.cm2./1e4), alpha=0.7)+
  scale_color_viridis_c(direction = -1)+
  scale_size("size_RES", range = c(0, 1.5), guide = F)+
  labs(col="BA")+
  ggtitle("Observed")

if (plot_to_file) png("traitspace.png", width=1618*1.5, height = 1196*1.5, res=300)
print(
cowplot::plot_grid(p2,p3,p4, align="hv")
)
if (plot_to_file) dev.off()

}


#### Sample results  ####
plot_sample=F

if (plot_sample){
par(mfrow=c(1,3), mar=c(5,6,4,1), oma=c(1,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(3.2,1,0), las=1)

with(dat %>% filter(YEAR<1200), matplot(y=cbind(GPP, NPP)*1e-3*365, x=YEAR, type="l", lty=1, col=c("green4", "green3"), ylab="GPP, NPP\n(kgC m-2 yr-1)", xlab="Year"))
# abline(h=c(3,3.5), col="grey")
add_hband(c(3,3.5), col=scales::alpha("black",0.2))
add_hband(c(1.31,1.4), col=scales::alpha("black", 0.4))#, col=scales::alpha("black",0.3))
# abline(h=c(1.31), col=scales::muted("green3"))
mtext(text = "CO2 Fluxes", side=3, line=1)

traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.005), las=0, main="", xlab="Wood density", col=NA, lwd=2)
traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
dat2 %>% select(YEAR, PID, BA) %>%
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
  filter(YEAR == 2000) %>%
  with(density(x =WD, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
mtext(text = "Sample\ntrait distribution", side=3, line=1)


matplot(y=cbind(as.numeric(dist_amb[gtools::mixedsort(names(dist_amb))][-1])
)*1e-2*1e4, # Convert stems m-1 m-2 --> stems cm-1 ha-1
x=x, type="l", log="y", lty=1, col=c("black", "yellow3"),
xlim=c(0.01, 1.2), ylim=c(1e-4, 1000), ylab="Density\n(stems cm-1 ha-1)", xlab="Diameter (m)", las=0)

# abline(v=1, col=scales::alpha("red", 0.2))

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
mtext(text = "Size distribution", side=3, line=1)
}




