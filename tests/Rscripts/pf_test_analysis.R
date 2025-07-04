library(tidyverse)
rm(list=ls())


input_dir = "~/codes/Plant-FATE/tests/data/"
output_dir = "~/codes/Plant-FATE/pspm_output_test"
# output_dir = "~/output_data/test_3spp_100yr"
expt_dir = "test_2spp_evol_p50" #_old_params"
# expt_dir = "test_3spp_100yr" #_old_params"

# expt_dir = "test_1spp_evol_p50" #_old_params"
# expt_dir = "test_4spp_evol_wd_hmat" #_old_params"
# expt_dir = "cont_test_main" #_old_params"
solver = "IEBT"

source("~/codes/Plant-FATE/R/process_outputs.R")

splinemax = function(x,y, plot=T, ...){
  f = splinefun(x=x, y=y)
  xnew = seq(min(x), max(x), length.out=1000)
  opt = xnew[f(xnew)==max(f(xnew))]

  if(plot){
    plot(y~x, type="p", ...)
    points(f(xnew)~xnew, type="l", col="red")
    abline(v=opt, col="grey")
  }
  opt
}

l1 = pf_read_outputs(input_dir, output_dir, expt_dir)

l1_slice = pf_slice_time(l1, -20000, 20200)

l = l1_slice # pf_subsample_outputs(l1_slice, interval = 11)

setwd(paste0(output_dir,"/",expt_dir))

plot_to_file = F
plot_trait_space = F

add_band = function(end = 30000, start=year(as.Date("2100-1-1"))){
  # polygon(x=c(start,end,end,start), y=c(-1e20,-1e20,1e20,1e20), border = NA, col=scales::alpha("yellow2",0.2))
}

add_hband = function(ylim, col="grey30", alpha=0.2, xlim=c(-1e20,2500)){
  polygon(y=c(ylim[1],ylim[2],ylim[2],ylim[1]), x=c(xlim[1],xlim[1],xlim[2],xlim[2]), border = NA, col=scales::alpha(col, alpha))
}

year_sq = 2000
year_fu = max(l$dat$YEAR)

# traits_obs = read.csv(file = paste0(input_dir, "/Amz_trait_orig.csv"))
# traits_used = read.csv(file = paste0(input_dir, "/Traits_random_HD2.csv"))

# To get avg size distribution, sum over species and average over years
dist_amb = l$dist %>% filter(YEAR > min(YEAR)) %>% filter(YEAR>min(year_sq-500,max(YEAR)-2) & YEAR<year_sq) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)
dist_ele = l$dist %>% filter(YEAR > min(YEAR)) %>% filter(YEAR>min(year_fu-500,max(YEAR)-2) & YEAR<year_fu) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)

# dist_amb = dist %>% filter(YEAR == 1100) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)
# dist_ele = dist %>% filter(YEAR == 1101) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)

n_species = l$dat2 %>% filter(!grepl("probe", .$PID)) %>% pull(PID) %>% unique() %>% length()
n_year = length(unique(l$dat2$YEAR))
col_species = rainbow(n = n_species, start = 0, end = 0.85, alpha = min(10/n_species, 1))
filter_span = span=30/n_year

calib_fluxes =
  read.csv(text = gsub(pattern = " ", replacement="",
                       "name, mean, min, max, type
      GPP,  NA, 3,  3.5, fluxes
      NPP, 1.31, NA, NA, fluxes
      AGB, NA, 16.9, 20.7, structure
      GS, 0.16, NA, NA, fluxes
      LAI, NA, 5.3, 6.2, structure
      CFR, NA, 0.48, 0.66, structure
      VCMAX, NA, 20, 50, fluxes"),
           header=T, sep=",", as.is = F) %>%
  as_tibble()

calib_ba =
  read.csv(text = gsub(pattern = " ", replacement="",
                       "name, mean, min, max, type
      BA,  NA, 30,  35, structure"))
 

calib = calib_fluxes %>% bind_rows(calib_ba)


plot_gauss_mix = function(x, means, sds, wts, add=T, col, add_polygon=F, ...){
  y = x*0
  if (length(sds)==1) sds = rep(sds, length(means))
  for (i in 1:length(means)){
    y = y + wts[i]*dnorm(x, mean = means[i], sd=sds[i])
  }
  if (add) points(y~x, col=col, ...)
  else plot(y~x, col=col, ...)
  if (add_polygon){
    polygon(y~x, col=scales::alpha(col, 0.5), border = NA)
  }
}


if (plot_to_file) png("master_plot.png", width=2412*1.5, height = 1472*1.5, res=300)

par(mfcol=c(4,6), mar=c(4.5,6,.5,1), oma=c(1,1,2,1), cex.lab=1.1, cex.axis=1.1, mgp=c(3.2,1,0), las=1)
seeds = l$dat2 %>% filter(!grepl("probe", .$PID)) %>% select(YEAR, PID, SEEDS) %>% spread(value = "SEEDS", key = "PID")
# seeds_smooth = seeds %>% pivot_longer(-YEAR) %>% group_by(name) %>% mutate(value = loess(value~YEAR, span=filter_span) %>% fitted()) %>% pivot_wider(names_from=name)
seeds_total = rowSums(seeds[,-1,drop=FALSE], na.rm=T)
matplot(seeds$YEAR, cbind(seeds[,-1], seeds_total), lty=1, col=scales::alpha(c(col_species, "black"), 0.7), type="l",
        las=1, xlab="Time (years)", ylab="Species seed\noutput", log="")
# matplot(seeds_smooth$YEAR, seeds_smooth[,-1], lty=1, type="l",
#         # col=col_species,
#         col=scales::alpha(scales::muted(col_species), alpha=min(30/n_species, 1)),
#         las=1, xlab="Time (years)", ylab="Species seed\noutput", log="",
#         add=T
#         )
mtext(line=0.5, side=3, text=expt_dir)
add_band()
# abline(v=3100, col="grey")

# matplot(seeds1$V1, seeds1[,-1], lty=1, col=rainbow(n = n_species+1, start = 0, end = 0.85), type="l",
#         las=1, xlab="Time (years)", ylab="Species Seed output", log="")
# mtext(line=0.5, side=3, text=expt_dir)

BA = l$dat2 %>% filter(!grepl("probe", .$PID)) %>% select(YEAR, PID, BA) %>% spread(value = "BA", key = "PID")
matplot(BA$YEAR, cbind(BA[,-1], rowSums(BA[,-1,drop=FALSE], na.rm=T))*1e4, lty=1, col=c(col_species, "black"), type="l",
        las=1, xlab="Time (years)", ylab="Basal area", log="")
add_hband(c(31.29, 31.29*1.02))
add_band()

matplot(l$Zp$YEAR, l$Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")
add_band()
# matplot(co$V1, co[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
#         las=1, xlab="Time (years)", ylab="Io")
# matplot(y=1:24, x=t(-lai_v[,3:26]+lai_v[,2:25]), lty=1, col=rainbow(n = n_year, start = 0, end = 0.85, alpha=0.05), type="l",
#         las=1, xlab="Leaf area density", ylab="Height")

matplot(y=cbind(l$dat$DPSI), x=l$dat$YEAR, type="l", lty=1, col=c("cyan4"), ylab="Dpsi\n(MPa)", xlab="Time (years)")
matlines(y=cbind(fitted(loess(l$dat$DPSI~l$dat$YEAR, span=filter_span))), x=l$dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
# add_hband(c(20,50)) #, col=scales::muted("green4"))
add_band()


matplot(y=1:25, x=t(l$lai_v[,2:26]), lty=1, col=rainbow(n = n_year, start = 0, end = 0.85, alpha=0.05), type="l",
        las=1, xlab="Cumulative LAI", ylab="Height")


plot(l$dat3$LAI~l$dat3$YEAR, type="l", col="red3", xlab="Time (years)", ylab="Total LAI")#, ylim=c(0,max(l$dat$LAI,8.5)))
matlines(y=cbind(fitted(loess(l$dat3$LAI~l$dat3$YEAR, span=filter_span))), x=l$dat3$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
# abline(h=c(5.3, 6.2), col=scales::muted("red"))
add_hband(c(5.3, 6.2))#, col=scales::alpha("red3", 0.2))
add_band()
# abline(h=c(3.5), col=scales::muted("grey100"))


matplot(y=cbind(l$dat$GPP, l$dat$NPP, l$dat$MORT)*365.2425, x=l$dat$YEAR, type="l", lty=1, col=c("green4", "green3", "brown"), ylab="GPP, NPP, MORT\n(kgC/m2/yr)", xlab="Time (years)")
matlines(y=cbind(fitted(loess(l$dat$GPP~l$dat$YEAR, span=filter_span)),
                 fitted(loess(l$dat$NPP~l$dat$YEAR, span=filter_span)))*365, x=l$dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
# points(y=l$dat$NPP/l$dat$GPP*4, x=l$dat$YEAR, type="l", lty=1, col=c("yellow1"))
# abline(h=c(3,3.5), col="grey")
add_hband(c(3,3.5))#, col=scales::alpha("black",0.3))
add_hband(c(1.31,1.3555))#, col=scales::alpha("black",0.3))
# abline(h=c(1.31), col=scales::muted("green3"))
add_band()


matplot(y=cbind(l$dat$GS), x=l$dat$YEAR, type="l", lty=1, col=c("cyan3"), ylab="Stomatal conductance\n(mol/m2/s)", xlab="Time (years)")
matlines(y=cbind(fitted(loess(l$dat$GS~l$dat$YEAR, span=filter_span))), x=l$dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
add_hband(c(0.16, 0.16555))#, col=scales::alpha("cyan4", 0.6))
# abline(h=c(0.16), col=scales::muted("cyan3"))
add_band()

agb = cbind(l$dat3$CL+l$dat3$CW)
matplot(y=agb, x=l$dat3$YEAR, type="l", lty=1, col=c("yellow4"), ylim=c(0,max(agb)), ylab="AGB\n(kgC/m2)", xlab = "Time (years)")
add_hband(c(16.9, 20.7))#, col=scales::alpha("yellow3", 0.3))
add_band()

matplot(y=cbind(l$dat3$CFR), x=l$dat3$YEAR, type="l", lty=1, col=c("brown"), ylab="C-FR\n(kgC/m2)", xlab = "Time (years)", ylim=c(0, max(l$dat$CFR/1e3,0.7)))
add_hband(c(0.48, 0.66))#, col=scales::alpha("brown", 0.3))
add_band()

matplot(y=cbind(l$dat$VCMAX), x=l$dat$YEAR, type="l", lty=1, col=c("green3"), ylab="Vcmax\n(umol/m2/s)", xlab="Time (years)", ylim=c(0,60))
matlines(y=cbind(fitted(loess(l$dat$VCMAX~l$dat$YEAR, span=filter_span))), x=l$dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
add_hband(c(20,50)) #, col=scales::muted("green4"))
add_band()

plot_size_dist = function(){
  matplot(y=cbind(as.numeric(dist_amb[gtools::mixedsort(names(dist_amb))][-1]),
                  as.numeric(dist_ele[gtools::mixedsort(names(dist_ele))][-1])
  )*1e-2*1e4, # Convert stems m-1 m-2 --> stems cm-1 ha-1
  x=l$x, type="l", log="y", lty=1, col=c("black", "yellow3"),
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

# traits_obs %>% select(Leaf.LMA..g.m2., Total.BasalArea_2017.cm2.) %>% drop_na %>%
#   # with(density(x =Leaf.LMA..g.m2., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.02), las=0, main="", xlab="LMA", col=NA, lwd=2)
#   with(plot_gauss_mix(x=seq(0,300, length.out=1000), means =Leaf.LMA..g.m2., wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.), sds=180*0.1, add=F, col="grey", add_polygon=T, type="l", lwd=1.5,  las=0, main="", xlab="LMA", ylab="density", ylim=c(0, 0.02)))
try(
  l$dat2 %>% select(YEAR, PID, BA) %>%
    left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(!grepl("probe", PID)) %>%
    mutate(YEAR = as.integer(YEAR)) %>%
    mutate(yeardiff = abs(YEAR-year_sq)) %>%
    filter(yeardiff == min(yeardiff)) %>%
    with(plot_gauss_mix(x=seq(0,300, length.out=1000), means =LMA*1000, wts=BA/sum(BA), sds=180*0.14, add=F, col="black", type="l", lwd=1.5, xlab="LMA"))
)

# traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>%
#   #with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.005), las=0, main="", xlab="Wood density", col=NA, lwd=2)
#   with(plot_gauss_mix(x=seq(200,1200, length.out=1000), means =meanWoodDensity..g.cm3.*1000, wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.), sds=800*0.1, add=F, col="grey", add_polygon=T, type="l", lwd=1.5,  las=0, main="", xlab="Wood density", ylab="density", ylim=c(0, 0.005)))
try(
  l$dat2 %>% select(YEAR, PID, BA) %>%
    left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(!grepl("probe", PID)) %>%
    mutate(YEAR = as.integer(YEAR)) %>%
    mutate(yeardiff = abs(YEAR-year_sq)) %>%
    filter(yeardiff == min(yeardiff)) %>%
    # with(density(x =WD, adjust=1, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
    with(plot_gauss_mix(x=seq(200,1200, length.out=1000), means =WD, wts=BA/sum(BA), sds=800*0.14, add=F, col="black", type="l", lwd=1.5, xlab="Wood density"))
)

# traits_obs %>% select(Height_Max.m., Total.BasalArea_2017.cm2.) %>% drop_na %>%
#   # with(density(x =Height_Max.m., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.18), las=0, main="", xlab="Max. height", col=NA, lwd=2)
#   with(plot_gauss_mix(x=seq(0,50, length.out=1000), means =Height_Max.m., wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.), sds=25*0.1, add=F, col="grey", add_polygon=T, type="l", lwd=1.5,  las=0, main="", xlab="Max. Height", ylab="density", ylim=c(0, 0.2)))
try(
  l$dat2 %>% select(YEAR, PID, BA) %>%
    left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(!grepl("probe", PID)) %>%
    mutate(YEAR = as.integer(YEAR)) %>%
    mutate(yeardiff = abs(YEAR-year_sq)) %>%
    filter(yeardiff == min(yeardiff)) %>%
    # with(density(x =HMAT, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
    with(plot_gauss_mix(x=seq(0,50, length.out=1000), means =HMAT, wts=BA/sum(BA), sds=25*0.14, add=F, col="black", type="l", lwd=1.5, xlab="Max. height"))
)

# traits_obs %>% select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% drop_na %>%
#   # with(density(x =Height_Max.m., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.18), las=0, main="", xlab="Max. height", col=NA, lwd=2)
#   with(plot_gauss_mix(x=seq(-6,-0.1, length.out=1000), means =P50..Mpa., wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.), sds=2*0.1, add=F, col="grey", add_polygon=T, type="l", lwd=1.5,  las=0, main="", xlab="P50 (MPa)", ylab="density", ylim=c(0, 1.5)))
try(
  l$dat2 %>% select(YEAR, PID, BA) %>%
    left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(!grepl("probe", PID)) %>%
    mutate(YEAR = as.integer(YEAR)) %>%
    mutate(yeardiff = abs(YEAR-year_sq)) %>%
    filter(yeardiff == min(yeardiff)) %>%
    # with(density(x =HMAT, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
    with(plot_gauss_mix(x=seq(-6,-0.1, length.out=1000), means =P50X, wts=BA/sum(BA), sds=2*0.31, add=F, col="black", type="l", lwd=1.5, xlab="P50"))
)


#
# traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.005), las=0, main="", xlab="Wood density", col=NA, lwd=2)
# traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
# l$dat2 %>% select(YEAR, PID, BA) %>%
#   filter(YEAR == min(2000, max(l$dat2$YEAR)-1)) %>%
#   left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>%
#   left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
#   # drop_na %>%
#   with(density(x =meanWoodDensity..g.cm3.*1000, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
# try(
#   l$dat2 %>% select(YEAR, PID, BA) %>%
#     filter(YEAR == min(2000, max(l$dat2$YEAR)-1)) %>%
#     left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>%
#     left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
#     #drop_na %>%
#     with(density(x =WD, weights=BA/sum(BA))) %>% points(col="yellow3", type="l", lwd=1.5)
# )
#
#

# cwm_wd = traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% summarise(cwm=sum(meanWoodDensity..g.cm3.*Total.BasalArea_2017.cm2.)/sum(Total.BasalArea_2017.cm2.))*1000
cwm_wd_pred = l$dat2 %>% select(YEAR, PID, BA) %>%
  left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
  drop_na %>% group_by(YEAR) %>%
  summarize(cwm_wd = sum(WD*BA)/sum(BA))

cwm_wd_pred %>% with(plot(cwm_wd~YEAR, type="l", ylim=c(580,840))) #, ylim=c(200,900)))
# add_hband(c(cwm_wd,cwm_wd+5))

# l$traits %>% select(YEAR, SPP, HMAT) %>% pivot_wider(names_from = "SPP", values_from = "HMAT") %>% with(matplot(x=.[,1], y=.[,2:5], lty=1, type="l"))
# l$traits %>% select(YEAR, SPP, WD) %>% pivot_wider(names_from = "SPP", values_from = "WD") %>% with(matplot(x=.[,1], y=.[,2:5], lty=1, type="l"))

hmat = l$traits %>% filter(!grepl("probe", .$SPP)) %>% select(YEAR, SPP, HMAT) %>% pivot_wider(names_from = "SPP", values_from = "HMAT")
wd = l$traits %>% filter(!grepl("probe", .$SPP)) %>% select(YEAR, SPP, WD) %>% pivot_wider(names_from = "SPP", values_from = "WD")
p50x = l$traits %>% filter(!grepl("probe", .$SPP)) %>% select(YEAR, SPP, P50X) %>% pivot_wider(names_from = "SPP", values_from = "P50X")

# matplot(y=hmat[,-1], x=wd[,-1], lty=1, type="o", pch=20, cex=0.5, col=col_species, xlab="Wood density", ylab="Max height")

matplot(x=wd[,1], y=wd[,-1], col=col_species, lty=1, type="l", ylab="Wood density", xlab="Year")
add_band()
matplot(x=hmat[,1], y=hmat[,-1], col=col_species, lty=1, type="l", ylab="Max. height", xlab="Year")
add_band()
matplot(x=p50x[,1], y=p50x[,-1], col=col_species, lty=1, type="l", ylab="P50x", xlab="Year")
add_band()


exp_fit = function(y,x, r=1e-4){
  y0 = first(y)
  x0 = first(x)
  # plot(y~x)
  a = last(y) # -1.5
  curve(y0 - (y0 - a)*(1-exp(-r*(x-x0))), add=T, col="red")
  points(I(y0 - (y0 - a)*(1-exp(-r*(x-x0))))~x, type="l", col="red")
  abline(v=x0, col="grey")
  mod = try(nls(formula = y~I(y0 - (y0 - a)*(1-exp(-r*(x-x0)))), start = list(a=a, r=r), algorithm = "default", trace=T)) # , upper = list(a=-0.1, r=1), lower=list(a=-6, r=1e-5))
  points(fitted(mod)~x, type="l", col="green", lwd=2) #, ylim=c(200,900)))
  a_opt = summary(mod)$parameters[,1]["a"]
  r_opt = summary(mod)$parameters[,1]["r"]
  a_opt
}


trait_dist = l$dat2 %>%
  select(YEAR, PID, BA, SEEDS) %>%
  mutate(BA=BA*1e4) %>%
  group_by(PID) %>%
  summarize(
    dSEEDS = lm(tail(SEEDS,100)~tail(YEAR,100))$coefficients[2],
    dBA = lm(tail(BA,10)~tail(YEAR,10))$coefficients[2],
    BA = last(BA),
    SEEDS = last(SEEDS),
  ) %>%
  ungroup() %>%
  left_join(l$traits %>% filter(YEAR == last(YEAR)), by = c("PID"="SPP")) %>%
  drop_na()

trait_dist_y = l$dat2 %>%
  select(YEAR, PID, BA, SEEDS) %>%
  mutate(BA=BA*1e4) %>%
  mutate(yr_int = cut(YEAR, breaks = seq(first(YEAR), last(YEAR), by=1000))) %>%
  group_by(PID, yr_int) %>%
  summarize(
    dSEEDS = lm(SEEDS~YEAR)$coefficients[2],
    dBA = lm(BA~YEAR)$coefficients[2],
    BA = last(BA),
    SEEDS = last(SEEDS),
    YEAR = first(YEAR)
  ) %>%
  left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
  drop_na()

# trait_dist_y %>% ggplot() +
#   geom_raster(aes(x=YEAR, y=WD, fill=dSEEDS))+
#   scale_fill_gradient2()

st_yr = last(l$dat$YEAR) - 10000
with(cwm_wd_pred, plot(y=cwm_wd, x= YEAR))
a_opt = try(
  cwm_wd_pred %>% filter(YEAR>st_yr) %>% with(exp_fit(cwm_wd, YEAR))
)
try(
  cat("cwm_wd* = ", a_opt, '\n')
)

wdopt_dseeds_y = 
try(
  trait_dist_y %>%
  group_by(YEAR) %>%
  summarize(
    WDopt_dseeds = splinemax(x=WD, y=dSEEDS, plot = F, ylab="dseeds_dt", xlab="wd")
  )
)

wdopt_dseeds = trait_dist %>% with(splinemax(x=WD, y=dSEEDS, ylab="dseeds_dt", xlab="wd"))
wdopt_seeds  = trait_dist %>% with(splinemax(x=WD, y=SEEDS, ylab="seeds", xlab="wd", plot=F))

wdopt_dseeds_y %>% with(plot(WDopt_dseeds~YEAR, type="o", ylab="Opt WD"))
wd_opt_dsseds_trend = wdopt_dseeds_y %>% filter(YEAR > st_yr) %>% with(exp_fit(WDopt_dseeds, YEAR))

wdopt_dba = trait_dist %>% with(splinemax(x=WD, y=dBA, ylab="dba_dt", xlab="wd", plot=F))
wdopt_ba  = trait_dist %>% with(splinemax(x=WD, y=BA, ylab="ba", xlab="wd", plot=F))

cat("Opt WD (dseeds/dt) = ", wdopt_dseeds, "\n")
cat("Opt WD (dseeds/dt Trend) = ", wd_opt_dsseds_trend, "\n")
cat("Opt WD (seeds) = ", wdopt_seeds, "\n")

cat("Opt WD (dba/dt) = ", wdopt_dba, "\n")
# cat("Opt WD (ba) = ", wdopt_ba, "\n")

# cwm_p50 = traits_obs %>% select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% drop_na %>% summarise(cwm=sum(P50..Mpa.*Total.BasalArea_2017.cm2.)/sum(Total.BasalArea_2017.cm2.))
# cwm_p50_pred = l$dat2 %>% select(YEAR, PID, BA) %>%
#   left_join(traits_obs, by = c("PID"="Species")) %>%
#   drop_na %>% group_by(YEAR) %>%
#   summarize(cwm_p50 = sum(P50..Mpa.*BA)/sum(BA))
#
# cwm_p50_pred %>% with(plot(cwm_p50~YEAR, type="l")) #, ylim=c(200,900)))
# add_hband(c(cwm_wd,cwm_wd+5))


# l$traits %>%
#   select(YEAR, SPP, r0_cesaro) %>%
#   pivot_wider(names_from = SPP, values_from = r0_cesaro) %>%
#   as.matrix() %>%
#   matplot(y=.[,-1], x=.[,1], type="l", lty=1, ylab="spp1_r0")
zz = l$traits %>%
  select(YEAR, SPP, r0_avg) %>%
  filter(!grepl(x = SPP, "probe")) %>%
  pivot_wider(names_from = SPP, values_from = r0_avg) %>%
  as.matrix()
matplot(y=tanh(zz[,-1]*20)/20, x=zz[,1], type="l", lty=1, ylab="r0", col=col_species)
abline(h=0, col="black", lwd=0.2)
# abline(v=1000, col="grey")

# Spp 1 Fitness gradient plots
# l$traits %>%
#   select(YEAR, SPP, r0_avg) %>%
#   pivot_wider(names_from = SPP, values_from = r0_avg) %>%
#   with(plot(I((`1_probe0`-`1`)/.001)~YEAR, type="l"))
# abline(h=0, col="red")
#
# l$traits %>%
#   select(YEAR, SPP, r0_avg) %>%
#   pivot_wider(names_from = SPP, values_from = r0_avg) %>%
#   with(plot(I((`1_probe1`-`1`)/.001)~YEAR, type="l"))
# abline(h=0, col="red")

# l$traits %>%
#   select(YEAR, SPP, r0_avg) %>%
#   pivot_wider(names_from = SPP, values_from = r0_avg) %>%
#   with(plot(I((`3_probe1`-`3`)/.001)~YEAR, type="l"))
# abline(h=0, col="red")




if (plot_to_file) dev.off()

#### Calib ####


#### Fittest traits

l$dat2 %>% select(YEAR, PID, BA) %>%
  left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
  mutate(BA = BA*1e4) %>%
  filter(YEAR == 2000) %>%
  arrange(desc(BA)) %>%
  slice(1:5)


#### TRAIT SPACE ####
if (plot_trait_space){

  if (plot_to_file) png("traitspace_HMAT_WD.png", width=1618, height = 1196, res=300)
  print(
    l$traits %>%
      left_join(l$dat2 %>% select(YEAR,PID,BA),
                by=c("SPP"="PID", "YEAR"="YEAR")) %>%
      # filter(YEAR > 1050) %>%
      filter(RES==T) %>%
      ggplot(aes(y=HMAT, x=WD))+
      theme_classic(base_size = 12)+
      geom_point(aes(size=YEAR, col=BA*1e4), alpha=0.4)+
      scale_color_viridis_c(direction = -1)+
      scale_size("YEAR", range = c(0, 1.5))
  )
  dev.off()

  # matplot(y=cbind(l$dat$ET), x=l$dat$YEAR, type="l", lty=1, col=c("blue"))


  p2 = l$dat2 %>% select(YEAR, PID, BA) %>%
    left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(YEAR > 1120 & YEAR < 2000) %>%
    # filter(RES==T) %>%
    ggplot(aes(y=LMA, x=WD))+
    theme_classic(base_size = 12)+
    geom_point(aes(col=BA*1e4, size=RES), alpha=0.7)+
    scale_color_viridis_c(direction = -1)+
    scale_size("size_RES", range = c(0, 1.5), guide = F)+
    labs(col="BA")+
    ggtitle("Ambient")

  p3 = l$dat2 %>% select(YEAR, PID, BA) %>%
    left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
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

  with(l$dat %>% filter(YEAR<1200), matplot(y=cbind(GPP, NPP)*1e-3*365, x=YEAR, type="l", lty=1, col=c("green4", "green3"), ylab="GPP, NPP\n(kgC m-2 yr-1)", xlab="Time (years)"))
  # abline(h=c(3,3.5), col="grey")
  add_hband(c(3,3.5), col=scales::alpha("black",0.2))
  add_hband(c(1.31,1.4), col=scales::alpha("black", 0.4))#, col=scales::alpha("black",0.3))
  # abline(h=c(1.31), col=scales::muted("green3"))
  mtext(text = "CO2 Fluxes", side=3, line=1)

  traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.005), las=0, main="", xlab="Wood density", col=NA, lwd=2)
  traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
  l$dat2 %>% select(YEAR, PID, BA) %>%
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


## Top 5 species by BA

l$dat2 %>% select(YEAR, PID, BA) %>%
  mutate(BA=BA*1e4) %>%
  left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
  filter(YEAR %in% c(max(YEAR)-1, round((max(YEAR)+min(YEAR))/2))) %>% arrange(desc(BA)) %>%
  group_by(YEAR) %>% slice(1:5)

l$dat2 %>% select(YEAR, PID, BA, SEEDS) %>%
  mutate(BA=BA*1e4) %>%
  left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
  filter(YEAR %in% c(as.integer(max(YEAR)/2), max(YEAR))) %>% arrange(desc(BA)) %>% group_by(YEAR) %>% slice(1:5)

l$dat2 %>% select(YEAR, PID, BA, SEEDS) %>%
  mutate(BA=BA*1e4) %>%
  left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
  filter(YEAR %in% c(as.integer(max(YEAR)/2), max(YEAR))) %>%
  ggplot(aes(x=WD, y=BA, group=YEAR, col=YEAR)) +
  geom_line()


# l$dat2 %>% select(YEAR, PID, BA) %>%
#   mutate(BA=BA*1e4) %>%
#   left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
#   #filter(YEAR %in% c(max(YEAR)-1, round((max(YEAR)+min(YEAR))/2))) %>% arrange(desc(BA)) %>%
#   group_by(YEAR) %>% slice(1:5) %>%
#   select(YEAR, PID, BA) %>%
#   spread(key = PID, value = BA) %>%
#   with(matplot(.[,1], .[,-1], type="l", lty=1))

# l$dat2 %>% select(YEAR, PID, BA) %>%
#   mutate(BA=BA*1e4) %>%
#   left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
#   filter(YEAR %in% c(max(YEAR)-1, round((max(YEAR)+min(YEAR))/2))) %>%
#   arrange(desc(BA)) %>% View()


## Dead species in first 400 years

traits_used %>% filter(Species %in% colnames(seeds)[seeds %>% slice(400) %>% is.na() %>% which()])

## DEBUG

# par(mfrow=c(1,1))
# # Plot GPP vs total seed production
# matplot(y=cbind(l$dat$GPP*365, seeds_total/900, ((seeds-seeds_smooth)[-1]+mean(seeds_total))/900), col=c("green3", "black", col_species), x=l$dat$YEAR, type="l", lty=1, lwd=c(2,2,rep(1, n_species)), ylab = "GPP/seeds", xlab="Time (years)")
# # plot(x=l$dat$GPP*365, y=seeds_total/900)
# plot(l$dat2$YEAR)

# # #### Cohort props ####
# #
# df_cohorts = readr::read_csv(paste0(output_dir,"/",expt_dir,"/cohort_props.csv"))
# #
# # library(patchwork)
# #
# df_cohorts %>%
#   filter (YEAR < 2500) %>%
#   # mutate(YRSEG = cut(YEAR, c(-Inf,-19800,1800,2000,19800,Inf), labels = c("start","trait_spin", "hist", "trait_spin_ele", "ele"))) %>%
#   # filter(YRSEG %in% c("start")) %>%
#   filter(cohortID != 0) %>%
#   ggplot(aes(x=YEAR, y=height, colour = cohortID, group=cohortID))+
#   geom_line(alpha=0.5)+
#   theme_bw()+
#   scale_color_viridis_c()
# # facet_wrap(~speciesID+YRSEG, scales="free_x")
