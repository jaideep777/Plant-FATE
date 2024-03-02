rm(list=ls())

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(scales)
library(zoo)
library(rphydro)

setwd("~/codes/phydro/vignettes/")

dat = read.csv(file="data/drying_experiments_meta-analysis_Joshi_et_al_2022.csv")
names = colnames(dat)
# [4] "Predawn.LWP..MPa."                                 
# [5] "A..umol.m.2.s.1."                                  
# [6] "gC..mol.m.2.s.1."                                  
# [7] "ca..ppm."                                          
# [8] "D..unitless...Pa...Pa.atm.."                       
# [9] "T..deg.C."                                                          
names[4:9] = c("LWP", "A", "gC", "ca", "D", "T")
names[13] = "Iabs_growth"
names[12] = "Iabs_used"
colnames(dat) = names

dpsi_df = read.csv(file = "data/drying_experiments_dpsi_Joshi_et_al_2022.csv")
par_data = read.csv("data/fitted_params_Joshi_et_al_2022.csv") 


source("utils.R")

plot_all = function(df_w_vol, varname, species, data, dpsi_data=NULL, analytical=F){
  # df_w_vol = df_w_vol[complete.cases(df_w_vol),]
  # View(df_w_vol)
  
  gx = log(df_w_vol$gs+1e-20)
  gy = df_w_vol$var
  # gx1 = seq(max(min(gx), -20), max(gx), length.out=100)
  f = splinefun(x = gx, y=gy)
  gs0 = df_w_vol$gs[which(df_w_vol$var==0)]
  psi88S = f(log(gs0*0.12))
  psi50S = f(log(gs0*0.50))
  psi12S = f(log(gs0*0.88))
  cat("psi88S = ", psi88S, "\n")
  cat("psi50S = ", psi50S, "\n")
  cat("psi12S = ", psi50S, "\n")
  
  dpx = df_w_vol$var
  dpy = df_w_vol$dpsi
  f1 = splinefun(dpy~dpx)
  dpx1 = seq(min(dpx), max(dpx), length.out=100)
  dp88S = f1(psi88S)
  dp50S = f1(psi50S)
  dp12S = f1(psi12S)
  cat("psiL88S = ", psi88S-dp88S, "\n")
  cat("psiL50S = ", psi50S-dp50S, "\n")
  cat("psiL12S = ", psi12S-dp12S, "\n")
  
  cat(psi88S, "\t", psi50S, "\t", psi12S, "\t", psi88S-dp88S, "\t", psi50S-dp50S, psi12S-dp12S, "\n")
  
  subdata = data %>% filter(LWP > -6)
  
  p1 <- df_w_vol %>% 
    ggplot() +
    geom_line(aes(x = var, y = vcmax), col="green3", size=1) +
    # geom_point(data=filter(dat, Species==species), aes(x=LWP, y=Vcmax))+
    # geom_vline(xintercept = par_plant_std$psi50, col="orange") +
    expand_limits(y=0)+
    xlab(varname)
  if (analytical) p1 = p1 + geom_line(aes(x = var, y = out_analytical_vcmax ), col="grey", size=1) 
  
  p2 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = dpsi), col="blue", size=1)+
    # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
    expand_limits(y=0)+
    xlab(varname)
  if (!is.null(dpsi_data)){
    cat("Adding points...\n")
    p2 <- p2 + 
      geom_point(data=dpsi_data, aes(x=SWP, y=Dpsi))
  }
  
  p3 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = gs), col="cyan2", size=1)+
    geom_point(data=subdata, aes(x=LWP, y=gC))+
    geom_vline(xintercept = psi88S, col="grey")+
    expand_limits(y=0)+
    xlab(varname)
  if (analytical) p3 = p3 + geom_line(aes(x = var, y = out_analytical_gs), col="grey", size=1)
  
  p4 <- df_w_vol %>%
    # mutate(chi = ci/out_hydraulics_ca) %>% 
    ggplot() +
    geom_line(aes(x = var, y = chi), col="magenta", size=1)+
    geom_point(data=subdata, aes(x=LWP, y=1-A/gC/ca))+
    # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
    expand_limits(y=0)+
    xlab(varname)
  if (analytical) p4 = p4 + geom_line(aes(x = var, y = out_analytical_chi), col="grey", size=1)

  p5 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = jmax), col="goldenrod1", size=1) +
    # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
    # geom_point(data=filter(dat, Species==species), aes(x=LWP, y=Jmax))+
    expand_limits(y=0)+
    xlab(varname)
  
  
  p6 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = a), col="green4", size=1) +
    geom_point(data=subdata, aes(x=LWP, y=A))+
    # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
    expand_limits(y=0)+
    xlab(varname)
  if (analytical) p6 = p6 + geom_line(aes(x = var, y = out_analytical_gpp), col="grey",size=1) 
  
  grid.arrange(p3,p6,p1,p2,p5,p4, ncol=2)
}



error_fun = function(x, data, dpsi_file="", plot=F, scale = 1, dpsi_calib=T, inst=F, k=7){
  x = x*scale
  # cat("start: ",x, "\n")
  
  data$Ciest = data$ca-data$A/data$gC
  
  par_plant_now = list(
    conductivity= x[1]*1e-16, #.7, 
    psi50 = x[2], #-.5, 
    b= 1 #x[3] #1.2
  )
  
  par_cost_now = list(
    alpha  = x[4], #0.13,       # cost of Jmax
    gamma = x[5] #.1,         # cost of hydraulic repair
  )

  if (dpsi_file != ""){
    dpsi_data = dpsi_df %>% filter(Species == dpsi_file) #read.csv(paste0("drying_experiments_meta/",dpsi_file,".csv"))
  }
  else{
    dpsi_data=NULL
  }
  
  if (inst == F){
    cat("Acc resp\n")
    lwp = seq(min(data$LWP), 0, length.out=20)
  dat1 = tibble(var = lwp) %>%
    mutate(p = purrr::map(var, ~pmodel_calibrate_analytical(tc = mean(data$T), ppfd = mean(data$Iabs_growth), vpd = mean(data$D*101325), co2 = mean(data$ca), elv = 0, fapar = .99, kphio = 0.087, psi_soil = ., rdark = 0.02, par_plant=par_plant_now, par_cost = par_cost_now, opt_hypothesis = "PM")) ) %>% unnest_wider(p) %>% unnest_wider(out_hydraulics) %>% unnest_wider(out_analytical, names_sep = "_")
  } else{
  cat("inst resp\n")
  ndays = mean(data$Drydown.days)
  psi_min = -5 #min(data$LWP)
  psi_max = 0 #max(data$LWP)
  # cat(ndays,"\n")

  lwp = seq(psi_min,0, length.out=10)
  day = ndays * (lwp-psi_max)/(psi_min-psi_max)

  lwp_day = function(day_num){
    psi_max + day_num/ndays * (psi_min-psi_max)
  }

  # k = 7
  lwp_week = rollmean(x = lwp_day(c(max(day):0, rep(0,k-1))), k = k, align = "right")

  spl = splinefun(x = max(day):0, y=lwp_week)

  # plot(x=day, y=lwp_day(day))
  # points(x=day, y=spl(day), col="red", type="o")
  # points(x=max(day):0, y=lwp_week, col="blue", type="l")

  dat_acc = tibble(var = spl(day)) %>% mutate(pmod = map(var, ~rphydro_analytical(tc = mean(data$T), ppfd = mean(data$Iabs_growth), vpd = mean(data$D*101325), co2 = mean(data$ca), elv = 0, fapar = .99, kphio = 0.087, psi_soil = ., rdark = 0.02, par_plant=par_plant_now, par_cost = par_cost_now, opt_hypothesis = "PM"))) %>% unnest_wider(pmod)
  dat1 = tibble(var = lwp, jmax_a=dat_acc$jmax, vcmax_a=dat_acc$vcmax) %>%
    mutate(p = purrr::pmap(list(var, jmax_a, vcmax_a), ~pmodel_calibrate_inst(tc = mean(data$T), ppfd = mean(data$Iabs_used), vpd = mean(data$D*101325), co2 = mean(data$ca), elv = 0, fapar = .99, kphio = 0.087, psi_soil = ..1, rdark = 0.02, par_plant=par_plant_now, par_cost = par_cost_now, jmax = ..2, vcmax = ..3, opt_hypothesis = "PM")) ) %>% unnest_wider(p) %>% unnest_wider(out_hydraulics) %>% unnest_wider(out_analytical, names_sep = "_")
  }

  # if(plot==T) dat_acc %>% plot_all(varname = "psi_soil", species=species, data = data, dpsi_data=dpsi_data, analytical = F)
  if(plot==T) dat1 %>% plot_all(varname = "psi_soil", species=species, data = data, dpsi_data=dpsi_data)
  
  a_spl = splinefun(x = lwp, y=dat1$a)
  g_spl = splinefun(x = lwp, y=dat1$gs)
  c_spl = splinefun(x = lwp, y=dat1$chi)

  w = rep(1, length(data$LWP))
  # w[data$LWP < -1.75] = 0
  y1 = mean((a_spl(data$LWP) - data$A)^2*w)/mean(data$A)^2
  y2 = mean((g_spl(data$LWP) - data$gC)^2*w)/mean(data$gC)^2 #*4000

  # w = exp(-data$LWP/par_plant_now$psi50)
  y4 = mean((c_spl(data$LWP) - data$Ciest/data$ca)^2*w)/mean(data$Ciest/data$ca)^2 #*800
  cat("c_spl:", c_spl(data$LWP), "\n")
  
  if (dpsi_calib == T){
    d_spl = splinefun(lwp, y=dat1$dpsi)
    w = 1 #0.1 + exp(dpsi_data$SWP)
    y3 = mean((d_spl(dpsi_data$SWP) - dpsi_data$Dpsi)^2*w)/mean(dpsi_data$Dpsi)^2 #*40
    cat("d_spl:", d_spl(dpsi_data$SWP), "\n")
  }else{
    # w=1
    # p50 = -1.5 #par_plant_now$psi50
    # y3 = mean((dat1$dpsi[lwp>p50] - -p50)^2*w)/mean(dat1$dpsi[lwp>p50])^2 #*20
    y3=0
  }
  y = (y1+y2+y3+y4)*100
  cat(x, "|", y1, "/ ", y2, " / ", y3, " / ", y4, " / ", y, "\n")
  cat(x, "|", y, "\n")
  

  x0 <<- x
  y
}






species = "Allocasuarina luehmannii"
data1=filter(dat, Species==species)

# data1 = data1 %>% arrange(desc(LWP)) %>% mutate(day = abs(LWP)/(max(abs(LWP)))*Drydown.days )
# data1 = data1 %>% filter(LWP > -1.5)

x0 = as.numeric((par_data %>% filter(Species == species) )[2:6])
# x0 = c(15.09, -0.35, 1.4, 0.028, 4.9977)

error_fun(x0, data1, plot=T, dpsi_file=species, dpsi_calib = T, inst=T)


# x0= c(0.1202499, -1.655592, 3.508305, 0.1070148, 0.08062944)
# x0 = c(0.4051118, -1.135582, 16.36669, 0.1050266, 1.306852)
# error_fun(x0, data1 %>% mutate(Drydown.days=38), plot=T, dpsi_file=species, dpsi_calib = T, inst=T)

optimr::optimr(par = x0/abs(x0),
               fn = error_fun,
               data=data1,
               scale=abs(x0),
               dpsi_file = species,  # Set to species name or "" if dpsi data is not available
               dpsi_calib = T,       # Set to F if dpsi data is not available
               inst=T,
               control = list(par_scale = 1, maxit=500)
               )
# error_fun(x0, data1, plot=T, dpsi_file="")

error_fun(x0, data1, plot=T, dpsi_file=species)


### Cross Validation
species_master = unique(dat$Species) %>% str_subset(pattern = "Mediterranean", negate = T) %>% str_subset("Glycine max2", negate=T)
res_master = data.frame()

K=5

species = species_master[1]
species_short = paste(substr(strsplit(species, " ")[[1]][-(5:8)], 0,2), collapse = ".")
data1=filter(dat, Species==species)

data1 = data1[sample(1:nrow(data1), replace = F), ]
N = ceiling(nrow(data1)/K)
if (nrow(data1) < 2*K-1) N=1
start_ids = seq(1, nrow(data1), N)
start_ids = c(start_ids, nrow(data1)+1)
cat("splits = ", start_ids, "\n")


ids_list = list()
for (i in 1:K){
  ids_list[[i]] = start_ids[i]:(start_ids[i+1]-1)
}

for (iter in 1:5){
  ids_test = unlist(ids_list[iter])
  ids_train = unlist(ids_list[-iter])
  
  dat_train = data1[ids_train,]
  dat_test  = data1[ids_test,]
  
  x0 = as.numeric((par_data %>% filter(Species == species) )[2:6])
  names(x0) = c("K.scalar", "P50", "b", "alpha" ,"yxp50")
  error_fun(x0, dat_train, plot=F, dpsi_file="")
  
  out = optimr::optimr(par = x0,
                 fn = error_fun,
                 data=dat_train,
                 # dpsi_file = species
                 control = list(par_scale = scale, maxit=100)
  )
  
  res = data.frame(as.list(out$par))
  res$height = mean(dat_train$Ht_used)
  res$species = species
  res$iter = iter
  res$err_train = error_fun(out$par, dat_train, plot=F, dpsi_file="")
  res$err_test  = error_fun(out$par, dat_test, plot=F, dpsi_file="")
  
  res_master = rbind(res_master, res)
}

write.csv(res_master, file = paste0("cross_validation/cross_validation_", species_short,".csv"))


# error_fun(out$par, dat_train, plot=T, dpsi_file="")
# error_fun(out$par, dat_test, plot=T, dpsi_file="")

# #### Analysis ####
# 
# library(MASS)
# library(tidyverse)
# d = read.csv("C:/Users/Jaideep/Documents/codes/rpmodel/vignettes_add/drying_experiments_meta/fitted_params_brd_all_3_p50pmin.csv")
# # d = d %>% filter(!(Species %in% c("Helianthus annuus", "Glycine max")))
# d = d %>% mutate(gamma = yxp50/P50^2)
# d = d %>% mutate(P88 = (log(0.12)/log(0.5))^(1/b)*P50 )
# d = d %>% mutate(K.scalar = K.scalar/H)
# d = d %>% dplyr::select(-Species, -Ref.1, -Ref2, -Ref, -b..Cavit., -b..SC., -Pclose, -Slope, -P12, -P50..stem., -PgS88, -PgS50, -pg12, -H)
# d = d %>% dplyr::select(-SLA_nopetiole, -SLA_wpetiole, -SLA1, -A.G)
# d = d %>% dplyr::select(-Pmin, -dpsi_int, -P88)
# d = d %>% dplyr::select(-pgL50, -Pgs90, -b)
# d = d %>% mutate(SM_skelton = -P50 + pgL12)
# # d = d %>% mutate(SM_choat = P50 - P88)
# # d1 = sapply(d, FUN = function(x){ifelse(x<0, -x, x)}, simplify = T)
# # d = d %>% mutate(logb = log(b)) %>% dplyr::select(-b)
# # d = d %>% mutate(logY = log(gamma)) %>% dplyr::select(-gamma)
# # d1 = d1 %>% mutate(logb = log(b)) %>% dplyr::select(-b)
# # d = d %>% mutate(logP50 = log(-P50)) %>% dplyr::select(-P50)
# # d = d %>% mutate(logK = log(K.scalar)) %>% dplyr::select(-K.scalar)
# 
# panel_func = function(x,y,...){
#   points(y~x,...);
#   mod = rlm(y~x, weights = rep(1, length(x)) );
#   if (abs(summary(mod)$coefficients[2,3]) > 1.96) {
#     abline(lm(y~x), col="black")
#   }
#   else {
#     if (abs(summary(mod)$coefficients[2,3]) > .9){
#       abline(lm(y~x), col="grey65")
#     }
#   }
# }
# 
# d %>% pairs(panel= panel_func, pch=20, cex=1.2, col="cyan4") #
# 
# d %>% dplyr::select(P88, pgL88, Ptlp) %>% pairs(panel= panel_func)
# 
# plot(d$pgL12~d$P50)
# abline(0,1)
# 
# par(mfrow=c(2,3), mar=c(4,4,1,1), cex.lab=1.2)
# 
# with(d1, plot(y=alpha, x=gamma, xlab="log(gamma)", ylab="log(alpha)"))
# abline(rlm((d1$alpha)~(d1$gamma), method="MM"), col="cyan3")
# 
# with(d1, plot(log(K.scalar)~log(P50)))
# abline(rlm(log(d1$K.scalar)~log(d1$P50), method="MM"), col="cyan3")
# 
# with(d1, plot((P50)~(P88)))
# abline(rlm((P50)~log(P88), method="MM"), col="cyan3")
# 
# with(d1, plot((P50)~(pgL12)))
# abline(rlm((d1$P50)~(d1$pgL12), method="MM"), col="cyan3")
# 
# with(d1, plot(log(K.scalar)~(alpha), xlab="log(alpha)"))
# abline(rlm(log(d1$K.scalar)~(d1$alpha), method="MM"), col="cyan3")
# 
# with(d1, plot((pgL50)~(X.tlp)))
# abline(rlm((d1$pgL50)~(d1$X.tlp), method="MM"), col="cyan3")
# 
