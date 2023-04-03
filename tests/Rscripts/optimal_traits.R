library(tidyverse)
library(PlantFATE)

# Ambient Env
zp_amb = c(19.22669, 15.71661, 3.15710, 0.00000)
co_amb = c(1.0000000, 0.4712540, 0.2218492, 0.0711068)

# zp_amb = c(as.numeric(colMeans(Zp[901:1000,-1], na.rm=T))[1:3], 0)
# co_amb = c(as.numeric(colMeans(co[901:1000,-1], na.rm=T))[1:4])

# eCO2 env
zp_ele = c(20.695389, 18.106550, 14.087510, 2.206985, 0.000000) 
co_ele = c(1.000000, 4.712543e-01, 2.220795e-01, 1.032055e-01, 2.763073e-02)

# zp_ele = c(as.numeric(colMeans(Zp[1901:2000,-1], na.rm=T))[1:4], 0)
# co_ele = c(as.numeric(colMeans(co[1901:2000,-1], na.rm=T))[1:5])


setwd("~/codes/Plant-FATE/")

fitness = function(par, modelpar, zeta= 0.2, zp, co, co2, pfile="tests/params/p.ini", ewd0, e_alpha, e_gamma){
  lho = new(LifeHistoryOptimizer, pfile)
  lho$set_metFile("tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv")
  # lho$traits0$lma = lma

  lho$traits0$lma          = 0.12 # par[1]
  lho$traits0$wood_density = par[1]
  lho$traits0$hmat         = par[2]
  lho$traits0$p50_xylem    = par[3]
  lho$traits0$zeta = zeta
  
  lho$par0$eWD       = ewd0
  lho$par0$eWD_alpha = e_alpha
  lho$par0$eWD_gamma = e_gamma
  
  lho$env$clim$co2 = co2
  lho$env$z_star = zp
  lho$env$canopy_openness = co
  
  lho$init()
  # lho$set_traits(c(traits$lma, traits$wd))
  # lho$printPlant()
  fitness1 = lho$calcFitness()
  cat(par, " | ", fitness1, "\n")
  -fitness1
}

# traits_obs = read.csv(file = "tests/data/Amz_trait_filled_HD.csv") %>% arrange(desc(Total.BasalArea_2017.cm2.))
traits_used = read.csv(file = "tests/data/Traits_random_HD2.csv") %>% arrange(desc(Total.BasalArea_2017.cm2.))

par0 = list(ewd0 = -0.5, e_alpha = -1.14, e_gamma = -1.8)

traits = c(600, 21, -2.5)

fitness(par = traits, zp=zp_amb, co=co_amb, co2 = 368, ewd0=-1, e_alpha = -1.14, e_gamma = -1.8)

calc_opt_traits = function(traits, zp, co, co2, ewd0, e_alpha, e_gamma){
  opt = optim(par = traits, fn = fitness, zp=zp_ele, co=co_ele, co2 = 568, method = "L-BFGS-B", control = list(parscale = c(600, 25, 2.5), maxit = 50), lower=c(250, 20, -6), upper = c(950, 25, -0.5), ewd0=ewd0, e_alpha = e_alpha, e_gamma = e_gamma)
  c(opt$par, opt$value)
}

df = 
  list(e_alpha = -1.1493, # seq(-2,-1, length.out = 5),
       e_gamma = -1.8392, # seq(-2,-1, length.out = 5),
       ewd0 = -1.0 # seq(-2,-0.5, length.out = 4),
  ) %>% 
  cross_df() %>% 
  mutate(fitness_a = pmap(.l = cur_data(), .f = calc_opt_traits, traits = traits, zp=zp_amb, co=co_amb, co2=368)) %>% 
  unnest_wider(fitness_a, names_sep = "_opt")


par_df = 
  list(e_alpha = seq(-2,-1, length.out = 5),
       e_gamma = seq(-2,-1, length.out = 5),
       ewd0    = seq(-2,-0.5, length.out = 4)
  ) %>% 
  cross_df()

df1 = par_df %>% 
  mutate(fitness_a = pmap(.l = cur_data(), .f = calc_opt_traits, traits = traits, zp=zp_amb, co=co_amb, co2=368)) %>% 
  unnest_wider(fitness_a, names_sep = "_opt")

df2 = par_df %>% 
  mutate(fitness_e = pmap(.l = cur_data(), .f = calc_opt_traits, traits = traits, zp=zp_ele, co=co_ele, co2=614)) %>% 
  unnest_wider(fitness_e, names_sep = "_opt")

df_all = left_join(df1, df2)
