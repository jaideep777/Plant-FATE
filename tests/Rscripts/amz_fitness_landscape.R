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

fitness_wd = function(wd, zp, co, co2, zeta=0.2, lma, hmat, p50, pfile="tests/params/p.ini", ewd0, e_alpha, e_gamma){
  lho = new(LifeHistoryOptimizer, pfile)
  lho$set_metFile("tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv")
  # lho$traits0$lma = lma
  lho$traits0$hmat = hmat
  lho$traits0$p50_xylem = p50
  lho$traits0$zeta = zeta
  
  lho$par0$eWD = ewd0
  lho$par0$eWD_alpha = e_alpha
  lho$par0$eWD_gamma = e_gamma

  lho$env$clim$co2 = co2
  lho$env$z_star = zp
  lho$env$canopy_openness = co
  
  lho$init()
  lho$set_traits(c(lma, wd))
  # lho$printPlant()
  lho$calcFitness()
}

# traits_obs = read.csv(file = "tests/data/Amz_trait_filled_HD.csv") %>% arrange(desc(Total.BasalArea_2017.cm2.))
traits_used = read.csv(file = "tests/data/Traits_random_HD2.csv") %>% arrange(desc(Total.BasalArea_2017.cm2.))

fitness_spp = function(i, zp, co, co2, ewd0, e_alpha, e_gamma, ...){
  cat("ewd0 =", ewd0, ", ea =", e_alpha, ", eg =", e_gamma, ", spp = ", i, "\n")
  fitness_wd(wd = traits_used$meanWoodDensity..g.cm3.[i]*1000, 
             zp =zp, co=co, co2=co2, 
             lma=traits_used$Leaf.LMA..g.m2.[i]/1000, 
             hmat=traits_used$Height_Max.m.[i], 
             p50 = traits_used$P50..Mpa.[i], 
             ewd0 = ewd0, e_alpha = e_alpha, e_gamma = e_gamma)
}

0.6251576
# wd = seq(350, 900, length.out=20)

N = 100

df = 
  list(e_alpha = -1.1493, # seq(-2,-1, length.out = 5),
       e_gamma = -1.8392, # seq(-2,-1, length.out = 5),
       ewd0 = -1.0, # seq(-2,-0.5, length.out = 4),
       i = 1:N
  ) %>% 
  cross_df() %>% 
  mutate(fitness_aCaE = pmap(.l = cur_data(), .f = fitness_spp, zp=zp_amb, co=co_amb, co2=368) %>% unlist(), 
         fitness_eCaE = pmap(.l = cur_data(), .f = fitness_spp, zp=zp_amb, co=co_amb, co2=614) %>% unlist(),
         fitness_eCeE = pmap(.l = cur_data(), .f = fitness_spp, zp=zp_ele, co=co_ele, co2=614) %>% unlist())

fitness_landscape = df$fitness_aCaE

# fitness_landscape = numeric(N)
# for (i in 1:N){
#   cat(i, "\n")
#   # fitness = fitness_wd(wd = traits_used$meanWoodDensity..g.cm3.[i]*1000, zp=zp_amb, co=co_amb, co2=368, lma=traits_used$Leaf.LMA..g.m2.[i]/1000, hmat=traits_used$Height_Max.m.[i], p50 = traits_used$P50..Mpa.[i], ewd0 = -1.0, e_alpha = -1.1493, e_gamma = -1.8392)
#   fitness = fitness_spp(i=i, zp=zp_amb, co=co_amb, co2=368, ewd0 = -1.0, e_alpha = -1.1493, e_gamma = -1.8392)
#   fitness_landscape[i] = fitness
# }

p1 = traits %>% filter(YEAR == 1000) %>% 
  ggplot(aes(y=HMAT, x=WD, col=((fitness_landscape^10)))) + 
  theme_classic(base_size = 14)+
  geom_point()+
  labs(col="Fitness")+
  scale_colour_viridis_c(direction = -1)+
  ggtitle("Expected")

p2 = dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR == 2000) %>% 
  # filter(RES==T) %>% 
  ggplot(aes(y=HMAT, x=WD))+
  theme_classic(base_size = 14)+
  geom_point(aes(col=BA*1e4), alpha=0.7)+
  scale_color_viridis_c(direction = -1)+
  # scale_size("size_RES", range = c(0, 1.5), guide = F)+
  labs(col="BA")+
  ggtitle("Emergent")


q1 = traits %>% filter(YEAR == 1000) %>% 
  ggplot(aes(y=P50X, x=WD, col=((fitness_landscape^10)))) + 
  theme_classic(base_size = 14)+
  geom_point()+
  labs(col="Fitness")+
  scale_colour_viridis_c(direction = -1)+
  ggtitle("Expected")

q2 = dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR == 2000) %>% 
  # filter(RES==T) %>% 
  ggplot(aes(y=P50X, x=WD))+
  theme_classic(base_size = 14)+
  geom_point(aes(col=BA*1e4), alpha=0.7)+
  scale_color_viridis_c(direction = -1)+
  # scale_size("size_RES", range = c(0, 1.5), guide = F)+
  labs(col="BA")+
  ggtitle("Emergent")


print(cowplot::plot_grid(p1,p2,q1,q2, align="hv"))


par(mfrow = c(1,1))
boxplot(cbind(df$fitness_aCaE, df$fitness_eCaE, df$fitness_eCeE))

wd_opt = function(wd, fitness){
  sum(wd*fitness)/sum(fitness)
}

wd_aCaE   = wd_opt(traits_used$meanWoodDensity..g.cm3.[1:N]*1000, ((df$fitness_aCaE^10)))
wd_eCaE   = wd_opt(traits_used$meanWoodDensity..g.cm3.[1:N]*1000, ((df$fitness_eCaE^10)))
wd_eCeE   = wd_opt(traits_used$meanWoodDensity..g.cm3.[1:N]*1000, ((df$fitness_eCeE^10)))

par(mfrow=c(1,1))
barplot(c(wd_aCaE,
          wd_eCaE,
          wd_eCeE
          ), ylim=c(700,900), xpd=F)



df_scan1 = 
  list(i = 1:N,
       e_alpha = seq(-2,-1, length.out = 3),
       e_gamma = seq(-2,-1, length.out = 3),
       ewd0 = seq(-2,-0.5, length.out = 4)
  ) %>% 
  cross_df()

df1 = df_scan1 %>% mutate(fitness_aCaE = pmap(.l = cur_data(), .f = fitness_spp, zp=zp_amb, co=co_amb, co2=368) %>% unlist())
df2 = df_scan1 %>% mutate(fitness_eCaE = pmap(.l = cur_data(), .f = fitness_spp, zp=zp_amb, co=co_amb, co2=614) %>% unlist())
df3 = df_scan1 %>% mutate(fitness_eCeE = pmap(.l = cur_data(), .f = fitness_spp, zp=zp_ele, co=co_ele, co2=614) %>% unlist())

df_scan = df1 %>% left_join(df2) %>% left_join(df3)
         
df_scan %>% write.csv("~/codes/Plant-FATE/fitness_scan.csv")

df_scan = read.csv("~/codes/Plant-FATE/fitness_scan.csv")


wd_diff = traits_used %>% 
  select(Species, meanWoodDensity..g.cm3.) %>% 
  mutate(i = 1:N) %>% 
  right_join(df_scan) %>% 
  group_by(e_alpha, e_gamma, ewd0) %>% summarize(wd_aCaE = wd_opt(wd = meanWoodDensity..g.cm3.*1e3, fitness = fitness_aCaE^10),
                                                 wd_eCaE = wd_opt(wd = meanWoodDensity..g.cm3.*1e3, fitness = fitness_eCaE^10),
                                                 wd_eCeE = wd_opt(wd = meanWoodDensity..g.cm3.*1e3, fitness = fitness_eCeE^10)) %>% 
  mutate(delta_wd_eCaE = wd_eCaE - wd_aCaE,
         delta_wd_eCeE = wd_eCeE - wd_aCaE)


wd_diff %>% ggplot(aes(x=e_alpha, y=e_gamma)) + 
  geom_raster(aes(fill=delta_wd_eCeE)) + 
  facet_wrap(~ewd0)+
  scale_fill_viridis_c(direction = -1, limits=c(-220,0))

wd_diff %>% ggplot(aes(x=e_alpha, y=e_gamma)) + 
  geom_raster(aes(fill=delta_wd_eCaE)) + 
  facet_wrap(~ewd0)+
  scale_fill_viridis_c(limits = c(0,270))

wd_diff %>% ggplot(aes(x=e_alpha, y=e_gamma)) + 
  geom_raster(aes(fill=wd_aCaE)) + 
  facet_wrap(~ewd0)+
  scale_fill_viridis_c(limits = c(400, 1000))


# 
# N = 100
# dat_amb = matrix(nrow=N, ncol=length(wd))
# for (i in 1:N){
#   cat(i, "\n")
#   fitness = wd %>% purrr::map_dbl(fitness_wd, zp=zp_amb, co=co_amb, co2=414, lma=traits_used$Leaf.LMA..g.m2.[i]/1000/1000, hmat=traits_used$Height_Max.m.[i], p50 = traits_used$P50..Mpa.)
#   dat_amb[i,] = fitness
# }
# matplot(y=t(dat_amb), x=wd, lty=1, type="l", col=scales::alpha(scales::colour_ramp(c("blue", "red"))(seq(0,1,length.out=N)), alpha = 0.3))
# dat_amb_opt = wd[apply(dat_amb[1:10,], MARGIN=1, FUN = function(x){which(x==max(x))})]


# dat_ele_base = wd %>% purrr::map_dbl(fitness_wd, zp=zp_amb, co=co_amb, co2=614, lma=)
# dat_ele = wd %>% purrr::map_dbl(fitness_wd, zp=zp_ele, co=co_ele, co2=614)
