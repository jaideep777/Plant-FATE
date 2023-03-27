library(PlantFATE)

# Ambient Env
zp_amb = c(19.22669, 15.71661, 3.15710, 0.00000)
co_amb = c(1.0000000, 0.4712540, 0.2218492, 0.0711068)

zp_amb = c(as.numeric(colMeans(Zp[901:1000,-1], na.rm=T))[1:3], 0)
co_amb = c(as.numeric(colMeans(co[901:1000,-1], na.rm=T))[1:4], 0)

# eCO2 env
zp_ele = c(20.695389, 18.106550, 14.087510, 2.206985, 0.000000) 
co_ele = c(1.000000, 4.712543e-01, 2.220795e-01, 1.032055e-01, 2.763073e-02)

setwd("~/codes/Plant-FATE/")

fitness_wd = function(wd, zp, co, co2, zeta=0.2, lma, hmat, p50, pfile="tests/params/p.ini"){
  lho = new(LifeHistoryOptimizer, pfile)
  lho$set_metFile("tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv")
  # lho$traits0$lma = lma
  lho$traits0$hmat = hmat
  lho$traits0$p50_xylem = p50
  lho$traits0$zeta = zeta
  lho$env$clim$co2 = co2
  lho$env$z_star = zp
  lho$env$canopy_openness = co
  lho$init()
  lho$set_traits(c(lma, wd))
  # lho$printPlant()
  lho$calcFitness()
}

traits_obs = read.csv(file = "tests/data/Amz_trait_filled_HD.csv") %>% arrange(desc(Total.BasalArea_2017.cm2.))
traits_used = read.csv(file = "tests/data/Traits_random_HD.csv") %>% arrange(desc(Total.BasalArea_2017.cm2.))

wd = seq(350, 900, length.out=20)

N = 100
fitness_landscape = numeric(N)
for (i in 1:N){
  cat(i, "\n")
  fitness = fitness_wd(wd = traits$WD[i], zp=zp_amb, co=co_amb, co2=414, lma=traits$LMA[i], hmat=traits$HMAT[i], p50 = traits$P50X[i])
  fitness_landscape[i] = fitness
}

p1 = traits %>% filter(YEAR == 1000) %>% 
  ggplot(aes(y=HMAT, x=WD, col=exp(exp(fitness_landscape)))) + 
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
  ggplot(aes(y=P50X, x=WD, col=exp(exp(fitness_landscape)))) + 
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

fitness_landscape_ele_nofeedback = numeric(N)
for (i in 1:N){
  cat(i, "\n")
  fitness = fitness_wd(wd = traits$WD[i], zp=zp_amb, co=co_ele, co2=414, lma=traits$LMA[i], hmat=traits$HMAT[i], p50 = traits$P50X[i])
  fitness_landscape_ele_nofeedback[i] = fitness
}

fitness_landscape_ele_wfeedback = numeric(N)
for (i in 1:N){
  cat(i, "\n")
  fitness = fitness_wd(wd = traits$WD[i], zp=zp_ele, co=co_ele, co2=414, lma=traits$LMA[i], hmat=traits$HMAT[i], p50 = traits$P50X[i])
  fitness_landscape_ele_wfeedback[i] = fitness
}

par(mfrow = c(1,1))
boxplot(cbind(fitness_landscape, fitness_landscape_ele_nofeedback, fitness_landscape_ele_wfeedback))

wd_opt = function(wd, fitness){
  sum(wd*fitness)/sum(fitness)
}

par(mfrow=c(1,1))
barplot(c(wd_opt(traits$WD[1:N],exp(exp(fitness_landscape))),
          wd_opt(traits$WD[1:N],exp(exp(fitness_landscape_ele_nofeedback))),
          wd_opt(traits$WD[1:N],exp(exp(fitness_landscape_ele_wfeedback)))
          ), ylim=c(500,1000), xpd=F)



# 
# N = 100
# dat_amb = matrix(nrow=N, ncol=length(wd))
# for (i in 1:N){
#   cat(i, "\n")
#   fitness = wd %>% purrr::map_dbl(fitness_wd, zp=zp_amb, co=co_amb, co2=414, lma=traits_used$Leaf.LMA..g.m2.[i]/1000, hmat=traits_used$Height_Max.m.[i], p50 = traits_used$P50..Mpa.)
#   dat_amb[i,] = fitness
# }
# matplot(y=t(dat_amb), x=wd, lty=1, type="l", col=scales::alpha(scales::colour_ramp(c("blue", "red"))(seq(0,1,length.out=N)), alpha = 0.3))
# dat_amb_opt = wd[apply(dat_amb[1:10,], MARGIN=1, FUN = function(x){which(x==max(x))})]


# dat_ele_base = wd %>% purrr::map_dbl(fitness_wd, zp=zp_amb, co=co_amb, co2=614, lma=)
# dat_ele = wd %>% purrr::map_dbl(fitness_wd, zp=zp_ele, co=co_ele, co2=614)
