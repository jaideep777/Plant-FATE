rm(list=ls())

library(rphydro)
library(tidyverse)
library(reshape2)

kphio = 0.087;        # quantum yield efficiency
ppfd = 300;          # umol/m2/s
vpd  = 1000;         # Pa
co2  = 400;          # ppm
elv  = 0;            # m.a.s.l.
fapar = 0.7;         # fraction
rdark = 0.015;
tc = 25;
vwind = 3;
netrad = ppfd/2
pa = calc_patm(elv)

par_cost = list(alpha=0.1, gamma=1)

options = list(gs_method = "GS_IGF", 
               et_method = "ET_DIFFUSION",
               ftemp_vj_method = "FV_kumarathunge19",
               ftemp_rd_method = "FR_heskel16",
               ftemp_br_method = "FB_atkin15",
               scale_alpha = F)

vpd_acc = 100

## calculation and plot with isohydric species (hydraulic acclimation) 
par_plant_acclim = list(conductivity=3e-17, 
                        psi50=-2, 
                        b=2)

acc1 = rphydro_analytical(tc, tc, ppfd, netrad, vpd_acc, co2, pa, fapar, kphio, 0, rdark, vwind, par_plant_acclim, par_cost, options)

dat1 = list(vpd = exp(seq(log(100),log(5000),length.out=50)),
           psi_soil = c(0, -2.5)) %>%
  cross_df() %>% 
  mutate(dat = purrr::map2(.x=psi_soil, .y=vpd, .f = ~rphydro_instantaneous_analytical(acc1$vcmax25, acc1$jmax25, tc, tc, ppfd, netrad, .y, co2, pa, fapar, kphio, .x, rdark, vwind, par_plant_acclim, par_cost, options))) %>% 
  unnest_wider(dat)


## calculation and plot with anisohydric species (no hydraulic acclimation) 
par_plant_no_acclim = list(conductivity=3e-17, 
                           psi50=-5, 
                           b=2)

acc2 = rphydro_analytical(tc, tc, ppfd, netrad, vpd_acc, co2, pa, fapar, kphio, 0, rdark, vwind, par_plant_no_acclim, par_cost, options)

dat2 = list(vpd = exp(seq(log(100),log(5000),length.out=50)),
           psi_soil = c(0, -2.5)) %>%
  cross_df() %>% 
  mutate(dat = purrr::map2(.x=psi_soil, .y=vpd, .f = ~rphydro_instantaneous_analytical(acc2$vcmax25, acc2$jmax25, tc, tc, ppfd, netrad, .y, co2, pa, fapar, kphio, .x, rdark, vwind, par_plant_no_acclim, par_cost, options))) %>% 
  unnest_wider(dat)

## Make the desired plot

dat1 %>% 
  select(vpd, psi_soil, gs, e) %>% 
  mutate(type = "With acclimation") %>% 
  rbind(dat2 %>% 
          select(vpd, psi_soil, gs, e) %>% 
          mutate(type = "No acclimation")) %>% 
  melt(c("vpd", "psi_soil", "type")) %>% 
  ggplot(aes(x=vpd, y=value, col=factor(psi_soil), group=factor(psi_soil))) + 
  facet_wrap(facets = c("variable", "type"), scales = "free", nrow = 2)+
  geom_line()
