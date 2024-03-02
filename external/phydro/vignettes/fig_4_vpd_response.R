rm(list=ls())
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(scales)
library(rphydro)

setwd("~/codes/phydro/vignettes/")
# P-model as per Wang et al 2017
# This needs PPFD in mol/m2/day

source("utils.R")

lm_eqn <- function(x,y, s="slope = "){
  m <- lm(y ~ x);
  # eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
  #      list(a = format(unname(coef(m)[1]), digits = 2),
  #           b = format(unname(coef(m)[2]), digits = 2),
  #          r2 = format(summary(m)$r.squared, digits = 3)))
  # as.character(as.expression(eq))
  cat(m$coefficients[2], "\n")
  D_exp <<- c(D_exp, m$coefficients[2])
  sprintf("%s%.2f", s,m$coefficients[2])
}

plot_slope = function(dat){
  
  dat = dat %>%
    mutate(chi = out_hydraulics_ci/out_analytical_ca) %>%
    mutate(logit_chi = log(chi/(1-chi))) %>%
    mutate(logit_chi_ana = log(out_analytical_chi/(1-out_analytical_chi)))
  
  dat %>% ggplot() +
    theme_classic()+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16),
          plot.tag.position = "topleft") +
    geom_smooth(aes(x=log(var), y=logit_chi), method="lm", se=F, color="magenta4")+
    geom_point(aes(x=log(var), y=logit_chi), size=1, color="magenta")+
    geom_text(x = 7, y = max(dat$logit_chi)-1., label = lm_eqn(log(dat$var), dat$logit_chi, "This study: "), parse = F, color="magenta4")+
    
    # geom_smooth(aes(x=log(var), y=logit_chi_ana), method="lm", se=F, color="grey40")+
    # geom_point(aes(x=log(var), y=logit_chi_ana), size=1, color="grey")+
    # geom_text(x = 6, y = max(dat$logit_chi)-.5, label = lm_eqn(log(dat$var), dat$logit_chi_ana, "Wang et al. (2017): "), parse = F, color="grey40")+
    ylab(expression("logit("*chi*")"))+
    xlab("log(D)")
  
}



dat = read.csv("data/drying_experiments_meta-analysis_Joshi_et_al_2022.csv")

cpar = read.csv("data/fitted_params_Joshi_et_al_2022.csv") 

spp = cpar$Species

D_exp = numeric(0)


for (species in spp){
  
  
  data=filter(dat, Species==species)
  x = as.numeric(cpar[which(cpar$Species == species), 2:6])
  
  par_plant_now = list(
    conductivity = x[1]*1e-16, #.7, 
    psi50 = x[2], #-.5, 
    b=x[3] #1.2
  )
  
  par_cost_now = list(
    alpha  = x[4], #0.13,       # cost of Jmax
    gamma = x[5] #.1,         # cost of hydraulic repair
  )
  
  dat1 = tibble(var=exp(seq(log(5),log(5000),length.out=50))) %>%
    mutate(out = purrr::map(var, ~pmodel_calibrate_analytical(tc = mean(data$T..deg.C.), ppfd = mean(data$Iabs.growth.in.growth.chamber..umol.m.2.s.1.), vpd = ., co2 = mean(data$ca..ppm.), elv = 0, fapar = .99, kphio = 0.087, psi_soil = 0, rdark = 0.02, par_plant = par_plant_now, par_cost = par_cost_now )) ) %>%    
    unnest_wider(out) %>%  
    unnest_wider(out_hydraulics, names_sep = "_") %>%  
    unnest_wider(out_analytical, names_sep = "_") 
  
  cat(species, "\t")
  p1 = dat1 %>% plot_slope() 
  plot(p1)
}

cpar$D_exp = D_exp

### VPD RESPONSE
library(tidyverse)

d = cpar #read.csv("C:/Users/Jaideep/Documents/codes/rpmodel/vignettes_add/drying_experiments_meta/fitted_params_brd_all_4_eucalyp.csv")
f = read.csv("hydricity.csv")
d = d %>% left_join(f, by = "Species")

d = d %>% filter(Species != "Eucalyptus pauciflora")


mytheme = function(){
  theme_classic()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          plot.tag.position = "topleft")  
}


poly = data.frame(x=c(-0.61, -0.61, -0.92, -0.92), y = c(0, 7.5, 7.5, 0), id = rep(1,4), value=rep(1,4))

p51 = d %>% ggplot() + 
  mytheme()+
  geom_polygon(data=poly, aes(x=x,y=y), size=1, fill=scales::alpha("green3", 0.2))+
  # geom_rect(xmin=-0.92, xmax=-0.61, ymin=-10, ymax=10,fill=scales::alpha("orange", 0.2))+
  geom_histogram(aes(x=D_exp), bins = 28, color="grey50", fill="grey80")+
  geom_vline(xintercept=-0.76, size=.8, col="green4")+
  labs(x=expression(atop("Slope of", "logit("*chi*") ~ log(D)")), 
       y="Frequency")+
  # geom_vline(xintercept = c(-0.61,-0.92), col="green4", size=.5)+
  geom_vline(xintercept = -0.5, col="orange2", size=.8)

p51

p61 = d %>% filter(Species != "Eucalyptus pauciflora") %>% ggplot(mapping = aes(x=1-dpsi_slope, y=D_exp)) +
  mytheme()+
  # scale_x_log10() + scale_y_log10() +
  geom_point(size=2) + 
  geom_smooth(method = lm, se = F) + 
  # geom_abline(slope = 1, intercept=0, color="grey")+
  ylab(expression(atop("Slope of", "logit("*chi*") ~ log(D)")))+
  xlab("Isohydricity") 

p61

p71 = d %>% filter(Species != "Eucalyptus pauciflora") %>% ggplot(mapping = aes(y=P50, x=D_exp)) +
  mytheme()+
  # scale_x_log10() + scale_y_log10() +
  geom_point(size=2) + 
  geom_smooth(method = lm, se = F) + 
  # geom_abline(slope = 1, intercept=0, color="grey")+
  xlab(expression(atop("Slope of", "logit("*chi*") ~ log(D)")))+
  ylab(expression(atop("Plant hydraulic", "vulnerability, "~psi["50"] ~ "(MPa)")))+
  scale_x_continuous(limits = c(-0.75, -0.65), breaks=c(-0.75, -0.7, -0.65), labels = c(-0.75, -0.7, -0.65))

p71


cairo_pdf(filename ="vpd_response_7_cairo.pdf", width=6.50, height=3.00)
cowplot::plot_grid(p51, p71, labels=LETTERS, label_size = 14, label_x = 0.1, label_colour = "grey50", hjust = -3, align = "hv", rel_widths = 1, ncol=2)
dev.off()

