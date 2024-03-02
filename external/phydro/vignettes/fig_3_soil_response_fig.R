library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(scales)
library(rphydro)
library(zoo)

rm(list=ls())

setwd("~/codes/phydro/vignettes/")

dat = read.csv(file="pred_obs_a_g_c.csv")
names = colnames(dat)
# [5] "Predawn.LWP..MPa."                                 
# [6] "A..umol.m.2.s.1."                                  
# [7] "gC..mol.m.2.s.1."                                  
# [8] "ca..ppm."                                          
# [9] "D..unitless...Pa...Pa.atm.."                       
# [10] "T..deg.C."                                         
names[5:10] = c("LWP", "A", "gC", "ca", "D", "T")
names[14] = "Iabs_growth"
names[13] = "Iabs_used"
colnames(dat) = names

dpsi_df = read.csv(file = "data/drying_experiments_dpsi_Joshi_et_al_2022.csv")
cpar = read.csv("data/fitted_params_Joshi_et_al_2022.csv") 


source("utils.R")


pmodel_fit = function(species, dpsi_data=NULL){
  
  data=filter(dat, Species==species)
  x = as.numeric(cpar[which(cpar$Species == species), -1])
  
  par_plant_now = list(
    # Ks0=1, 
    # v_huber=1, 
    conductivity= x[1]*1e-16, # in Joshi et al 2021, K is in  10^-16 m 
    # height=1, 
    psi50 = x[2],  
    b=x[3] 
  )
  
  par_cost_now = list(
    alpha  = x[4],          # cost of Jmax
    gamma = x[5]    # cost of hydraulic repair.  
  )
  
  
  if (cpar[which(cpar$Species == species), "inst"] == 0){
    cat("Acclimated response: days = ", data$Drydown.days, "\n")
    dat1 = tibble(var = seq(-6,0, length.out = 30)) %>% 
      mutate(p = purrr::map(var, ~pmodel_calibrate_analytical(tc = mean(data$T), ppfd = mean(data$Iabs_growth), vpd = mean(data$D*101325), co2 = mean(data$ca), elv = 0, fapar = .99, kphio = 0.087, psi_soil = ., rdark = 0.02, par_plant = par_plant_now, par_cost = par_cost_now)) ) %>% unnest_wider(p) %>% unnest_wider(out_hydraulics) %>% unnest_wider(out_analytical, names_sep = "_") 
  } else {
    ndays = mean(data$Drydown.days)
    psi_min = min(data$LWP)
    psi_max = max(data$LWP)
    # cat(ndays,"\n")
    
    lwp = seq(-6,0, length.out=30)
    day = ndays * (lwp-psi_max)/(psi_min-psi_max)
    
    lwp_day = function(day_num){
      psi_max + day_num/ndays * (psi_min-psi_max)
    } 
    
    k = 7
    lwp_week = rollmean(x = lwp_day(c(max(day):0, rep(0,k-1))), k = k, align = "right")
    
    cat("Instantaneous response with slow acclimation: days = ", data$Drydown.days, ", k = ", k, "\n")
    
    spl = splinefun(x = max(day):0, y=lwp_week)
    
    # plot(x=day, y=lwp_day(day))
    # points(x=day, y=spl(day), col="red", type="o")
    # points(x=max(day):0, y=lwp_week, col="blue", type="l")
    
    dat_acc = tibble(var = spl(day)) %>% mutate(pmod = map(var, ~rphydro_analytical(tc = mean(data$T), ppfd = mean(data$Iabs_growth), vpd = mean(data$D*101325), co2 = mean(data$ca), elv = 0, fapar = .99, kphio = 0.087, psi_soil = ., rdark = 0.02, par_plant=par_plant_now, par_cost = par_cost_now, opt_hypothesis = "PM"))) %>% unnest_wider(pmod)
    dat1 = tibble(var = lwp, jmax_a=dat_acc$jmax, vcmax_a=dat_acc$vcmax) %>%
      mutate(p = purrr::pmap(list(var, jmax_a, vcmax_a), ~pmodel_calibrate_inst(tc = mean(data$T), ppfd = mean(data$Iabs_used), vpd = mean(data$D*101325), co2 = mean(data$ca), elv = 0, fapar = .99, kphio = 0.087, psi_soil = ..1, rdark = 0.02, par_plant=par_plant_now, par_cost = par_cost_now, jmax = ..2, vcmax = ..3, opt_hypothesis = "PM")) ) %>% unnest_wider(p) %>% unnest_wider(out_hydraulics) %>% unnest_wider(out_analytical, names_sep = "_")
  }
  

  p88 = (log(0.12)/log(.5))^(1/par_plant_now$b)*par_plant_now$psi50
  cat("Profit vector = ", dat1$profit[dat1$var > p88], "\n")
  cat("Mean Profit (P88) = ", mean(dat1$profit[dat1$var > p88]), "\n")
  cat("Mean Profit (P50) = ", mean(dat1$profit[dat1$var > par_plant_now$psi50]), "\n")
  
  # p88 = (log(0.12)/log(.5))^(1/par_plant_now$b)*par_plant_now$psi50
  # profit = dat2$a - par_cost_now$c*dat2$jmax - par_cost_now$a1*(dat2$dpsi)^2/par_plant_now$psi50^2  
  # cat("Profit @ 1200 PPFD = ", mean(profit[dat2$var > par_plant_now$psi50]) )
  
  dat1
}





plot_list = function(n, results, varname, specieslist, dpsilist, nulllist=NULL){
  
  p1 = ggplot()
  p2 = ggplot()
  p3 = ggplot()
  p4 = ggplot()
  p5 = ggplot()
  p6 = ggplot()
  
  ptype = c(15,17,16)
  
  colsa = colorRampPalette(colors = c("chartreuse3", "chartreuse4"))(n)
  colsgs = colorRampPalette(colors = c("cyan2", "cyan4"))(n)
  colsx = colorRampPalette(colors = c("magenta", "magenta4"))(n)
  colspsi = colorRampPalette(colors = c("skyblue", "skyblue4"))(n)
  colsvcmax = colorRampPalette(colors = c("seagreen1", "seagreen4"))(n)
  colsjmax = colorRampPalette(colors = c("goldenrod1", "goldenrod4"))(n)
  
  
  
  
  for (i in 1:n){
    
    df = results[[i]]
    species = specieslist[[i]]
    dpsi_data = dpsilist[[i]]
    
    if (!is.null(nulllist)) dfnull = nulllist[[i]]
    
    subdata = dat %>% filter(Species == species) %>% filter(LWP > -6)
    
    mytheme = function(){
      theme_classic()+
        theme(axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              plot.tag.position = "topleft")
    }
    
    colsa = viridis_pal(begin=.2, end=.9)(5)
    cols_p = scales::alpha(colsa[seq(1,n*2,2)], 0.75)
    cols_l = scales::alpha(colsa[seq(2,n*2,2)], 1)
    p1 <- p1 +
      mytheme() + 
      geom_line(data = df, mapping = aes(x = var, y = vcmax), col=cols_l[i], size=1) +
      expand_limits(y=0)+
      xlab(varname)+
      # ylab(expression(V["cmax"]~"("*mu*"mol m"^-2~"s"^"-1"*")"))
      ylab(expression(atop(atop(
        textstyle("Carboxylation"), 
        textstyle("capacity,")), 
        textstyle(V["cmax"]~"("*mu*"mol m"^-2~"s"^-1*")") 
      )))  
    if (!is.null(nulllist))
      p1 <- p1 + geom_line(data = dfnull, mapping = aes(x = var, y = vcmax), col=(cols_l[i]), size=0.2)
    
    p2 <- p2 + 
      mytheme() + 
      geom_line(data = df, mapping = aes(x = var, y = dpsi), col=cols_l[i], size=1)+
      # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
      expand_limits(y=0)+
      xlab(varname)+
      # ylab(expression(Delta*psi~"(MPa)"))
      ylab(expression(atop(atop(
        textstyle("Soil-leaf water-"), 
        textstyle("potential difference,")), 
        textstyle(Delta*psi~"(MPa)") 
      ))) 
    if (!is.null(dpsi_data)){
      p2 <- p2 +
        geom_point(data=dpsi_data, aes(x=SWP, y=Dpsi), shape=ptype[i], col=cols_p[i])
    }
    if (!is.null(nulllist))
      p2 = p2 + geom_line(data = dfnull, mapping = aes(x = var, y = dpsi), col=cols_l[i], size=0.2)
    
    p3 <- p3 +
      mytheme() + 
      geom_line(data = df, mapping = aes(x = var, y = gs), col=cols_l[i], size=1)+
      # geom_line(data = df, mapping = aes(x = var, y = out_analytical_gs), col="grey", size=1)+
      geom_point(data=subdata, aes(x=LWP, y=gC), shape=ptype[i], col=cols_p[i])+
      expand_limits(y=0)+
      xlab(varname)+
      # ylab(expression(g[s]~"(mol m"^-2~"s"^"-1"*")"))
      ylab(expression(atop(atop(
        textstyle("Stomatal"), 
        textstyle("conductance,")), 
        textstyle(g[s]~"(mol m"^-2~"s"^-1*")") 
      ))) 
    
    if (!is.null(nulllist))
      p3 = p3 + geom_line(data = dfnull, mapping = aes(x = var, y = gs), col=cols_l[i], size=0.2)
    
    p4 <- p4 +
      mytheme() + 
      geom_line(data = df, mapping = aes(x = var, y = chi), col=cols_l[i], size=1)+
      # geom_line(data = df, mapping = aes(x = var, y = out_analytical_chi), col="grey", size=1)+
      geom_point(data=subdata, aes(x=LWP, y=Ciest/ca), shape=ptype[i], col=cols_p[i])+
      # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
      expand_limits(y=0)+
      xlab(varname)+
      # ylab(expression(chi))
      ylab(expression(atop(atop(
        textstyle("Leaf internal-to-"), 
        textstyle("external CO"[2]*" ratio,")), 
        textstyle(chi) 
      ))) 
    if (!is.null(nulllist))
      p4 = p4 + geom_line(data = dfnull, mapping = aes(x = var, y = chi), col=cols_l[i], size=0.2)
    
    p5 <- p5 +
      mytheme() + 
      geom_line(data = df, mapping = aes(x = var, y = jmax), col=cols_l[i], size=1) +
      # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
      # geom_point(data=filter(dat, Species==species), aes(x=LWP, y=Jmax))+
      expand_limits(y=0)+
      xlab(varname)+
      # ylab(expression(J["max"]~"("*mu*"mol m"^-2~"s"^"-1"*")"))
      ylab(expression(atop(atop(
        textstyle("Electron-transport"), 
        textstyle("Capacity,")), 
        textstyle(J["max"]~"("*mu*"mol m"^-2~"s"^-1*")") 
      ))) 
    if (!is.null(nulllist))
      p5 = p5 + geom_line(data = dfnull, mapping = aes(x = var, y = jmax), col=cols_l[i], size=0.2) 
    
    p6 <- p6 +
      mytheme() + 
      geom_line(data = df, mapping = aes(x = var, y = a), col=cols_l[i], size=1) +
      # geom_line(data = df, mapping = aes(x = var, y = out_analytical_gpp), col="grey",size=1) +
      geom_point(data=subdata, aes(x=LWP, y=A), shape=ptype[i], col=cols_p[i])+
      # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
      expand_limits(y=0)+
      xlab(varname)+
      # ylab(expression(A~"("*mu*"mol m"^-2~"s"^"-1"*")"))
      ylab(expression(atop(atop(
        textstyle("Assimilation"), 
        textstyle("rate,")), 
        textstyle(A~"("*mu*"mol m"^-2~"s"^-1*")") 
      ))) 
    if (!is.null(nulllist))
      p6 <- p6 + geom_line(data = dfnull, mapping = aes(x = var, y = a), col=cols_l[i], size=0.2)
    
  }
  
  # # p6
  # # grid.arrange(p1,p2,p3,p4,p5,p6, ncol=2)
  # cowplot::plot_grid(p6, p3, p4, p2, p1, p5, labels="AUTO", label_size = 16, label_colour = "grey50", label_x = 0.3, hjust = 0, align = "hv", rel_widths = 1, ncol=2)
  
  
  list(p6,p3,p4, p2,p1,p5)  
  # grid.arrange(p3,p6,p1,p2,p5,p4, ncol=2)
}




make_plots = function(spp_list, dpsi_files_list){

  dat_list = list()
  for (spp in spp_list){
    dat_list = c(dat_list, list(pmodel_fit(spp)))
  }
  
  dpsi_list = list()
  for (dpsifile in dpsi_files_list){
    if (dpsifile != "") dpsi = read.csv(paste0("drying_experiments_meta/",dpsifile,".csv"))
    else dpsi = NULL
    dpsi_list = c(dpsi_list, list(dpsi))
  }
  
  l <<- c(l, plot_list(n = length(spp_list), results = dat_list, 
                     varname = expression(atop(
                       textstyle("Soil water potential,"), 
                       textstyle(psi[s]~"(MPa)") 
                     )),
                     specieslist = spp_list, dpsilist = dpsi_list))
  
}


#### Main text fig ####

l = list()

dat1 = pmodel_fit("Eucalyptus pilularis")
dat2 = pmodel_fit("Eucalyptus populnea")
dat1$chi[dat1$var < -3.5]=NA

dpsi1 = dpsi_df %>% filter(Species == "Eucalyptus pilularis")
dpsi2 = dpsi_df %>% filter(Species == "Eucalyptus populnea")

l = c(l, plot_list(n = 2, results = list(dat1, dat2), 
                   varname = expression(atop(
                     textstyle("Soil water potential,"), 
                     textstyle(psi[s]~"(MPa)") 
                   )),
                   specieslist = list("Eucalyptus pilularis", "Eucalyptus populnea"), dpsilist = list(dpsi1, dpsi2)))
# dev.off()

cairo_pdf(filename = "eucalyptus_final_cairo.pdf", width = 7.4, height=7.9)
cowplot::plot_grid(plotlist = l, labels="AUTO", label_size = 16, label_colour = "grey50", label_x = 0.38, hjust = 0, align = "hv", rel_widths = 1, ncol=2)
dev.off()



#### SI figs ####

plot_list = function(n, results, varname, specieslist, dpsilist, nulllist=NULL){
  
  p1 = ggplot()
  p2 = ggplot()
  p3 = ggplot()
  p4 = ggplot()
  p5 = ggplot()
  p6 = ggplot()
  
  ptype = c(15,17,16)
  
  colsa = colorRampPalette(colors = c("chartreuse3", "chartreuse4"))(n)
  colsgs = colorRampPalette(colors = c("cyan2", "cyan4"))(n)
  colsx = colorRampPalette(colors = c("magenta", "magenta4"))(n)
  colspsi = colorRampPalette(colors = c("skyblue", "skyblue4"))(n)
  colsvcmax = colorRampPalette(colors = c("seagreen1", "seagreen4"))(n)
  colsjmax = colorRampPalette(colors = c("goldenrod1", "goldenrod4"))(n)
  
  
  
  
  for (i in 1:n){
    cat(i, "\n")
    df = results[[i]]
    species = specieslist[[i]]
    dpsi_data = dpsilist[[i]]
    
    if (!is.null(nulllist)) dfnull = nulllist[[i]]
    
    subdata = dat %>% filter(Species == species) %>% filter(LWP > -6)
    
    mytheme = function(){
      theme_classic()+
        theme(axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              plot.tag.position = "topleft")
    }
    
    colsa = viridis_pal(begin=.2, end=.9)(6)
    cols_p = scales::alpha(colsa[seq(1,n*2,2)], 0.75)
    cols_l = scales::alpha(colsa[seq(2,n*2,2)], 1)
    p1 <- p1 +
      mytheme() + 
      geom_line(data = df, mapping = aes(x = var, y = vcmax), col=cols_l[i], size=1) +
      expand_limits(y=0)+
      xlab(varname)+
      # ylab(expression(V["cmax"]~"("*mu*"mol m"^-2~"s"^"-1"*")"))
      ylab(expression(atop(atop(
        textstyle("Carboxylation"), 
        textstyle("capacity,")), 
        textstyle(V["cmax"]~"("*mu*"mol m"^-2~"s"^-1*")") 
      )))  
    if (!is.null(nulllist))
      p1 <- p1 + geom_line(data = dfnull, mapping = aes(x = var, y = vcmax), col=(cols_l[i]), size=0.2)
    
    p2 <- p2 + 
      mytheme() + 
      geom_line(data = df, mapping = aes(x = var, y = dpsi), col=cols_l[i], size=1)+
      # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
      expand_limits(y=0)+
      xlab(varname)+
      # ylab(expression(Delta*psi~"(MPa)"))
      ylab(expression(atop(atop(
        textstyle("Soil-leaf water-"), 
        textstyle("potential difference,")), 
        textstyle(Delta*psi~"(MPa)") 
      ))) 
    if (!is.null(dpsi_data)){
      p2 <- p2 +
        geom_point(data=dpsi_data, aes(x=SWP, y=Dpsi), shape=ptype[i], col=cols_p[i])
    }
    if (!is.null(nulllist))
      p2 = p2 + geom_line(data = dfnull, mapping = aes(x = var, y = dpsi), col=cols_l[i], size=0.2)
    
    p3 <- p3 +
      mytheme() + 
      geom_line(data = df, mapping = aes(x = var, y = gs), col=cols_l[i], size=1)+
      # geom_line(data = df, mapping = aes(x = var, y = out_analytical_gs), col="grey", size=1)+
      geom_point(data=subdata, aes(x=LWP, y=gC), shape=ptype[i], col=cols_p[i])+
      expand_limits(y=0)+
      xlab(varname)+
      # ylab(expression(g[s]~"(mol m"^-2~"s"^"-1"*")"))
      ylab(expression(atop(atop(
        textstyle("Stomatal"), 
        textstyle("conductance,")), 
        textstyle(g[s]~"(mol m"^-2~"s"^-1*")") 
      ))) 
    
    if (!is.null(nulllist))
      p3 = p3 + geom_line(data = dfnull, mapping = aes(x = var, y = gs), col=cols_l[i], size=0.2)
    
    p4 <- p4 +
      mytheme() + 
      geom_line(data = df, mapping = aes(x = var, y = chi), col=cols_l[i], size=1)+
      # geom_line(data = df, mapping = aes(x = var, y = out_analytical_chi), col="grey", size=1)+
      geom_point(data=subdata, aes(x=LWP, y=Ciest/ca), shape=ptype[i], col=cols_p[i])+
      # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
      expand_limits(y=0)+
      xlab(varname)+
      # ylab(expression(chi))
      ylab(expression(atop(atop(
        textstyle("Leaf internal-to-"), 
        textstyle("external CO"[2]*" ratio,")), 
        textstyle(chi) 
      ))) 
    if (!is.null(nulllist))
      p4 = p4 + geom_line(data = dfnull, mapping = aes(x = var, y = chi), col=cols_l[i], size=0.2)
    
    p5 <- p5 +
      mytheme() + 
      geom_line(data = df, mapping = aes(x = var, y = jmax), col=cols_l[i], size=1) +
      # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
      # geom_point(data=filter(dat, Species==species), aes(x=LWP, y=Jmax))+
      expand_limits(y=0)+
      xlab(varname)+
      # ylab(expression(J["max"]~"("*mu*"mol m"^-2~"s"^"-1"*")"))
      ylab(expression(atop(atop(
        textstyle("Electron-transport"), 
        textstyle("capacity,")), 
        textstyle(J["max"]~"("*mu*"mol m"^-2~"s"^-1*")") 
      ))) 
    if (!is.null(nulllist))
      p5 = p5 + geom_line(data = dfnull, mapping = aes(x = var, y = jmax), col=cols_l[i], size=0.2) 
    
    p6 <- p6 +
      mytheme() + 
      geom_line(data = df, mapping = aes(x = var, y = a), col=cols_l[i], size=1) +
      # geom_line(data = df, mapping = aes(x = var, y = out_analytical_gpp), col="grey",size=1) +
      geom_point(data=subdata, aes(x=LWP, y=A), shape=ptype[i], col=cols_p[i])+
      # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
      expand_limits(y=0)+
      xlab(varname)+
      # ylab(expression(A~"("*mu*"mol m"^-2~"s"^"-1"*")"))
      ylab(expression(atop(atop(
        textstyle("Assimilation"), 
        textstyle("rate,")), 
        textstyle(A~"("*mu*"mol m"^-2~"s"^-1*")") 
      ))) 
    if (!is.null(nulllist))
      p6 <- p6 + geom_line(data = dfnull, mapping = aes(x = var, y = a), col=cols_l[i], size=0.2)
    
  }
  
  # # p6
  # # grid.arrange(p1,p2,p3,p4,p5,p6, ncol=2)
  # cowplot::plot_grid(p6, p3, p4, p2, p1, p5, labels="AUTO", label_size = 16, label_colour = "grey50", label_x = 0.3, hjust = 0, align = "hv", rel_widths = 1, ncol=2)
  
  
  list(p6,p3,p4, p2,p1,p5)  
  # grid.arrange(p3,p6,p1,p2,p5,p4, ncol=2)
}

l = list()

spp_list = list("Eucalyptus pilularis", "Eucalyptus populnea")
dpsi_files_list = spp_list
make_plots(spp_list, dpsi_files_list)

spp_list = list("Quercus coccifera", "Quercus ilex", "Quercus suber")
dpsi_files_list = spp_list
make_plots(spp_list, dpsi_files_list)

spp_list = list("Broussonetia papyrifera (Linnaeus) L_He ritier ex Ventenat", "Platycarya longipes Wu", "Pteroceltis tatarinowii Maximowicz")
dpsi_files_list = list("","","")
make_plots(spp_list, dpsi_files_list)




spp_list = list("Olea europaea var. Meski", "Olea europaea var. Chemlali")
dpsi_files_list = list("", "")
make_plots(spp_list, dpsi_files_list)

spp_list = list("Helianthus annuus", "Glycine max")
dpsi_files_list = list("", "Glycine max")
make_plots(spp_list, dpsi_files_list)

spp_list = list("Cedrus atlantica", "Pseudotzuga menziesii")
dpsi_files_list = list("", "Pseudotzuga menziesii")
make_plots(spp_list, dpsi_files_list)



spp_list = list("Rosa cymosa Trattinnick", "Ficus tikoua")
dpsi_files_list = list("", "")
make_plots(spp_list, dpsi_files_list)

spp_list = list("Allocasuarina luehmannii")
dpsi_files_list = list("Allocasuarina luehmannii")
make_plots(spp_list, dpsi_files_list)

spp_list = list("Cinnamomum bodinieri H. Leveille")
dpsi_files_list = list("")
make_plots(spp_list, dpsi_files_list)



nsets = length(l)/6
ncols = 3

for (i in seq(1,ceiling(nsets/ncols), 1)){
  lsub = l[(ncols*6*(i-1)+1) : min((ncols*6*(i)), length(l))]
  ncol = length(lsub)/6
  png(paste0("species_response_2_",i,".png"), height = 1200*3.5, width = 300*ncol*3.5, res=300)
  cowplot::plot_grid(plotlist = lsub, labels="", label_size = 16, label_colour = "grey50", label_x = 0.3, hjust = 0, align = "hv", rel_widths = 1, nrow=6, byrow=F)
  dev.off()
}
