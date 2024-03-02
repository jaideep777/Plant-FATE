library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(scales)

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
colnames(dat) = names

dpsi_df = read.csv(file = "data/drying_experiments_dpsi_Joshi_et_al_2022.csv")
cpar = read.csv("data/fitted_params_Joshi_et_al_2022.csv") 

# # Remove because < 10 data points
# dpsi_df = dpsi_df %>% filter(!(Species %in% c("Quercus coccifera", "Quercus suber")))

dat = left_join(x = dat, y= cpar, by="Species")
dpsi_df = left_join(x = dpsi_df, y= cpar, by="Species")

# remove pauciflora
dat = dat %>% filter(Species != "Eucalyptus pauciflora")
dpsi_df = dpsi_df %>% filter(Species != "Eucalyptus pauciflora")
# # 
# dat = dat %>% filter(Species != "Helianthus annuus")
# dpsi_df = dpsi_df %>% filter(Species != "Helianthus annuus")


dat$chi = dat$Ciest/dat$ca

dpsi_df = dpsi_df %>% mutate(source = factor(source, levels = c(0,1), labels = c("Experiment", "Literature") ))

makeTransparent = function(col, alpha=0.7){
  rgb(t(col2rgb(col)/255), alpha = alpha)
}

mytheme = function(){
  theme_classic()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          legend.text=element_text(size=11),
          plot.tag.position = "topleft") 
    
}

p1 =
  dat %>% ggplot(mapping = aes(x=a_pred, y=A)) + 
  mytheme()+
  geom_point(aes(fill=-LWP/pgL88), shape=21, color=makeTransparent("black",0)) + 
  geom_smooth(method = "lm", se = F, col="black")+
  geom_abline(slope = 1, intercept=0, color="grey") +
  # scale_fill_gradient(oob=squish, limits = c(-4,0), low = makeTransparent("goldenrod1"), high = makeTransparent("blue"))+
  scale_fill_viridis_c(oob=squish, limits = c(-1,0), direction = -1, begin=.25)+
  labs(fill = expression(psi[s]*"/|"*psi[g88]*"|"))+
  annotate(size=4, geom="text", label = paste("italic(R)^2 == ", sprintf("%.2f", cor(dat$A, dat$a_pred, use = "pairwise.complete.obs")^2)), x=Inf, y=-Inf, hjust=1.2, vjust=-1.2, parse=T)+
  xlab(expression(textstyle("Predicted A")))+
  ylab(expression(textstyle("Observed A")))+
  ggtitle(expression(atop("Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))

p1 


p2 = dat %>% ggplot(mapping = aes(x=g_pred, y=gC)) +   
  mytheme()+
  geom_point(aes(fill=-LWP/pgL88), shape=21, color=makeTransparent("black",0)) + 
  geom_smooth(method = "lm", se = F, col="black")+
  geom_abline(slope = 1, intercept=0, color="grey") +
  # scale_fill_gradient(oob=squish, limits = c(-4,0), low = makeTransparent("goldenrod1"), high = makeTransparent("blue"))+
  scale_fill_viridis_c(oob=squish, limits = c(-1,0), direction = -1, begin=.25)+
  labs(fill = expression(psi[s]*"/|"*psi[g88]*"|"))+
  annotate(size=4, geom="text", label = paste("italic(R)^2 == ", sprintf("%.2f", cor(dat$gC, dat$g_pred, use = "p")^2)), x=Inf, y=-Inf, hjust=1.2, vjust=-1.2, parse=T)+
  xlab(expression("Predicted g"[s]))+
  ylab(expression("Observed g"[s]))+
  ggtitle(expression(atop("Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))
  # ylim(0, 0.3)

# scale_x_log10()+
  # scale_y_log10()
p2 



p3 = dat %>% ggplot(mapping = aes(x=c_pred, y=Ciest/ca)) +   
  mytheme()+
  geom_point(aes(fill=-LWP/pgL88), shape=21, color=makeTransparent("black",0)) + 
  geom_smooth(data = dat[which(dat$LWP/dat$pgL88 < 1),], method = "lm", se = F, col="black")+
  geom_abline(slope = 1, intercept=0, color="grey") +
  # scale_fill_gradient(oob=squish, limits = c(-4,0), low = makeTransparent("goldenrod1"), high = makeTransparent("blue"))+
  scale_fill_viridis_c(oob=squish, limits = c(-1,0), direction = -1, begin=.25)+
  annotate(size=4, geom="text", label = paste("italic(R)^2 == ", sprintf("%.2f", cor(dat$c_pred[-dat$LWP/dat$pgL88>-1], (dat$Ciest/dat$ca)[-dat$LWP/dat$pgL88>-1], use = "pairwise.complete.obs")^2)), x=Inf, y=-Inf, hjust=1.2, vjust=-1.2, parse=T)+
  labs(fill = expression(psi[s]*"/|"*psi[g88]*"|"))+
  xlab(expression("Predicted"~chi))+
  ylab(expression("Observed"~chi))+
  ggtitle(expression(atop("Leaf internal-to-","ambient CO"[2]~"ratio,"~italic(chi))))

p3



p4 = dpsi_df %>% ggplot() +   
  mytheme()+
  geom_point(mapping = aes(x=d_pred, y=Dpsi, color=-SWP/pgL88, shape=source, size = source)) + 
  geom_smooth(mapping=aes(x=d_pred, y=Dpsi), method = "lm", se = F, col="black")+
  geom_abline(slope = 1, intercept=0, color="grey") +
  # scale_color_gradient(low = makeTransparent("goldenrod1"), high = makeTransparent("blue"))+
  scale_color_viridis_c(oob=squish, limits = c(-1,0), direction = -1, begin=.25)+
  scale_size_manual(values=c(1.5,1))+
  annotate(size=4, geom="text", label = paste("italic(R)^2 == ", sprintf("%.2f", cor(dpsi_df$Dpsi, dpsi_df$d_pred, use = "pairwise.complete.obs")^2)), x=Inf, y=-Inf, hjust=1.2, vjust=-1.2, parse=T)+
  xlab(expression("Predicted"~Delta*psi))+
  ylab(expression("Observed"~Delta*psi)) + 
  ggtitle(expression(atop("Soil-leaf water-potential","difference,"~italic(Delta*psi)~"(MPa)")))+
  labs(color=expression(psi[s]*"/|"*psi[g88]*"|"), size="Source", shape="Source")+
  ylim(c(0,NA))
p4

cairo_pdf(filename = "calibration_validation_viridis_wtrend_g88_4_cairo.pdf", width = 8.20, height = 6.50)
cowplot::plot_grid(p1, p2, p3, p4, labels=LETTERS, label_size = 16, label_colour = "grey50", label_x = 0.05, align = "hv", rel_widths = 1, ncol=2)
dev.off()


# 
# x = log(dat$gC)
# y = log(dat$g_pred)
# cor(x,y, use = "p")

