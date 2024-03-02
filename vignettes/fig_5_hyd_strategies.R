rm(list=ls())
library(MASS)
library(tidyverse)

setwd("~/codes/phydro/vignettes")

d = read.csv("data/fitted_params_Joshi_et_al_2022.csv")
d = d %>% mutate(yxp50 = gamma*P50^2)
d = d %>% mutate(P88 = (log(0.12)/log(0.5))^(1/b)*P50 )
d = d %>% dplyr::select(-Ref.1, -Ref2, -Ref, -Ref.SLA, -b..Cavit., -b..SC., -Pclose, -Slope, -P12, -P50..stem., -PgS88, -PgS50, -pg12, -H)
d = d %>% dplyr::select(-SLA_nopetiole, -SLA_wpetiole, -SLA1, -E.D)
d = d %>% dplyr::select(-Pmin, -dpsi_int)

# d = d %>% filter(!(Species %in% c("Eucalyptus pauciflora"))) 

f = read.csv("hydricity.csv")

d = d %>% left_join(f, by = "Species")

d = d %>% left_join(d %>% group_by(A.G) %>% summarise(npoints = length(A.G)), by="A.G")

# d = d %>% filter(dpsi_cal == 1)

#### Fig 4 ####


my_theme = function(){
  theme_classic()+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          legend.position="none",
          plot.tag.position = "topleft")
}

spp_cols = scales::hue_pal()(5)  #c("red", "green4", "yellowgreen", "cyan2", "blue")


datme = read.csv(file = "hydricity.csv")
datmm = read.csv(file = "C:/Users/Jaideep/Documents/zotero_data/storage/NXGJVVYW/data.csv")

q1 = ggplot() + my_theme() + 
  geom_histogram(data=datmm, aes(x=s..MPa.MPa.1.), alpha=0.4, bins=10, fill="green3") + 
  geom_histogram(data=datme, aes(x=1-dpsi_slope), fill="grey50", col="black", bins=10, alpha=0.7)+
  ylab("Number of species")+
  xlab("Anisohydricity")
q1


weird <- scales::trans_new("signed_log",
                                transform=function(x) -log((-x)),
                                inverse=function(x) -exp(-(x)))

q2 = d %>% filter(K.scalar<10) %>% ggplot(mapping = aes(y=K.scalar, x=-P50, color=A.G)) +
  my_theme()+
  # scale_x_continuous(trans="log", breaks = c(0.5, 1, 2, 4)) + scale_y_log10() + 
  geom_point(aes(color=A.G, shape = factor(dpsi_cal)), size=3, stroke=1) + 
  scale_shape_manual(values = c(21, 16)) +
  scale_color_manual(values = spp_cols)+
  # geom_smooth(aes(group=A.G), method = lm, se = F) + 
  geom_smooth(method = lm, se = F, col="black") + 
  # xlab(expression(-psi[50]))+
  # ylab(expression(K)) +
  xlab(expression(atop(
    textstyle("Plant hydraulic"), 
    textstyle("vulnerability,"~-psi[50]~"(MPa)") 
  ))) +
  ylab(expression(atop(
    textstyle("Plant conductance,"), 
    textstyle(K[p]~"(10"^"-16"~m*")") 
  ))) #+
  # expand_limits(y=10)
q2

q3 = d %>% ggplot(mapping = aes(x=pgL88, y=P50X)) +
  my_theme()+
  geom_point(aes(colour=A.G, shape=factor(dpsi_cal)), size=3, stroke=1) + 
  geom_smooth(method = rlm, se = F, col="black") + 
  scale_shape_manual(values = c(21, 16)) +
  scale_color_manual(values = spp_cols)+
  geom_abline(slope = 1, intercept=0, color="grey")+
  # ylab(expression(tilde(psi)["50X"]))+
  # xlab(expression(psi["g88"])) 
  xlab(expression(atop(
    textstyle("Stomatal closure point,"), 
    textstyle(psi["g88"]~"(MPa)") 
  ))) +
  ylab(expression(atop(
    textstyle("Xylem hydraulic"), 
    textstyle("vulnerability,"~tilde(psi)["50x"]~"(MPa)") 
  ))) 
# geom_abline(slope = 1, intercept = 0, col="grey")
q3


# subd = d %>% filter(!is.na(Ptlp)) %>% arrange(desc(pgL88))
# plot(subd$Ptlp, ylim=c(-8,0))
# points(subd$P50, type="l", col="red")
# points(subd$pgL88, type="l", col="green")

subd = d %>% filter(!is.na(Ptlp)) #%>% arrange(desc(Species)) 

l = lapply(subd$Species, FUN = function(x){unlist(strsplit(x, " "))[1:2]})
l1 = transpose(l) %>% map(unlist)
genus = paste0(substr(l1[[1]], 1,1), ".")
spp = substr(l1[[2]],1,3)
spp_short = paste0(genus, " ", spp, ".")

q4 = subd %>% 
  ggplot() + 
  my_theme() + 
  theme(axis.text.x=element_text(angle = 90, hjust = 1))+
  geom_ribbon(aes(x=1:7, ymin=pgL88, ymax=P50), fill="grey60", alpha=0.5)+
  geom_line(aes(x=1:7, y=pgL88)) + 
  geom_line(aes(x=1:7, y=P50)) + 
  geom_line(aes(x=1:7, y=Ptlp), col="green3", size=1.5) + 
  ylab(expression(atop(
    textstyle("Turgor loss point,"),
    textstyle(tilde(psi)["tlp"]~"(MPa)")
  ))) + 
  xlab("Species") + 
  scale_x_continuous(breaks=1:7, labels = spp_short)
q4    

q4_new= subd %>% 
  ggplot() + 
  my_theme() + 
  theme(axis.text.x=element_text(angle = 90, hjust = 1))+
  geom_crossbar(aes(x=1:7, y=Ptlp, ymin=pgL88, ymax=P50), fill="grey60", alpha=0.5, fatten=0, col=NA)+
  geom_crossbar(aes(x=1:7, y=Ptlp, ymin=Ptlp, ymax=Ptlp), fill="grey60", alpha=0.5, fatten=3, col="green3")+
  geom_crossbar(aes(x=1:7, y=P50, ymin=P50, ymax=P50), fill="grey60", alpha=0.5, fatten=2, col="black")+
  geom_crossbar(aes(x=1:7, y=pgL88, ymin=pgL88, ymax=pgL88), fill="grey60", alpha=0.5, fatten=2, col="brown")+
  ylab(expression(atop(
    textstyle("Turgor loss point,"),
    textstyle(tilde(psi)["tlp"]~"(MPa)")
  ))) + 
  xlab("Species") + 
  scale_x_continuous(breaks=1:7, labels = spp_short)

q4_new
  
cairo_pdf(filename = "hyd_strategies_fig_9_cairo.pdf", height=5*1.7, width=9)
cowplot::plot_grid(q1, q4_new, q2, q3, labels="AUTO", label_size = 18, label_colour = "grey40", label_x = 0.2, align = "hv", rel_widths = 1)
dev.off()


q1.1 = ggplot() + my_theme() + 
  geom_histogram(data=datmm, aes(x=L..MPa.), alpha=0.4, bins=10, fill="green3") + 
  geom_histogram(data=datme, aes(x=-dpsi_intercept), fill="grey50", col="black", bins=10, alpha=0.7)+
  ylab("Number of species")+
  xlab("Leaf water potential\nunder wet conditions")
q1.1

png(filename = "dpsi_intercept_SI_fig.png", height=550*3, width=500*3, res=300)
q1.1
dev.off()






mytheme = function(){
  theme_classic()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          plot.tag.position = "topleft")  
}


p1 = d %>% filter(Species != "Helianthus annuus") %>% ggplot(mapping = aes(x=alpha, y=gamma, group=A.G, color=A.G)) +
  my_theme()+
  # scale_x_log10() + 
  scale_y_log10() +
  geom_point(aes(shape=factor(dpsi_cal)), size=3, stroke=1) + 
  geom_hline(data = d %>% group_by(A.G) %>% summarise(ym=mean(gamma)), mapping = aes(yintercept=ym, color=A.G))+
  scale_shape_manual(values = c(21, 16)) +
  scale_color_manual(values = spp_cols)+
  # geom_smooth(method = "rlm", se = F) + 
  # xlab(expression(gamma%*%psi[50]^"2"))+
  # ylab(expression(alpha))
  ylab(expression(atop(atop(
    textstyle("Unit cost of"), 
    textstyle("hydraulic capacity,")), 
    textstyle(gamma~"(mol m"^-2*"s"^-1*"MPa"^-2*")") 
  ))) +
  xlab(expression(atop(atop(
    textstyle("Unit cost of"), 
    textstyle("photosynthetic capacity,")), 
    textstyle(alpha) 
  ))) +
  scale_x_continuous(limits = c(0.065, 0.125), breaks=c(0.07, 0.095, 0.12))

p1


# %>% filter(dpsi_cal == 1 & A.G == "Angiosperm")
p5 = d %>% filter(Species != "Helianthus annuus") %>% ggplot(mapping = aes(x=K.scalar, y=SLA, colour=A.G)) +
  my_theme()+
  # scale_x_log10() + 
  # scale_y_log10() + 
  geom_point(aes(shape=factor(dpsi_cal)), size=3, stroke=1) + 
  scale_shape_manual(values = c(21, 16)) +
  scale_color_manual(values = spp_cols)+
  geom_smooth(aes(group=A.G), method = rlm, se = F) + 
  geom_smooth(method = rlm, se = F, col="black") + 
  # ylab(expression(tilde("SLA")))+
  # xlab(expression("K")) 
  xlab(expression(atop(atop(
    textstyle("Leaf specific"), 
    textstyle("plant conductance,")), 
    textstyle(K[p]~"(10"^-16~"m)") 
  ))) +
  ylab(expression(atop(atop(
    textstyle("Specific"), 
    textstyle("leaf area,")), 
    textstyle(tilde("SLA")~"(m"^2~"kg"^-1*")") 
  ))) 
# p5


png(file="paramters_from_traits_SI.png", width=650*1.3*3, height=300*1.3*3, res=300)
cowplot::plot_grid(p1, p5, labels=LETTERS, label_size = 14, label_x = 0.1, label_colour = "grey50", hjust = -3, align = "hv", rel_widths = 1, ncol=2)
dev.off()




# p1 = d %>% filter(Species != "Helianthus annuus") %>% ggplot(mapping = aes(y=gamma, x=A.G, color=A.G)) +
#   my_theme()+
#   # scale_x_log10() + 
#   scale_y_log10() +
#   geom_boxplot()
# + 
#   geom_hline(data = d %>% group_by(A.G) %>% summarise(ym=mean(gamma)), mapping = aes(yintercept=ym, color=A.G))+
#   scale_shape_manual(values = c(21, 16)) +
#   scale_color_manual(values = spp_cols)+
#   # geom_smooth(method = "rlm", se = F) + 
#   # xlab(expression(gamma%*%psi[50]^"2"))+
#   # ylab(expression(alpha))
#   ylab(expression(atop(atop(
#     textstyle("Hydraulic"), 
#     textstyle("cost")), 
#     textstyle(gamma) 
#   ))) +
#   xlab(expression(atop(atop(
#     textstyle("Photosynthetic"), 
#     textstyle("cost")), 
#     textstyle(alpha) 
#   ))) 
# 
# p1
# 
# weird <- scales::trans_new("signed_log",
#                            transform=function(x) -log((-x)),
#                            inverse=function(x) -exp(-(x)))
# 
# p2 = d %>% ggplot(mapping = aes(y=K.scalar, x=-P50, color=A.G)) +
#   my_theme()+
#   scale_x_continuous(trans="log", breaks = c(0.5, 1, 2, 4)) + scale_y_log10() + 
#   geom_point(aes(color=A.G, shape = factor(dpsi_cal)), size=2) + 
#   scale_shape_manual(values = c(21, 16)) +
#   scale_color_manual(values = spp_cols)+
#   # geom_smooth(aes(group=A.G), method = lm, se = F) + 
#   geom_smooth(method = lm, se = F, col="black") + 
#   # xlab(expression(-psi[50]))+
#   # ylab(expression(K)) +
#   xlab(expression(atop(atop(
#     textstyle("Plant hydraulic"), 
#     textstyle("vulnerability")), 
#     textstyle(-psi[50]~"(MPa)") 
#   ))) +
#   ylab(expression(atop(atop(
#     textstyle("Plant"), 
#     textstyle("conductivity")), 
#     textstyle(K~"(x10"^"-16"~m*")") 
#   ))) +
#   expand_limits(y=10)
# # p2
# 
# p3 = d %>% ggplot(mapping = aes(y=gamma, x=P50, color=A.G)) +
#   my_theme()+
#   scale_x_continuous(trans=weird, breaks = c(-0.5, -1, -2)) + scale_y_log10() + 
#   geom_point(aes(colour=A.G, shape=factor(dpsi_cal)), size=2) + 
#   geom_hline(data = d %>% group_by(A.G) %>% summarise(ym=mean(gamma)), mapping = aes(yintercept=ym, color=A.G))+
#   # geom_smooth(aes(group=A.G), method = rlm, se = F) + 
#   scale_shape_manual(values = c(21, 16)) +
#   scale_color_manual(values = spp_cols)+
#   # xlab(expression(psi[50]))+
#   # ylab(expression(gamma)) 
#   xlab(expression(atop(atop(
#     textstyle("Plant hydraulic"), 
#     textstyle("vulnerability")), 
#     textstyle(psi[50]~"(MPa)") 
#   ))) +
#   ylab(expression(atop(atop(
#     textstyle("Hydraulic"), 
#     textstyle("cost")), 
#     textstyle(gamma) 
#   ))) 
# p3
# 
# 
# p4 = d %>% ggplot(mapping = aes(x=pgL88, y=P50X)) +
#   my_theme()+
#   geom_point(aes(colour=A.G, shape=factor(dpsi_cal)), size=2) + 
#   geom_smooth(method = rlm, se = F, col="black") + 
#   scale_shape_manual(values = c(21, 16)) +
#   scale_color_manual(values = spp_cols)+
#   geom_abline(slope = 1, intercept=0, color="grey")+
#   # ylab(expression(tilde(psi)["50X"]))+
#   # xlab(expression(psi["g88"])) 
#   xlab(expression(atop(atop(
#     textstyle("Stomatal closure"), 
#     textstyle("point")), 
#     textstyle(psi["g88"]~"(MPa)") 
#   ))) +
#   ylab(expression(atop(atop(
#     textstyle("Xylem hydraulic"), 
#     textstyle("vulnerability")), 
#     textstyle(tilde(psi)["50X"]~"(MPa)") 
#   ))) 
# # geom_abline(slope = 1, intercept = 0, col="grey")
# p4
# 
# 
# p4.1 = d %>% ggplot(mapping = aes(x=P50, y=P50X)) +
#   my_theme()+
#   geom_point(aes(colour=A.G, shape=factor(dpsi_cal)), size=2) + 
#   geom_smooth(method = rlm, se = F, col="black") + 
#   scale_shape_manual(values = c(21, 16)) +
#   scale_color_manual(values = spp_cols)+
#   geom_abline(slope = 1, intercept=0, color="grey")+
#   # ylab(expression(tilde(psi)["50X"]))+
#   # xlab(expression(psi["g88"])) 
#   xlab(expression(atop(atop(
#     textstyle("Plant hydraulic"), 
#     textstyle("vulnerability")), 
#     textstyle(psi["50"]~"(MPa)") 
#   ))) +
#   ylab(expression(atop(atop(
#     textstyle("Xylem hydraulic"), 
#     textstyle("vulnerability")), 
#     textstyle(tilde(psi)["50X"]~"(MPa)") 
#   ))) 
# # geom_abline(slope = 1, intercept = 0, col="grey")
# p4.1
# 
# 
# 
# 
# p6 = d %>% ggplot(mapping = aes(y=pgL88, x=P88)) +
#   my_theme()+
#   geom_point(aes(colour=A.G, shape=factor(dpsi_cal)), size=2) + 
#   geom_smooth(method = rlm, se = F, col="black") + 
#   scale_shape_manual(values = c(21, 16)) +
#   scale_color_manual(values = spp_cols)+
#   geom_abline(slope = 1, intercept=0, color="grey")+
#   # ylab(expression(psi["g12"]))+
#   # xlab(expression(psi["50"])) 
#   xlab(expression(atop(atop(
#     textstyle("Plant hydraulic"), 
#     textstyle("vulnerability")), 
#     textstyle(psi["88"]~"(MPa)") 
#   ))) +
#   ylab(expression(atop(atop(
#     textstyle("Stomatal closure-"), 
#     textstyle("point")), 
#     textstyle(psi["g88"]~"(MPa)") 
#   ))) 
# p6
# 
# 
# p7 = d %>% filter(K.scalar < 10) %>% 
#   ggplot(aes(x=-pgL88, y=K.scalar)) + 
#   my_theme()+
#   geom_point(aes(colour=A.G)) + 
#   geom_smooth(method = 'lm', se=F)
# 
# p7
# 
# p7.1 = d %>% filter(K.scalar < 10) %>% 
#   ggplot(aes(x=-pgL50, y=K.scalar)) + 
#   my_theme()+
#   geom_point(aes(colour=A.G)) + 
#   geom_smooth(method = 'lm', se=F)
# 
# p7.1
# 
# p8.1 = d %>% filter(K.scalar < 10) %>% 
#   ggplot(aes(x=-P50, y=K.scalar)) + 
#   my_theme()+
#   geom_point(aes(colour=A.G)) + 
#   geom_smooth(method = 'lm', se=F)
# 
# p8.1
# 
# p8 = d %>% filter(K.scalar < 10) %>% 
#   ggplot(aes(x=-pgL50, y=P50)) + 
#   my_theme()+
#   geom_point(aes(colour=A.G)) + 
#   geom_smooth(method = 'lm', se=F)
# 
# p8
# 
# 
# 
# # png(filename = "C:/Users/Jaideep/Documents/codes/Phydro_paper/hyd_strategies_fig_4.png", height=550*3*1.5, width=850*3*1.5, res=300)
# # cowplot::plot_grid(p6, p2, p4, p3, p1, p5, labels="AUTO", label_size = 18, label_colour = "grey40", label_x = 0.4, align = "hv", rel_widths = 1)
# # dev.off()
# 
# 
# 
# #### Pairs plot #####
# 
# 
# 
# panel_func = function(x,y,...){
#   points(y~x,...);
#   mod = lm(y~x, weights = rep(1, length(x)) );
#   if (abs(summary(mod)$coefficients[2,3]) > 1.96) {
#     abline(lm(y~x), col="black")
#     mtext(text = sprintf("r = %.2f", cor(x,y, use="pairwise.complete.obs")), side=3, line=0., cex=.45)
#   }
#   else {
#     if (abs(summary(mod)$coefficients[2,3]) > .9){
#       abline(lm(y~x), col="grey65")
#     }
#   }
# }
# 
# 
# d %>% 
#   dplyr::filter(!(Species %in% c("Helianthus annuus"))) %>% 
#   dplyr::select(-Species, -yxp50, -b,
#                 -pgL50, -Pgs90, 
#                 -dpsi_cal, -inst,
#                 -psil_slope,  -A.G, 
#                 -X,
#                 -P88) %>%
#   dplyr::select(-D_exp_std, -D_exp) %>% 
#   pairs(panel= panel_func, pch=20, cex=1.2, col="cyan4") #
# 
# 
# 
# #### Strategies fig ####
# 
# library(ggrepel)
# 
# mytheme = function(){
#   theme_classic()+
#     theme(axis.text=element_text(size=14),
#           axis.title=element_text(size=14),
#           plot.tag.position = "topleft")  
# }
# 
# d %>% filter(Species != "Helianthus annuus") %>%  
#   ggplot(aes(x=K.scalar, y=SLA)) + 
#   mytheme()+
#   geom_point(aes(colour=A.G, shape=factor(dpsi_cal)), size=4) + 
#   geom_text_repel(aes(label=Species), hjust=-0.1, vjust=-0.5, size=3)+
#   scale_shape_manual(values=c(21,16))+
#   # scale_color_gradient2(low = "red", high="blue", mid="grey90", midpoint = 1)+
#   xlab("Higher photosynthetic capacity -->")+
#   ylab("Stronger hydraulics -->")+
#   labs(shape = "Group", color = "Anisohydricity")
# 
# d %>% #filter(Species != "Helianthus annuus") %>%  
#   ggplot(aes(x=P50, y=SLA)) + 
#   mytheme()+
#   geom_point(aes(colour=A.G, shape=factor(dpsi_cal)), size=4) + 
#   geom_text_repel(aes(label=Species), hjust=-0.1, vjust=-0.5, size=3)+
#   scale_shape_manual(values=c(21,16))+
#   # scale_color_gradient2(low = "red", high="blue", mid="grey90", midpoint = 1)+
#   xlab("Higher photosynthetic capacity -->")+
#   ylab("Stronger hydraulics -->")+
#   labs(shape = "Group", color = "Anisohydricity")
# 
# 
# 
# 
# 
# par(mfrow=c(2,2), mar=c(4,4,1,1))
# 
# p1 = d %>% ggplot(aes(y=Ptlp, x=P50X)) +
#   mytheme()+
#   geom_point()+
#   geom_abline(slope=1, col="grey")+
#   ylab(expression(psi[tlp]))+
#   xlab(expression(psi["50X"]))
# 
# p2 = d %>% ggplot(aes(y=Ptlp, x=pgL88)) +
#   mytheme()+
#   geom_point()+
#   geom_abline(slope=1, col="grey")+
#   ylab(expression(psi[tlp]))+
#   xlab(expression(psi["g88"]))
# 
# p3 = d %>% ggplot(aes(y=Ptlp, x=pgL12)) +
#   mytheme()+
#   geom_point()+
#   geom_abline(slope=1, col="grey")+
#   ylab(expression(psi[tlp]))+
#   xlab(expression(psi["g12"]))
# 
# p4 = d %>% ggplot(aes(y=Ptlp, x=pgL50)) +
#   mytheme()+
#   geom_point()+
#   geom_abline(slope=1, col="grey")+
#   ylab(expression(psi[tlp]))+
#   xlab(expression(psi["g50"]))
# 
# 
# png(file="C:/Users/Jaideep/Documents/codes/rpmodel/vignettes_add/tlp_fig.png", width=500*3, height=450*3, res=300)
# cowplot::plot_grid(p1,p2,p3,p4, labels=LETTERS, label_size = 16, label_x = 0.32, label_colour = "grey50", hjust = 0, align = "hv", rel_widths = 1, ncol=2)
# dev.off()
# 
# 
# 
# d = d %>% dplyr::select(-pgL50, -pgL88, -Pgs90,  -Ptlp)
# d = d %>% mutate(SM = (-P50 + pgL12))
# 
# mytheme = function(){
#   theme_classic()+
#     theme(axis.text=element_text(size=0),
#           axis.title=element_text(size=14),
#           plot.tag.position = "topleft")  
# }
# 
# 
# 
# 
# d = g
# rownames(d) = d$Species
# d = d %>% dplyr::select(-Species) 
# d = d %>% dplyr::select(-A.G, -dpsi_intercept) 
# d = d %>% dplyr::select(-pgL50, -Pgs90, -Ptlp, -X) 
# 
# 
# hist(d$b)
# 
# # d = d %>% mutate(SM_choat = P50 - P88)
# # d = sapply(d, FUN = function(x){ifelse(x<0, -x, x)}, simplify = T)
# # d = log(d)
# # d = d %>% mutate(logb = log(b)) %>% dplyr::select(-b)
# # d = d %>% mutate(logY = log(gamma)) %>% dplyr::select(-gamma)
# # d = d %>% mutate(logP50 = log(-P50)) %>% dplyr::select(-P50)
# # d = d %>% mutate(logK = log(K.scalar)) %>% dplyr::select(-K.scalar)
# 
# 
# 
# # par(mfrow=c(2,3), mar=c(4,4,1,1), cex.lab=1.2)
# 
# 
# 
# 
# library(MASS)
# library(tidyverse)
# 
# # rownames(d) = d$Species
# 
# g = d %>% 
#   dplyr::select(-yxp50, -b,
#                 -pgL50, -Pgs90, 
#                 -dpsi_cal, -inst,
#                 -psil_slope,  -A.G, 
#                 -X,
#                 -P88) %>%
#   dplyr::select(-D_exp_std) %>% 
#   dplyr::select(-P50X, -Ptlp, -SLA) %>% 
#   column_to_rownames("Species")
# 
# 
# pp2 = g
# pp1 = d %>% dplyr::select(K.scalar, P50, b, alpha, gamma, yxp50)
# # rownames(pp) = d$Species
# # pp = pp %>% dplyr::select(-Species)
# require(ggfortify)
# # rownames(pp1) = d$Species
# outliers = c("Helianthus annuus")
# outliers = c()
# # pp2_noo = pp2[!rownames(pp2) %in% ,]
# # pp1_noo = pp1[!rownames(pp1) %in% c("Helianthus annuus", "Eucalyptus pilularis"),]
# 
# require(pca3d)
# pp3 = log(abs(pp2))
# pca1 = prcomp(pp2[!(rownames(pp2) %in% outliers),], center = T, scale. = T)  
# pca3d(pca1, fancy=T, show.labels = T,labels.col = "grey60", biplot=T, legend = "bottomleft") #, group = d[!(d$Species %in% outliers),"A.G"])
# snapshotPCA3d(file="C:/Users/Jaideep/Documents/codes/rpmodel/vignettes_add/PCA3D.png", res=300)
# 
# mytheme = function(){
#   theme_classic()+
#     theme(axis.text=element_text(size=0),
#           axis.title=element_text(size=14),
#           plot.tag.position = "topleft")  
# }
# 
# 
# 
# p1 = pca1 %>% 
#   autoplot(data=d %>% filter(!(Species %in% outliers)), label = T, loadings=T, loadings.label=T, label.repel=T, loadings.label.colour = "black", loadings.colour = "black", loadings.label.repel=T, x=1, y=2, colour = "A.G") +
#   mytheme()
# p2 = pca1 %>% 
#   autoplot(data=d %>% filter(!(Species %in% outliers)), label = T, loadings=T, loadings.label=T, label.repel=T, loadings.label.colour = "black", loadings.colour = "black",  loadings.label.repel=T, x=2, y=3, colour = "A.G") +
#   mytheme()
# p3 = pca1 %>% 
#   autoplot(data=d %>% filter(!(Species %in% outliers)), label = T, loadings=T, loadings.label=T, label.repel=T, loadings.label.colour = "black", loadings.colour = "black",  loadings.label.repel=T, x=1, y=3, colour = "A.G") +
#   mytheme()
# 
# png(file="C:/Users/Jaideep/Documents/codes/rpmodel/vignettes_add/pca2D.png", width=1000*3.5, height=650*3.5, res=300)
# cowplot::plot_grid(p1, p2, p3, labels=LETTERS, label_size = 16, label_x = 0.1, label_colour = "grey50", hjust = 0, align = "hv", rel_widths = 1, ncol=2)
# dev.off()
# 
# 
# prcomp(pp3, center = T, scale. = T)  %>% autoplot(data=g, label = T, loadings=T, loadings.label=T, label.repel=T, loadings.label.repel=T, x=2, y=3, colour='A.G')
# 
# 
# 
# 
# 
# #### Dpsi histgrams
# 
# datme = read.csv(file = "hydricity.csv")
# datmm = read.csv(file = "C:/Users/Jaideep/Documents/zotero_data/storage/NXGJVVYW/data.csv")
# 
# par(mfrow=c(2,2), mar=c(4,4,0,0))
# 
# hist(datmm$L..MPa., breaks = seq(-3.5,0,length.out=10), main = "", xlab=expression(psi[l]~"at"~psi[s]*"=0"), col=scales::alpha("skyblue", .5))
# hist(datmm$s..MPa.MPa.1., breaks = seq(0,1.5,length.out=10), main = "", xlab = expression("slope of "*psi[l]*"~"*psi[s]~"relationship"), col=scales::alpha("skyblue",.5))
# hist(-datme$dpsi_intercept, breaks = seq(-3.5,0,length.out=10), main = "", xlab=expression(psi[l]~"at"~psi[s]*"=0"), col="skyblue3")
# hist(1-datme$dpsi_slope, breaks = seq(0,1.5,length.out=10), main = "", xlab=expression("slope of "*psi[l]*"~"*psi[s]~"relationship"), col="skyblue3")
# 
# 
#   
# with(d %>% filter(A.G == "Angiosperm"), plot(log(K.scalar)~I(1/SLA)))
# with(d %>% filter(A.G == "Angiosperm"), abline(lm(log(K.scalar)~I(1/SLA))))
# with(d %>% filter(A.G == "Angiosperm"), (lm(log(K.scalar)~I(1/SLA))))
# 
# par(mfrow=c(2,1), mar=c(4,4,1,1))
# with(d %>% filter(A.G == "Angiosperm"), plot((K.scalar)~I(1/SLA)))
# curve(exp(1.71-8.628*x), add=T)
# 
# with(d %>% filter(A.G == "Angiosperm"), plot(log(K.scalar)~I(1/SLA)))
# curve((1.71-8.628*x), add=T)
# 
# with(d, plot(I(-P50)~I(1/SLA)))
# with(d, (abline(lm(I(-P50)~I(1/SLA)))))
# with(d, (abline(lm(I(-P50)~I(1/SLA)))))
# curve((1.71-8.628*x), add=T)
