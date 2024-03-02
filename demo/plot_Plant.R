
plot_dists = function(folder, title, mtext = F){
  
  up1 = read.delim(paste0(dir, "/", folder, "/species_0_u.txt"), header=F)
  up2 = read.delim(paste0(dir, "/", folder, "/species_1_u.txt"), header=F)
  up3 = read.delim(paste0(dir, "/", folder, "/species_2_u.txt"), header=F)
  
  up1 = as.matrix(up1[,-ncol(up1)])
  up2 = as.matrix(up2[,-ncol(up2)])
  up3 = as.matrix(up3[,-ncol(up3)])
  up1[up1<0]=0
  up2[up2<0]=0
  up3[up3<0]=0
  
  hp1 = read.delim(paste0(dir, "/", folder, "/species_0_X.txt"), header=F)
  hp2 = read.delim(paste0(dir, "/", folder, "/species_1_X.txt"), header=F)
  hp3 = read.delim(paste0(dir, "/", folder, "/species_2_X.txt"), header=F)
  
  hp1 = as.matrix(hp1[,-ncol(hp1)])
  hp2 = as.matrix(hp2[,-ncol(hp2)])
  hp3 = as.matrix(hp3[,-ncol(hp3)])
  
  times = as.numeric(hp1[,1])

  col1 = scales::colour_ramp(colors = c('white','#f1eef6','#a6bddb','#2b8cbe','#045a8d', "cyan", "white"))(seq(0,1,0.01))  
  col2 = scales::viridis_pal()(100)
  image(x=times, y=hp1[1,-1],  z=log(1+100*up1[,-1]), col=col1, ylab="Height", xlab="Time (years)")
  mtext(title, side=3, line=1)
  if(mtext) mtext("lma = 0.0825", side=2, line=5)
  image(x=times, y=hp2[1,-1],  z=log(1+100*up2[,-1]), col=col1, ylab="Height", xlab="Time (years)")
  if(mtext) mtext("lma = 0.2625", side=2, line=5)
  image(x=times, y=hp3[1,-1],  z=log(1+100*up3[,-1]), col=col1, ylab="Height", xlab="Time (years)")
  if(mtext) mtext("lma = 0.4625", side=2, line=5)
}

##### Fixed input seed rain mode
dir = "~/codes/libpspm/demo/Plant_model"
setwd(dir)

png("../size_dists_noFeedback.png", width = 1400*4, height=750*3, res=300)
par(mar=c(4,4,1,1), oma = c(1,4,2,1), cex.lab=1.2, cex.axis=1.2)
layout(mat = matrix(1:24, nrow=3, byrow=F))

plot_dists("outputs/fmu", "FMU", T)
plot_dists("outputs/ifmu", "IFMU")
# plot_dists("outputs/ifmu2", "ILUD")
plot_dists("outputs/ebt", "EBT")
plot_dists("outputs/iebt", "IEBT")
plot_dists("outputs/cm", "CM")
plot_dists("outputs/icm", "ICM")
plot_dists("outputs/abm", "ABM")
dev.off()



setwd(dir)
seeds_fmu = read.delim("outputs/fmu/seed_rains.txt", header = F)
seeds_ifmu = read.delim("outputs/ifmu/seed_rains.txt", header = F)
# seeds_ifmu2 = read.delim("outputs/ifmu2/seed_rains.txt", header = F)
seeds_ebt = read.delim("outputs/ebt/seed_rains.txt", header = F)
seeds_iebt = read.delim("outputs/iebt/seed_rains.txt", header = F)
seeds_cm = read.delim("outputs/cm/seed_rains.txt", header = F)
seeds_icm = read.delim("outputs/icm/seed_rains.txt", header = F)
seeds_abm = read.delim("outputs/abm/seed_rains.txt", header = F)

# cols_m = c("purple", "green3", "mediumspringgreen", "darkgoldenrod2", "red3", "pink", "#2b8cbe")
cols_m = c("darkgreen", 
           "yellowgreen",
           # "green3", 
           "magenta", 
           "purple", 
           "darkgoldenrod2", 
           "darkgoldenrod3", 
           "turquoise2")
names = c("FMU", 
          "IFMU",
          # "ILUD", 
          "EBT", 
          "IEBT", 
          "CM", 
          "ICM", 
          "ABM")

plot_seeds = function(y, title, ...){
  matplot(y = y, x=seeds_fmu$V1, type="l", lty=1, col=scales::alpha(cols_m, alpha=0.7), ylab="Seed rain", ...)
  mtext(title, line=1)
}

cairo_pdf("../seed_rains_noFeedback.pdf", width = 6.6, height=7.66)
par(mfrow=c(3,1), mar = c(4,4,1,1), oma = c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
plot_seeds(cbind(seeds_fmu$V2, 
                 seeds_ifmu$V2,
                 # seeds_ifmu2$V2, 
                 seeds_ebt$V2, 
                 seeds_iebt$V2, 
                 seeds_cm$V2, 
                 seeds_icm$V2, 
                 seeds_abm$V2), 
           "Species 1", xlab="", lwd=2, ylim=c(00,400))
legend(x = 60, y=400, legend = names[1:4], col=cols_m[1:4], lwd=2, bty = "n", cex=1.3)
legend(x = 80, y=400, legend = names[5:8], col=cols_m[5:8], lwd=2, bty = "n", cex=1.3)
plot_seeds(cbind(seeds_fmu$V3, 
                 seeds_ifmu$V3,
                 # seeds_ifmu2$V3, 
                 seeds_ebt$V3, 
                 seeds_iebt$V3, 
                 seeds_cm$V3, 
                 seeds_icm$V3, 
                 seeds_abm$V3), 
           "Species 2", xlab="", lwd=2, ylim=c(00,1000))
plot_seeds(cbind(seeds_fmu$V4, 
                 seeds_ifmu$V4,
                 # seeds_ifmu2$V4, 
                 seeds_ebt$V4, 
                 seeds_iebt$V4, 
                 seeds_cm$V4, 
                 seeds_icm$V4, 
                 seeds_abm$V4), 
           "Species 3", xlab="Time (years)", lwd=2, ylim=c(00,1200))
dev.off()




##### Feedback mode

# 
# dir = "~/codes/libpspm/demo/Plant_model"
# # dir = "C:/Users/Jaideep/OneDrive - IIASA/libpspm_paper/demo/Plant_model/"
# setwd(dir)

# png("../size_dists_withFeedback.png", width = 1200*3, height=750*3, res=300)
# par(mar=c(4,4,1,1), oma = c(1,4,2,1), cex.lab=1.2, cex.axis=1.2)
# layout(mat = matrix(1:18, nrow=3, byrow=F))
#   
# plot_dists("outputs/fmu_f/", "FMU", T)
# plot_dists("outputs/ifmu_f", "IFMU")
# plot_dists("outputs/ifmu2_f", "ILUD")
# plot_dists("outputs/ebt_f/", "EBT")
# plot_dists("outputs/iebt_f/", "IEBT")
# plot_dists("outputs/cm_f/", "CM")
# dev.off()
# 
# # 
# setwd(dir)
# seeds_fmu = read.delim("outputs/fmu_f/seed_rains.txt", header = F)
# seeds_ifmu = read.delim("outputs/ifmu_f/seed_rains.txt", header = F)
# seeds_ifmu2 = read.delim("outputs/ifmu2_f/seed_rains.txt", header = F)
# seeds_ebt = read.delim("outputs/ebt_f/seed_rains.txt", header = F)
# seeds_iebt = read.delim("outputs/iebt_f/seed_rains.txt", header = F)
# seeds_cm = read.delim("outputs/cm_f/seed_rains.txt", header = F)
# # seeds_abm = read.delim("outputs/abm_f/seed_rains.txt", header = F)
#   
# cols_m1 = cols_m[-c(3,length(cols_m))]
# 
# plot_seeds = function(y, title, ...){
#   matplot(y = y, x=seeds_fmu$V1, type="l", lty=1, col=scales::alpha(cols_m1, alpha=0.7), ylab="Seed rain", ...)
#   mtext(title, line=1)
# }
# 
# cairo_pdf("../seed_rains_withFeedback.pdf", width = 6.6, height=7.66)
# par(mfrow=c(3,1), mar = c(4,4,1,1), oma = c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
# plot_seeds(cbind(seeds_fmu$V2, seeds_ifmu$V2, seeds_ebt$V2, seeds_iebt$V2, seeds_cm$V2), "Species 1", xlab="", lwd=2)
# legend(x = 150, y=65, legend = c("FMU", "IFMU", "EBT", "IEBT", "CM"), col=cols_m1, lwd=2, bty = "n", cex=1.3)
# plot_seeds(cbind(seeds_fmu$V3, seeds_ifmu$V3, seeds_ebt$V3, seeds_iebt$V3, seeds_cm$V3), "Species 2", xlab="", lwd=2)
# plot_seeds(cbind(seeds_fmu$V4, seeds_ifmu$V4, seeds_ebt$V4, seeds_iebt$V4, seeds_cm$V4), "Species 3", xlab="Time (years)", lwd=2)
# dev.off()


# par(mfrow=c(3,1), mar = c(4,4,1,1), oma = c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
# plot_seeds(cbind(seeds_fmu$V2, seeds_ebt$V2), "Species 1", xlab="", log="")
# plot_seeds(cbind(seeds_fmu$V3, seeds_ebt$V3), "Species 2", xlab="", log="")
# plot_seeds(cbind(seeds_fmu$V4, seeds_ebt$V4), "Species 3", xlab="Time (years)", log="")
# 
# 
# par(mfrow=c(3,1), mar = c(4,4,1,1), oma = c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
# plot_seeds(cbind(seeds_fmu$V2), "Species 1", xlab="")
# plot_seeds(cbind(seeds_fmu$V3), "Species 2", xlab="")
# plot_seeds(cbind(seeds_fmu$V4), "Species 3", xlab="Time (years)")
# 
# # png("../seed_rains.png", width = 660*3, height=766*3, res=300)
# par(mfrow=c(3,1), mar = c(4,4,1,1), oma = c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
# plot_seeds(cbind(seeds_fmu$V2, seeds_ebt$V2, seeds_cm$V2), "Species 1", xlab="")
# plot_seeds(cbind(seeds_fmu$V3, seeds_ebt$V3, seeds_cm$V3), "Species 2", xlab="")
# plot_seeds(cbind(seeds_fmu$V4, seeds_ebt$V4, seeds_cm$V4), "Species 3", xlab="Time (years)")
# # dev.off()




setwd(dir)
seeds_fmu = read.delim("outputs/fmu_f3/seed_rains.txt", header = F)
seeds_ifmu = read.delim("outputs/ifmu_f3/seed_rains.txt", header = F)
# seeds_ifmu2 = read.delim("outputs/ifmu2_f/seed_rains.txt", header = F)
seeds_ebt = read.delim("outputs/ebt_f3/seed_rains.txt", header = F)
seeds_iebt = read.delim("outputs/iebt_f3/seed_rains.txt", header = F)
seeds_cm = read.delim("outputs/cm_f3/seed_rains.txt", header = F)
# seeds_icm = read.delim("outputs/icm_f3/seed_rains.txt", header = F)
# seeds_abm = read.delim("outputs/abm_f3/seed_rains.txt", header = F)

cols_m1 = cols_m[c(1:5, 7)]
names1 = names[c(1:5, 7)]

plot_seeds = function(y, title, ...){
  matplot(y = y, x=seeds_fmu$V1, type="l", lty=1, col=scales::alpha(cols_m1, alpha=0.7), ylab="Seed rain", ...)
  mtext(title, line=1)
}

cairo_pdf("../seed_rains_withFeedback_t400.pdf", width = 6.6, height=7.66)
par(mfrow=c(3,1), mar = c(4,4,1,1), oma = c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
plot_seeds(cbind(seeds_fmu$V2, 
                 seeds_ifmu$V2,
                 seeds_ebt$V2, 
                 seeds_iebt$V2, 
                 seeds_cm$V2
                 # seeds_abm$V2
                 ), 
           "Species 1", xlab="", lwd=2)
legend(x = 320, y=60, legend = names1, col=cols_m1, lwd=2, bty = "n", cex=1.3)
plot_seeds(cbind(seeds_fmu$V3, 
                 seeds_ifmu$V3,
                 seeds_ebt$V3, 
                 seeds_iebt$V3, 
                 seeds_cm$V3
                 # seeds_abm$V3
                 ), 
           "Species 2", xlab="", lwd=2)
plot_seeds(cbind(seeds_fmu$V4, 
                 seeds_ifmu$V4,
                 seeds_ebt$V4, 
                 seeds_iebt$V4, 
                 seeds_cm$V4
                 # seeds_abm$V4
                 ), 
           "Species 3", xlab="Time (years)", lwd=2)
dev.off()



dir = "~/codes/libpspm/demo/Plant_model"

png("../size_dists_withFeedback_u0.png", width = 600/1.2*3, height=750/1.2*3, res=300)
par(mar=c(4,4,1,1), oma = c(1,4,2,1), cex.lab=1.2, cex.axis=1.2)
layout(mat = matrix(1:6, nrow=3, byrow=F))

plot_dists("outputs/ifmu_f3_ic", "IFMU")
plot_dists("outputs/iebt_f3_ic", "IEBT")
dev.off()

setwd(dir)
seeds_ifmu = read.delim("outputs/ifmu_f3_ic/seed_rains.txt", header = F)
# seeds_ifmu2 = read.delim("outputs/ifmu2_f/seed_rains.txt", header = F)
seeds_iebt = read.delim("outputs/iebt_f3_ic/seed_rains.txt", header = F)
# seeds_abm = read.delim("outputs/abm_f3/seed_rains.txt", header = F)

cols_m2 = cols_m[c(2,4)]

plot_seeds = function(y, title, ...){
  matplot(y = y, x=seeds_iebt$V1, type="l", lty=1, col=scales::alpha(cols_m2, alpha=0.7), ylab="Seed rain", ...)
  mtext(title, line=1)
}

cairo_pdf("../seed_rains_withFeedback_u0.pdf", width = 6.6, height=7.66)
par(mfrow=c(3,1), mar = c(4,4,1,1), oma = c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
plot_seeds(cbind(
  seeds_ifmu$V2,
  seeds_iebt$V2
  ), "Species 1", xlab="", lwd=2, xlim=c(0,200))
plot_seeds(cbind(
  seeds_ifmu$V3,
  seeds_iebt$V3
  ), "Species 2", xlab="", lwd=2, xlim=c(0,200))
legend(x = 160, y=275, legend = c("IFMU", "IEBT"), col=cols_m2, lwd=2, bty = "n", cex=1.3)
plot_seeds(cbind(
  seeds_ifmu$V4,
  seeds_iebt$V4), "Species 3", xlab="Time (years)", lwd=2, xlim=c(0,200))
dev.off()



my = c(1.3653132665544065e-216, 1.3816055871431628e-216, 1.6681228668550405e-216,         3.5287923779251789e-216, 2.2090043248436949e-215, 6.8864724381916476e-214, 1.6652782834571183e-211,                             4.3411376729180807e-208, 1.4815511691232643e-203, 7.0509524976408433e-198, 4.4324332226239083e-191,                             3.1774966411448391e-183, 2.0760371873055674e-174, 9.3593164345929964e-165, 2.1294696170414097e-154,                             1.7640096088421506e-143, 3.8172111284339624e-132, 1.5582032270059151e-120, 8.782550672582997e-109,                              5.1039164383764004e-97, 2.3211761073295252e-85, 1.3153690385356825e-79, 6.3610550205392561e-74, 2.5436575657998198e-68,         8.1512264975096469e-63, 2.0275240781180787e-57, 3.7865032676920237e-52, 5.1280508290828182e-47, 4.8538412134174207e-42,         3.0855952472940975e-37, 1.2616197239893738e-32, 2.1300015399314637e-30, 3.1646753082270224e-28, 4.1070202015112156e-26,     4.6254213940286272e-24, 4.4886679872028685e-22, 3.7233000133439702e-20, 2.621536511236298e-18, 1.5526082430102352e-16,          7.67925275558447e-15, 3.1422280412950867e-13, 1.0554320818787243e-11, 2.8844417396961996e-10, 6.3598167239283952e-09,           1.1226552578108941e-07, 1.573527013335508e-06, 1.7406200281955848e-05, 5.2857392251496907e-05, 0.00015095613454644836,          0.00040531793927071642, 0.0010227556268993675, 0.0024243304057254404, 0.005397137699403501, 0.011288613279157828,               0.022189538215720291, 0.041012675584603635, 0.071305984324963545, 0.11675521479503821, 0.18031056027813858,                     0.21933541788913755, 0.26306881603935967, 0.3111917862623601, 0.36329907338504513, 0.47597389857715683,                         0.59345969018859446, 0.70662916076972893, 0.80671768176832093, 0.82870123755352465, 0.84943819321263148,                        0.86883022537606591, 0.88685750123206419, 0.91824944178368395, 0.94384895254129109, 0.95460453238172061,                        0.96393692449885526, 0.97179706594746118, 0.97848822720319861, 0.98866945669034745, 0.99486069639423769,                        0.99676559362338035, 0.99811678057512432, 0.99907407073534371, 0.99967420374378846, 0.99998340031355171, 1)
mx = c(0, 0.625, 1.25, 1.875, 2.5, 3.125, 3.75, 4.375, 5, 5.625, 6.25, 6.875, 7.5,        8.125, 8.75, 9.375, 10, 10.625, 11.25, 11.875, 12.5, 12.8125, 13.125, 13.4375, 13.75, 14.0625, 14.375, 14.6875, 15,             15.3125, 15.625, 15.78125, 15.9375, 16.09375, 16.25, 16.40625, 16.5625, 16.71875, 16.875, 17.03125, 17.1875, 17.34375,          17.5, 17.65625, 17.8125, 17.96875, 18.125, 18.203125, 18.28125, 18.359375, 18.4375, 18.515625, 18.59375, 18.671875,             18.75, 18.828125, 18.90625, 18.984375, 19.0625, 19.1015625, 19.140625, 19.1796875, 19.21875, 19.296875, 19.375,                 19.453125, 19.53125, 19.55078125, 19.5703125, 19.58984375, 19.609375, 19.6484375, 19.6875, 19.70703125, 19.7265625,             19.74609375, 19.765625, 19.8046875, 19.84375, 19.86328125, 19.8828125, 19.90234375, 19.921875, 19.9609375, 20)

