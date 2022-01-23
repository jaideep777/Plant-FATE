setwd("/home/jaideep/codes/tmodel_cpp/pspm_output/trait_100_trial_4_outputTestNew/")

plot_to_file = F
n_species = 5
n = 101

hp   = read.delim(paste0("species_",0,"_height.txt"), header=F, col.names = paste0("V", 1:n))
Dp   = read.delim(paste0("species_",0,"_X.txt"), header=F, col.names = paste0("V", 1:n))
Zp   = read.delim("z_star.txt", header=F, col.names = paste0("V", 1:50))
COp   = read.delim("canopy_openness.txt", header=F, col.names = paste0("V", 1:50))
BAp   = read.delim("basal_area.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
#LAIp   = read.delim("LAI.txt", header=F, col.names = paste0("V", 1:10))
seeds = read.delim("seeds.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
#cwmt = read.delim("cwmt.txt", header=F, col.names = paste0("V", 1:20))

amzY = read.delim("AmzFACE_Y.txt", header=T)
amzD = read.delim("AmzFACE_D.txt", header=T)
tD = amzD$YEAR + amzD$DOY/365
tY = amzY$YEAR

ids_x = 2:31
ids_h = which(diff(as.numeric(hp[1,ids_x])) > 0)[-1]
t = Dp[,1]

Up   = read.delim(paste0("species_",0,"_u.txt"), header=F, col.names = paste0("V", 1:n))
Up_all = Up*0
for (i in 0:(n_species-1)){
  Up   = read.delim(paste0("species_",i,"_u.txt"), header=F, col.names = paste0("V", 1:n))
  ids_x = 2:31
  ids_h = which(diff(as.numeric(hp[1,ids_x])) > 0)[-1]
  Up_all[,-1] = Up_all[,-1] + Up[,-1]
}

if (plot_to_file) png("emergent_props.png", height=750*3, width=630*3, res=300)
par(mfrow=c(4,2), mar=c(4,5,1,.5))

cols = rainbow(n = length(ids_x)+5, start = 0, end = 0.85)
cols_t = scales::alpha(cols, 0.5)


matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")

matplot(COp$V1, COp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Canopy openness")

matplot(BAp$V1, BAp[,-1]*1e4, lty=1, col=rainbow(n = n_species, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Basal area (m2 / Ha)\n", log="")#, ylim=c(exp(-9), exp(3.5)))

matplot(tD, amzD$LAI, lty=1, col=rainbow(n = n_species, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Community LAI")
  
matplot(seeds$V1, seeds[,-1], lty=1, col=rainbow(n = n_species, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Species Seed output", log="")

plot(tY, amzY$BA*1e4, lty=1, type="l",
     las=1, xlab="Time (years)", ylab="Total BA", log="")


matplot(x=(t(Dp[,-1])), y=(t(Up_all[,-1]))*1e4/100, lty=1, col=scales::alpha(rainbow(n = nrow((Dp)), start = 0, end = 0.75),10/length(t)), type="l",
        las=1, xlab = "Diameter", ylab="Density (Ind/cm/Ha)\n", log="xy", ylim=c(exp(-10), exp(10)))
Up_mean = colMeans(Up_all[t>1000, ids_x]*1e4/100)
points(as.numeric(Up_mean)~as.numeric(Dp[1,ids_x]))

matplot(x=(t(hp[,-1])), y=(t(Up_all[,-1]))*1e4/100, lty=1, col=scales::alpha(rainbow(n = nrow((Dp)), start = 0, end = 0.75),10/length(t)), type="l",
        las=1, xlab = "Height", ylab="Density (Ind/cm/Ha)\n", log="xy", ylim=c(exp(-10), exp(10)))
Up_mean = colMeans(Up_all[t>1000, ids_h]*1e4/100)
points(as.numeric(Up_mean)~as.numeric(hp[1,ids_h]))

if (plot_to_file) dev.off()


if (plot_to_file) png("CWM_traits.png", height=680*3, width=680*3, res=300)

par(mfrow=c(4,2), mar=c(4,5,1,.5))
matplot(tD, cbind(amzD$GPP, amzD$NPP, amzD$RAU), lty=1, col=c("green4", "green3", "brown"), type="l",
        las=1, xlab="Time (years)", ylab="GPP, NPP, Auto.Resp (gC/m2/d)", log="")

matplot(tD, amzD$ET, lty=1, col="blue", type="l",
        las=1, xlab="Time (years)", ylab="Transpiration (mm/d)", log="")

matplot(tD, cbind(amzD$CL, amzD$CW, amzD$CCR, amzD$CFR), lty=1, col=c("green4", "brown", "yellow3", "yellow2"), type="l",
        las=1, xlab="Time (years)", ylab="Biomass pools (gC/m2)", log="")

matplot(tY, amzY$DE*1e4, lty=1, col=rainbow(n = n_species, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Number of individuals / Ha", log="")

matplot(tY, amzY$TB, lty=1, col=rainbow(n = n_species, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Biomass (kg / m2)", log="")

matplot(tY, amzY$SLA*1000, lty=1, col=rainbow(n = n_species, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="CWM LMA", log="")

matplot(tY, amzY$WD, lty=1, col=rainbow(n = n_species, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="CWM Wood density", log="")

matplot(tY, amzY$P50, lty=1, col=rainbow(n = n_species, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="CWM Xylem P50", log="")


if (plot_to_file) dev.off()
# 
# par(mfrow = c(4,3), mar=c(4,5,1,1), oma=c(1,1,1,1))
# for (i in 0:9){
#   Dp   = read.delim(paste0("species_",i,"_X.txt"), header=F, col.names = paste0("V", 1:n))
#   Up   = read.delim(paste0("species_",i,"_u.txt"), header=F, col.names = paste0("V", 1:n))
# 
#   matplot(x=(t(Dp[,-1])), y=(t(Up[,-1]))*1e4/100, lty=1, col=scales::alpha(rainbow(n = nrow((Dp)), start = 0, end = 0.75),0.05), type="l",
#         las=1, xlab = "Diameter", ylab="Density (Ind/cm/Ha)\n", log="xy", ylim=c(exp(-15), exp(10)))
#   Up_mean = colMeans(Up[t>500, ids_x]*1e4/100)
#   points(as.numeric(Up_mean)~as.numeric(Dp[1,ids_x]))
# }
# 
# par(mfrow = c(4,3), mar=c(4,5,1,1), oma=c(1,1,1,1))
# for (i in 0:9){
#   hp   = read.delim(paste0("species_",i,"_height.txt"), header=F, col.names = paste0("V", 1:n))
#   Up   = read.delim(paste0("species_",i,"_u.txt"), header=F, col.names = paste0("V", 1:n))
# 
#   matplot(x=(t(hp[,-1])), y=(t(Up[,-1]))*1e4/100, lty=1, col=scales::alpha(rainbow(n = nrow((hp)), start = 0, end = 0.75),0.05), type="l",
#           las=1, xlab = "Height", ylab="Density (Ind/cm/Ha)\n", log="y", ylim=c(exp(-10),exp(10)))
# }
# 
# 
par(mfrow = c(5,6), mar=c(4,5,.5,.5), oma=c(1,1,1,1))
for (i in 0:(n_species-1)){
  hp   = read.delim(paste0("species_",i,"_height.txt"), header=F, col.names = paste0("V", 1:n))
  Dp   = read.delim(paste0("species_",i,"_X.txt"), header=F, col.names = paste0("V", 1:n))
  Up   = read.delim(paste0("species_",i,"_u.txt"), header=F, col.names = paste0("V", 1:n))
  Lp   = read.delim(paste0("species_",i,"_lai.txt"), header=F, col.names = paste0("V", 1:n))
  Gp   = read.delim(paste0("species_",i,"_g.txt"), header=F, col.names = paste0("V", 1:n))
  Ap   = read.delim(paste0("species_",i,"_gpp.txt"), header=F, col.names = paste0("V", 1:n))
  ids_x = 2:31
  ids_h = which(diff(as.numeric(hp[1,ids_x])) > 0)[-1]

  image(x=hp[,1], y=as.numeric(hp[1,ids_h]), z=log(1+1e8*pmax(as.matrix(Up[,ids_h])/max(as.matrix(Up[,ids_h]),na.rm=T),0)), log="", xlab="Time", ylab = "Height distribution, Z*", col=scales::viridis_pal()(100))
  matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
          las=1, xlab="Time (years)", ylab="Z*", add=T)

  image(x=Dp[,1], y=as.numeric(Dp[1,ids_x]), z=log(1+1e8*pmax(as.matrix(Up[,ids_x])/max(as.matrix(Up[,ids_x]),na.rm=T),0)), log="y", xlab="Time", ylab = "Diameter distribution", col=scales::viridis_pal()(100))

  image(x=hp[,1], y=as.numeric(hp[1,ids_h]), z=as.matrix(Lp[,ids_h]), zlim=c(0,8), log="", xlab="Time", ylab = "LAI distribution, Z*", col=scales::viridis_pal()(100))
  matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
          las=1, xlab="Time (years)", ylab="Z*", add=T)

  image(x=hp[,1], y=as.numeric(hp[1,ids_h]), z=log(1+1e4*pmax(as.matrix(Gp[,ids_h])/max(as.matrix(Gp[,ids_h]),na.rm=T),0)), log="", xlab="Time", ylab = "R. Growth rate dist, Z*", col=scales::viridis_pal()(100))
  matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
          las=1, xlab="Time (years)", ylab="Z*", add=T)

  image(x=hp[,1], y=as.numeric(hp[1,ids_h]), z=log(1+1e4*pmax(as.matrix(Ap[,ids_h])/max(as.matrix(Ap[,ids_h]),na.rm=T),0)), log="", xlab="Time", ylab = "GPP dist, Z*", col=scales::viridis_pal()(100))
  matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
          las=1, xlab="Time (years)", ylab="Z*", add=T)

  matplot(x=(t(Dp[,-1])), y=(t(Up[,-1]))*1e4/100, lty=1, col=scales::alpha(rainbow(n = nrow((Dp)), start = 0, end = 0.75),0.05), type="l",
        las=1, xlab = "Diameter", ylab="Density (Ind/cm/Ha)\n", log="xy", ylim=c(exp(-15), exp(10)))
  Up_mean = colMeans(Up[t>100, ids_x]*1e4/100)
  points(as.numeric(Up_mean)~as.numeric(Dp[1,ids_x]))

}




#### TRAIT SPACE ####
trait = read.csv("../../tests/data/trait_100_filled.csv", header = TRUE)

trait$Leaf.LMA..g.m2.[is.na(trait$Leaf.LMA..g.m2.)] <- 119.37
trait$meanWoodDensity..g.cm3.[is.na(trait$meanWoodDensity..g.cm3.)] <- 0.68
trait$Height_Max.m.[is.na(trait$Height_Max.m.)] <- 23.99

dominant = which(BAp[nrow(BAp),-1]*1e4 >= 0.01)
length(dominant)
T1 <- trait[dominant,]

BAdominant <- as.numeric(BAp[nrow(BAp),-1])[dominant]*1e4
BAsimulated <- as.numeric(BAp[nrow(BAp),-1])[1:100]*1e4

t50 <- trait[1:100,]
p <- t50 %>% ggplot(aes(x=Leaf.LMA..g.m2., y=meanWoodDensity..g.cm3., colour = log(BA))) + geom_point() + theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1)) + geom_point(size = 3) + xlab("LMA") + ylab("WD")
p1 <- p + scale_colour_viridis_c(limits = log(c(0.01, max(t50$BA))))

p <- t50 %>% ggplot(aes(x=Leaf.LMA..g.m2., y=Height_Max.m., colour = log(BA))) + geom_point() + theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1)) + geom_point(size = 3) + xlab("LMA") + ylab("HT")
p2 <- p + scale_colour_viridis_c(limits = log(c(0.01, max(t50$BA))))

p <- t50 %>% ggplot(aes(x=Leaf.LMA..g.m2., y=P50..Mpa., colour = log(BA))) + geom_point() + theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1)) + geom_point(size = 3) + xlab("LMA") + ylab("HT")
p3 <- p + scale_colour_viridis_c(limits = log(c(0.01, max(t50$BA))))

q <- t50 %>% ggplot(aes(x=Leaf.LMA..g.m2., y=meanWoodDensity..g.cm3., colour = log(BAsimulated))) + geom_point() + theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1)) + geom_point(size = 3) + xlab("LMA") + ylab("WD") 
q1 <- q + scale_colour_viridis_c(limits = log(c(0.01, max(BAsimulated))))

q <- t50 %>% ggplot(aes(x=Leaf.LMA..g.m2., y=Height_Max.m., colour = log(BAsimulated))) + geom_point() + theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1)) + geom_point(size = 3) + xlab("LMA") + ylab("HT")
q2 <- q + scale_colour_viridis_c(limits = log(c(0.01, max(BAsimulated))))

q <- t50 %>% ggplot(aes(x=Leaf.LMA..g.m2., y=P50..Mpa., colour = log(BAsimulated))) + geom_point() + theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1)) + geom_point(size = 3) + xlab("LMA") + ylab("P50")
q3 <- q + scale_colour_viridis_c(limits = log(c(0.01, max(BAsimulated))))

if (plot_to_file) png("trait_space.png", height=500*3, width=1450*3, res=300)
cowplot::plot_grid(p1, p2, p3, q1, q2, q3, nrow = 2, align = "hv")
if (plot_to_file) dev.off()

  
plot(T1$BA~BAdominant, log="xy")
abline(0,1)
plot(t50$BA~BAsimulated, log="xy")
abline(0,1)


