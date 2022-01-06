setwd("/home/jaideep/codes/tmodel_cpp/pspm_output/amazontraitreduceBA6trait")

n_species = 2
n = 201

hp   = read.delim(paste0("species_",0,"_height.txt"), header=F, col.names = paste0("V", 1:n))
Dp   = read.delim(paste0("species_",0,"_X.txt"), header=F, col.names = paste0("V", 1:n))
Zp   = read.delim("z_star.txt", header=F, col.names = paste0("V", 1:50))
COp   = read.delim("canopy_openness.txt", header=F, col.names = paste0("V", 1:50))
BAp   = read.delim("basal_area.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
LAIp   = read.delim("LAI.txt", header=F, col.names = paste0("V", 1:10))
seeds = read.delim("seeds.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
  
ids_x = 2:31
ids_h = which(diff(as.numeric(hp[1,ids_x])) > 0)[-1]
t = Dp[,1]


par(mfrow=c(3,2), mar=c(4,5,1,.5))

cols = rainbow(n = length(ids_x)+5, start = 0, end = 0.85)
cols_t = scales::alpha(cols, 0.5)


matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")

matplot(COp$V1, COp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Canopy openness")

matplot(BAp$V1, BAp[,-1]*1e4, lty=1, col=rainbow(n = n_species, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Basal area (m2 / Ha)\n", log="")#, ylim=c(exp(-9), exp(4)))

matplot(LAIp$V1, LAIp[,-1], lty=1, col=rainbow(n = n_species, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Community LAI")

matplot(seeds$V1, seeds[,-1], lty=1, col=rainbow(n = n_species, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Species Seed output", log="")

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

  

