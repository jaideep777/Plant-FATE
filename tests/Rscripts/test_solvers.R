output_dir = "pspm_output_3"
prefix = "timestep_test"

solver = "IEBT0.1"
setwd(paste0("~/codes/Plant-FATE/",output_dir,"/",prefix,"_",solver))

plot_to_file = F
n_species = 100
n = 101

seeds = read.delim("seeds.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
Zp   = read.delim("z_star.txt", header=F, col.names = paste0("V", 1:50))
BA  =     read.delim("basal_area.txt", header=F, col.names = paste0("V", 1:(n_species+2)))

par(mfcol=c(2,3), mar=c(6,6,1,1), oma=c(1,1,2,1), cex.lab=1.3, cex.axis=1.2, mgp=c(4,1,0))
matplot(seeds$V1, seeds[,-1], lty=1, col=rainbow(n = n_species+1, start = 0, end = 0.65), type="l",
        las=1, xlab="Time (years)", ylab="Species Seed output", log="")
mtext(line=0.5, side=3, text=solver)
matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")

matplot(BA$V1, cbind(BA[,-1], rowSums(BA[,-1], na.rm=T))*1e4, lty=1, col=c(rainbow(n = n_species+1, start = 0, end = 0.85), "black"), type="l",
        las=1, xlab="Time (years)", ylab="Basal area", log="")


solver = "ABM0.1"
setwd(paste0("~/codes/Plant-FATE/",output_dir,"/",prefix,"_",solver))

plot_to_file = F
n_species = 3
n = 102

seeds = read.delim("seeds.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
Zp   = read.delim("z_star.txt", header=F, col.names = paste0("V", 1:50))


#par(mfrow=c(2,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
matplot(seeds$V1, seeds[,-1], lty=1, col=rainbow(n = n_species, start = 0, end = 0.65), type="l",
        las=1, xlab="Time (years)", ylab="Species Seed output", log="")
mtext(line=0.5, side=3, text=solver)
matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")
    


solver = "ABM0.1_20"
setwd(paste0("~/codes/Plant-FATE/",output_dir,"/",prefix,"_",solver))

plot_to_file = F
n_species = 3
n = 102

seeds = read.delim("seeds.txt", header=F, col.names = paste0("V", 1:(n_species+2)))
Zp   = read.delim("z_star.txt", header=F, col.names = paste0("V", 1:50))


#par(mfrow=c(2,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
matplot(seeds$V1, seeds[,-1], lty=1, col=rainbow(n = n_species, start = 0, end = 0.65), type="l",
        las=1, xlab="Time (years)", ylab="Species Seed output", log="")
mtext(line=0.5, side=3, text=solver)
matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")


# 
# 
# # 
# # par(mfrow = c(5,6), mar=c(4,5,.5,.5), oma=c(1,1,1,1))
# for (i in 0:(n_species-1)){
# #   # hp   = read.delim(paste0("species_",i,"_height.txt"), header=F, col.names = paste0("V", 1:n))
#   Dp   = read.delim(paste0("species_",i,"_X.txt"), header=F, col.names = paste0("V", 1:n))
#   Up   = read.delim(paste0("species_",i,"_u.txt"), header=F, col.names = paste0("V", 1:n))
# #   # Lp   = read.delim(paste0("species_",i,"_lai.txt"), header=F, col.names = paste0("V", 1:n))
# #   # Gp   = read.delim(paste0("species_",i,"_g.txt"), header=F, col.names = paste0("V", 1:n))
# #   # Ap   = read.delim(paste0("species_",i,"_gpp.txt"), header=F, col.names = paste0("V", 1:n))
# #   ids_x = 2:101
# #   ids_h = which(diff(as.numeric(hp[1,ids_x])) > 0)[-1]
# #   
# #   image(x=hp[,1], y=as.numeric(hp[1,ids_h]), z=log(1+1e8*pmax(as.matrix(Up[,ids_h])/max(as.matrix(Up[,ids_h]),na.rm=T),0)), log="", xlab="Time", ylab = "Height distribution, Z*", col=scales::viridis_pal()(100))
# #   matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
# #           las=1, xlab="Time (years)", ylab="Z*", add=T)
# #   
#   image(x=Dp[,1], y=as.numeric(Dp[1,ids_x]), z=log(1+1e8*pmax(as.matrix(Up[,ids_x])/max(as.matrix(Up[,ids_x]),na.rm=T),0)), log="y", xlab="Time", ylab = "Diameter distribution", col=scales::viridis_pal()(100))
# #   
# #   image(x=hp[,1], y=as.numeric(hp[1,ids_h]), z=as.matrix(Lp[,ids_h]), zlim=c(0,8), log="", xlab="Time", ylab = "LAI distribution, Z*", col=scales::viridis_pal()(100))
# #   matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
# #           las=1, xlab="Time (years)", ylab="Z*", add=T)
# #   
# #   image(x=hp[,1], y=as.numeric(hp[1,ids_h]), z=log(1+1e4*pmax(as.matrix(Gp[,ids_h])/max(as.matrix(Gp[,ids_h]),na.rm=T),0)), log="", xlab="Time", ylab = "R. Growth rate dist, Z*", col=scales::viridis_pal()(100))
# #   matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
# #           las=1, xlab="Time (years)", ylab="Z*", add=T)
# #   
# #   image(x=hp[,1], y=as.numeric(hp[1,ids_h]), z=log(1+1e4*pmax(as.matrix(Ap[,ids_h])/max(as.matrix(Ap[,ids_h]),na.rm=T),0)), log="", xlab="Time", ylab = "GPP dist, Z*", col=scales::viridis_pal()(100))
# #   matplot(Zp$V1, Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
# #           las=1, xlab="Time (years)", ylab="Z*", add=T)
# #   
# #   matplot(x=(t(Dp[,-1])), y=(t(Up[,-1]))*1e4/100, lty=1, col=scales::alpha(rainbow(n = nrow((Dp)), start = 0, end = 0.75),0.05), type="l",
# #           las=1, xlab = "Diameter", ylab="Density (Ind/cm/Ha)\n", log="xy", ylim=c(exp(-15), exp(10)))
# #   Up_mean = colMeans(Up[t>100, ids_x]*1e4/100)
# #   points(as.numeric(Up_mean)~as.numeric(Dp[1,ids_x]))
# #   
# }
