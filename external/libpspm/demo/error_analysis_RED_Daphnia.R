library(tidyverse)

setwd("~/codes/libpspm/demo/")

# Error analysis
methods = c("FMU", 
            "IFMU",
            "ILUD", 
            "EBT", 
            "IEBT", 
            "CM", 
            "ICM", 
            "ABM")
cols = scales::alpha(c("purple", "green3", "mediumspringgreen", "darkgoldenrod2", "red3", "pink", "#2b8cbe"), alpha=0.7)
cols = scales::alpha(c("darkgreen", 
                       "yellowgreen",
                       "green3",
                       "magenta",
                       "purple", 
                       "darkgoldenrod2",
                       "darkgoldenrod4",
                       "turquoise2"), alpha=0.7)

cairo_pdf("error_analysis_RED_Daphnia.pdf", height = 7.80, width = 8.00)

par(mfcol = c(2,2), mar=c(4,4,1,1), oma=c(1,1,3,1), cex.lab=1.2, cex.axis=1.2)

err_ebt = read.delim("Daphnia_model/ebt_error_analysis.txt")
err_iebt = read.delim("Daphnia_model/iebt_error_analysis.txt")
err_fmu = read.delim("Daphnia_model/fmu_error_analysis.txt")
err_ifmu = read.delim("Daphnia_model/ifmu_error_analysis.txt")
# err_ifmu2 = read.delim("Daphnia_model/ifmu2_error_analysis.txt")
err_cm = read.delim("Daphnia_model/cm_error_analysis.txt")
err_icm = read.delim("Daphnia_model/icm_error_analysis.txt")
err_abm = read.delim("Daphnia_model/abm_error_analysis.txt")

plot(x=1, y=NA, xlim=c(5,20000), ylim=c(1e-6, 1e2), log="xy", ylab = "Biomass relative error", xlab = "Resolution")
err_fmu %>% with(points(Eb~Nf, type="o", col=cols[1], lwd=2))
err_ifmu %>% with(points(Eb~Nf, type="o", col=cols[2], lwd=2))
# err_ifmu2 %>% with(points(Eb~Nf, type="o", col=cols[3], lwd=2))
err_ebt %>% with(points(Eb~Nf, type="o", col=cols[4], lwd=2))
err_iebt %>% with(points(Eb~Nf, type="o", col=cols[5], lwd=2))
err_cm %>% with(points(Eb~Nf, type="o", col=cols[6], lwd=2))
err_icm %>% with(points(Eb~Nf, type="o", col=cols[7], lwd=2))
err_abm %>% with(points(Eb~Nf, type="o", col=cols[8], lwd=2))
# legend(x = 6, y = 1e-4, legend = methods, col=cols, lty=1, lwd=2)
mtext("Daphnia model", line=1)

plot(x=1, y=NA, xlim=c(5,200000), ylim=c(1e-6, 1e2), log="xy", ylab = "Biomass relative error", xlab = "Execution time (ms)")
err_fmu %>% with(points(Eb~tsys, type="o", col=cols[1], lwd=2))
err_ifmu %>% with(points(Eb~tsys, type="o", col=cols[2], lwd=2))
# err_ifmu2 %>% with(points(Eb~tsys, type="o", col=cols[3], lwd=2))
err_ebt %>% with(points(Eb~tsys, type="o", col=cols[4], lwd=2))
err_iebt %>% with(points(Eb~tsys, type="o", col=cols[5], lwd=2))
err_cm %>% with(points(Eb~tsys, type="o", col=cols[6], lwd=2))
err_icm %>% with(points(Eb~tsys, type="o", col=cols[7], lwd=2))
err_abm %>% with(points(Eb~tsys, type="o", col=cols[8], lwd=2))


## RED

err_ebt = read.delim("RED_model/ebt_error_analysis.txt")
err_iebt = read.delim("RED_model/iebt_error_analysis.txt")
err_fmu = read.delim("RED_model/fmu_error_analysis.txt")
err_ifmu = read.delim("RED_model/ifmu_error_analysis.txt")
# err_ifmu2 = read.delim("RED_model/ifmu2_error_analysis.txt")
err_cm = read.delim("RED_model/cm_error_analysis.txt")
err_icm = read.delim("RED_model/icm_error_analysis.txt")
err_abm = read.delim("RED_model/abm_error_analysis.txt")


plot(x=1, y=NA, xlim=c(1,50000), ylim=c(.5, 150000)/13000, log="xy", ylab = "Biomass relative error", xlab = "Resolution")
err_fmu %>% with(points(Eb~Nf, type="o", col=cols[1], lwd=2))
err_ifmu %>% with(points(Eb~Nf, type="o", col=cols[2], lwd=2))
# err_ifmu2 %>% with(points(Eb~Nf, type="o", col=cols[3], lwd=2))
err_ebt %>% with(points(Eb~Nf, type="o", col=cols[4], lwd=2))
err_iebt %>% with(points(Eb~Nf, type="o", col=cols[5], lwd=2))
err_cm %>% with(points(Eb~Nf, type="o", col=cols[6], lwd=2))
err_icm %>% with(points(Eb~Nf, type="o", col=cols[7], lwd=2))
err_abm %>% with(points(Eb~Nf, type="o", col=cols[8], lwd=2))
# legend(x = 1, y = 1e-2, legend = methods, col=cols, lty=1, lwd=2)
mtext("RED model", line=1)
legend(x = 1, y = 3e-2, legend = methods[-c(3)], col=cols[-c(3)], lty=1, lwd=2, bty = "n")

plot(x=1, y=NA, xlim=c(10,100000), ylim=c(.5, 150000)/13000, log="xy", ylab = "Biomass relative error", xlab = "Execution time (ms)")
err_fmu %>% with(points(Eb~tsys, type="o", col=cols[1], lwd=2))
err_ifmu %>% with(points(Eb~tsys, type="o", col=cols[2], lwd=2))
# err_ifmu2 %>% with(points(Eb~tsys, type="o", col=cols[3], lwd=2))
err_ebt %>% with(points(Eb~tsys, type="o", col=cols[4], lwd=2))
err_iebt %>% with(points(Eb~tsys, type="o", col=cols[5], lwd=2))
err_cm %>% with(points(Eb~tsys, type="o", col=cols[6], lwd=2))
err_icm %>% with(points(Eb~tsys, type="o", col=cols[7], lwd=2))
err_abm %>% with(points(Eb~tsys, type="o", col=cols[8], lwd=2))

dev.off()

