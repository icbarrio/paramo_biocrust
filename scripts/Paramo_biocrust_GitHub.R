setwd("~/R/Collaborations/Tundra_Paramo")
getwd()
dir()



# =================================== Stats ====================================

Paramo_biocrust_cover_stats<- read.csv(file="Paramo_biocrust_cover_data_GitHub.csv", sep=",", header=T) 
attach(Paramo_biocrust_cover_stats)
head(Paramo_biocrust_cover_stats)

Paramo_biocrust_stability_stats<- read.csv(file="Paramo_biocrust_stability_data_GitHub.csv", sep=",", header=T) 
attach(Paramo_biocrust_stability_stats)
head(Paramo_biocrust_stability_stats)

Paramo_biocrust_roughness_stats<- read.csv(file="Paramo_biocrust_roughness_data_GitHub.csv", sep=",", header=T) 
attach(Paramo_biocrust_roughness_stats)
head(Paramo_biocrust_roughness_stats)


library(lme4)
library(lmerTest)


# Cover stats
# ============WITH non-linear term
Cover_VP_B_BG_nonlin <- lmer(Cover ~ poly(Elevation_long, 2) * Cover_type + (1|Grid),
                   data = Paramo_biocrust_cover_stats)
summary(Cover_VP_B_BG_nonlin)
AIC(Cover_VP_B_BG_nonlin)

# ============WITHOUT non-linear term
Cover_VP_B_BG <- lmer(Cover ~ Elevation_long * Cover_type + (1|Grid), 
                      data = Paramo_biocrust_cover_stats)
summary(Cover_VP_B_BG)
AIC(Cover_VP_B_BG)


# ==== Cover models by cover type ==========

# Vasc plant non-linear
Data_Vasc_plant_cover <- subset(Paramo_biocrust_cover_stats, Cover_type == "Vasc_plant")
Cover_Vasc_plant_cover_nonlin <- lmer(Cover ~ poly(Elevation_long, 2) + (1|Grid),
                                data = Data_Vasc_plant_cover)
summary(Cover_Vasc_plant_cover_nonlin)
AIC(Cover_Vasc_plant_cover_nonlin) # 112

# Vasc plant cover linear
Cover_Vasc_plant_cover_lin <- lmer(Cover ~ Elevation_long + (1|Grid),
                             data = Data_Vasc_plant_cover)
summary(Cover_Vasc_plant_cover_lin)
AIC(Cover_Vasc_plant_cover_lin) # 130


# Biocrust cover non-linear
Data_bioc_cover <- subset(Paramo_biocrust_cover_stats, Cover_type == "Biocrust")
Cover_Bioc_cover_nonlin <- lmer(Cover ~ poly(Elevation_long, 2) + (1|Grid),
                                data = Data_bioc_cover)
summary(Cover_Bioc_cover_nonlin)
AIC(Cover_Bioc_cover_nonlin) # 123

# Biocrust cover linear
Cover_Bioc_cover_lin <- lmer(Cover ~ Elevation_long + (1|Grid),
                             data = Data_bioc_cover)
summary(Cover_Bioc_cover_lin)
AIC(Cover_Bioc_cover_lin) # 144

# Moss non-linear
Data_Moss_cover <- subset(Paramo_biocrust_cover_stats, Cover_type == "Moss")
Cover_Moss_cover_nonlin <- lmer(Cover ~ poly(Elevation_long, 2) + (1|Grid),
                                      data = Data_Moss_cover)
summary(Cover_Moss_cover_nonlin)
AIC(Cover_Moss_cover_nonlin) # 124

# Moss cover linear
Cover_Moss_cover_lin <- lmer(Cover ~ Elevation_long + (1|Grid),
                                   data = Data_Moss_cover)
summary(Cover_Moss_cover_lin)
AIC(Cover_Moss_cover_lin) # 144


# Lichen non-linear
Data_Lichen_cover <- subset(Paramo_biocrust_cover_stats, Cover_type == "Lichen")
Cover_Lichen_cover_nonlin <- lmer(Cover ~ poly(Elevation_long, 2) + (1|Grid),
                                data = Data_Lichen_cover)
summary(Cover_Lichen_cover_nonlin)
AIC(Cover_Lichen_cover_nonlin) # 119

# Lichen cover linear
Cover_Lichen_cover_lin <- lmer(Cover ~ Elevation_long + (1|Grid),
                             data = Data_Lichen_cover)
summary(Cover_Lichen_cover_lin)
AIC(Cover_Lichen_cover_lin) # 139


# Bare_ground non-linear
Data_Bare_ground_cover <- subset(Paramo_biocrust_cover_stats, Cover_type == "Bare_ground")
Cover_Bare_ground_cover_nonlin <- lmer(Cover ~ poly(Elevation_long, 2) + (1|Grid),
                                  data = Data_Bare_ground_cover)
summary(Cover_Bare_ground_cover_nonlin)
AIC(Cover_Bare_ground_cover_nonlin) # 122

# Bare_ground cover linear
Cover_Bare_ground_cover_lin <- lmer(Cover ~ Elevation_long + (1|Grid),
                               data = Data_Bare_ground_cover)
summary(Cover_Bare_ground_cover_lin)
AIC(Cover_Bare_ground_cover_lin) # 150



# Roughness stats  
# ============WITH non-linear term
Roughness_nonlin <- lmer(Roughness ~ poly(Elevation, 2) + (1|Grid),
                             data = Paramo_biocrust_roughness_stats)
summary(Roughness_nonlin)
AIC(Roughness_nonlin)

# ============WITHOUT non-linear term
Roughness_lin <- lmer(Roughness ~ Elevation + (1|Grid),
                         data = Paramo_biocrust_roughness_stats)
summary(Roughness_lin)
AIC(Roughness_lin)


# Stability stats
# ============WITH non-linear term
Stability_nonlin <- lmer(Stability ~ poly(Elevation, 2) + (1|Grid),
                         data = Paramo_biocrust_stability_stats)
summary(Stability_nonlin)
AIC(Stability_nonlin)

# ============WITHOUT non-linear term
Stability_non <- lmer(Stability ~ Elevation + (1|Grid),
                         data = Paramo_biocrust_stability_stats)
summary(Stability_non)
AIC(Stability_non)




# =================================== Plot =====================================

Paramo_biocrust_plot<- read.csv(file="tundra_paramo_biocrust_data_plot.csv", sep=",", header=T) 
attach(Paramo_biocrust_plot)
head(Paramo_biocrust_plot)

# NOTICE: lowest elevation is 4359 m.a.s.l. Plotted as 4400 m.a.s.l., so the two highest y axis values appear in the plot 
Data_elev_4359 <- subset(Paramo_biocrust_plot, Elevation_masl == 4400)
Data_elev_4648 <- subset(Paramo_biocrust_plot, Elevation_masl == 4648)
Data_elev_4660 <- subset(Paramo_biocrust_plot, Elevation_masl == 4660)
Data_elev_4690 <- subset(Paramo_biocrust_plot, Elevation_masl == 4690)
Data_elev_4698 <- subset(Paramo_biocrust_plot, Elevation_masl == 4698)




tiff('Fig_Paramo_biocrust.tiff', units="in", width=8, 
     height=7, res=600)

par(mfrow=c(1,3))


# ===================== Cover ============================================

par(mgp=c(2.5,1,0), 
    mar=c(5,4,2,0))
par(family = "serif")

#Vascular plants
plot(Cover_vasc_plant_mean, Elevation_masl, 
     type="b", col="darkgreen", xlim=c(0,110),
     xlab="Cover (%)", ylab="",
     cex.axis=1.6, cex.lab=1.6,
     pch=17, cex=1.6,
     yaxt = "n")

# Define manual tick positions (elevations)
tick_positions <- c(4400, 4648, 4660, 4690, 4698)

axis(side = 4, 
     at = tick_positions,
     labels = FALSE) 

# Right axis with elevation data

# Define manual tick positions (elevations)
#tick_positions <- c(4500, 4648, 4660, 4690, 4698)

# Define the labels you want at those positions (years)
tick_labels <- c(1850, 1960, 1990, 2010, 2020)

# Add right axis with manual ticks
axis(2, at = tick_positions, labels = tick_labels, 
     cex.axis=1.5, las=2)

mtext("Year glacier edge", side=2, line=1.5, 
      cex=1.1)



segments((Data_elev_4359$Cover_vasc_plant_mean + Data_elev_4359$Cover_vasc_plant_se), 4400, 
         (Data_elev_4359$Cover_vasc_plant_mean - Data_elev_4359$Cover_vasc_plant_se), 4400, col="darkgreen")
segments((Data_elev_4648$Cover_vasc_plant_mean + Data_elev_4648$Cover_vasc_plant_se), 4648, 
         (Data_elev_4648$Cover_vasc_plant_mean - Data_elev_4648$Cover_vasc_plant_se), 4648, col="darkgreen")
segments((Data_elev_4660$Cover_vasc_plant_mean + Data_elev_4660$Cover_vasc_plant_se), 4660, 
         (Data_elev_4660$Cover_vasc_plant_mean - Data_elev_4660$Cover_vasc_plant_se), 4660, col="darkgreen")
segments((Data_elev_4690$Cover_vasc_plant_mean + Data_elev_4690$Cover_vasc_plant_se), 4690, 
         (Data_elev_4690$Cover_vasc_plant_mean - Data_elev_4690$Cover_vasc_plant_se), 4690, col="darkgreen")
segments((Data_elev_4698$Cover_vasc_plant_mean + Data_elev_4698$Cover_vasc_plant_se), 4698, 
         (Data_elev_4698$Cover_vasc_plant_mean - Data_elev_4698$Cover_vasc_plant_se), 4698, col="darkgreen")


points(Cover_rock_bare_ground_mean , Elevation_masl, type="b", col="grey",
       pch=15, cex=1.6)
segments((Data_elev_4359$Cover_rock_bare_ground_mean + Data_elev_4359$Cover_rock_bare_ground_se), 4400, 
         (Data_elev_4359$Cover_rock_bare_ground_mean - Data_elev_4359$Cover_rock_bare_ground_se), 4400, col="grey")
segments((Data_elev_4648$Cover_rock_bare_ground_mean + Data_elev_4648$Cover_rock_bare_ground_se), 4648, 
         (Data_elev_4648$Cover_rock_bare_ground_mean - Data_elev_4648$Cover_rock_bare_ground_se), 4648, col="grey")
segments((Data_elev_4660$Cover_rock_bare_ground_mean + Data_elev_4660$Cover_rock_bare_ground_se), 4660, 
         (Data_elev_4660$Cover_rock_bare_ground_mean - Data_elev_4660$Cover_rock_bare_ground_se), 4660, col="grey")
segments((Data_elev_4690$Cover_rock_bare_ground_mean + Data_elev_4690$Cover_rock_bare_ground_se), 4690, 
         (Data_elev_4690$Cover_rock_bare_ground_mean - Data_elev_4690$Cover_rock_bare_ground_se), 4690, col="grey")
segments((Data_elev_4698$Cover_rock_bare_ground_mean + Data_elev_4698$Cover_rock_bare_ground_se), 4698, 
         (Data_elev_4698$Cover_rock_bare_ground_mean - Data_elev_4698$Cover_rock_bare_ground_se), 4698, col="grey")
par(xpd = TRUE)  
text(1, 4717, "A", cex = 1.5)


# Moss cover
points(Cover_moss_mean , Elevation_masl, type="b", col="palegreen2",
       pch=18, cex=1.6)
segments((Data_elev_4359$Cover_moss_mean + Data_elev_4359$Cover_moss_se), 4400, 
         (Data_elev_4359$Cover_moss_mean - Data_elev_4359$Cover_moss_se), 4400, col="palegreen2")
segments((Data_elev_4648$Cover_moss_mean + Data_elev_4648$Cover_moss_se), 4648, 
         (Data_elev_4648$Cover_moss_mean - Data_elev_4648$Cover_moss_se), 4648, col="palegreen2")
segments((Data_elev_4660$Cover_moss_mean + Data_elev_4660$Cover_moss_se), 4660, 
         (Data_elev_4660$Cover_moss_mean - Data_elev_4660$Cover_moss_se), 4660, col="palegreen2")
segments((Data_elev_4690$Cover_moss_mean + Data_elev_4690$Cover_moss_se), 4690, 
         (Data_elev_4690$Cover_moss_mean - Data_elev_4690$Cover_moss_se), 4690, col="palegreen2")
segments((Data_elev_4698$Cover_moss_mean + Data_elev_4698$Cover_moss_se), 4698, 
         (Data_elev_4698$Cover_moss_mean - Data_elev_4698$Cover_moss_se), 4698, col="palegreen2")

# Lichen cover
points(Cover_lichen_mean , Elevation_masl, type="b", col="gold",
       pch=25, bg="gold",cex=1.6)
segments((Data_elev_4359$Cover_lichen_mean + Data_elev_4359$Cover_lichen_se), 4400, 
         (Data_elev_4359$Cover_lichen_mean - Data_elev_4359$Cover_lichen_se), 4400, col="gold")
segments((Data_elev_4648$Cover_lichen_mean + Data_elev_4648$Cover_lichen_se), 4648, 
         (Data_elev_4648$Cover_lichen_mean - Data_elev_4648$Cover_lichen_se), 4648, col="gold")
segments((Data_elev_4660$Cover_lichen_mean + Data_elev_4660$Cover_lichen_se), 4660, 
         (Data_elev_4660$Cover_lichen_mean - Data_elev_4660$Cover_lichen_se), 4660, col="gold")
segments((Data_elev_4690$Cover_lichen_mean + Data_elev_4690$Cover_lichen_se), 4690, 
         (Data_elev_4690$Cover_lichen_mean - Data_elev_4690$Cover_lichen_se), 4690, col="gold")
segments((Data_elev_4698$Cover_lichen_mean + Data_elev_4698$Cover_lichen_se), 4698, 
         (Data_elev_4698$Cover_lichen_mean - Data_elev_4698$Cover_lichen_se), 4698, col="gold")


#Biocrust
points(Cover_biocrust_mean, Elevation_masl, type="b", col=("darkgoldenrod4"),
       pch=16, cex=1.6)
segments((Data_elev_4359$Cover_biocrust_mean + Data_elev_4359$Cover_biocrust_se), 4400, 
         (Data_elev_4359$Cover_biocrust_mean - Data_elev_4359$Cover_biocrust_se), 4400, col="darkgoldenrod4")
segments((Data_elev_4648$Cover_biocrust_mean + Data_elev_4648$Cover_biocrust_se), 4648, 
         (Data_elev_4648$Cover_biocrust_mean - Data_elev_4648$Cover_biocrust_se), 4648, col="darkgoldenrod4")
segments((Data_elev_4660$Cover_biocrust_mean + Data_elev_4660$Cover_biocrust_se), 4660, 
         (Data_elev_4660$Cover_biocrust_mean - Data_elev_4660$Cover_biocrust_se), 4660, col="darkgoldenrod4")
segments((Data_elev_4690$Cover_biocrust_mean + Data_elev_4690$Cover_biocrust_se), 4690, 
         (Data_elev_4690$Cover_biocrust_mean - Data_elev_4690$Cover_biocrust_se), 4690, col="darkgoldenrod4")
segments((Data_elev_4698$Cover_biocrust_mean + Data_elev_4698$Cover_biocrust_se), 4698, 
         (Data_elev_4698$Cover_biocrust_mean - Data_elev_4698$Cover_biocrust_se), 4698, col="darkgoldenrod4")


par(xpd = TRUE)  
text(1, 4717, "A", cex = 1.5)



# ================= Roughness (R) =======================

par(mgp=c(2.5,1,0), 
    mar=c(5,2,2,2))
par(family = "serif")

plot(R_sd_mean, Elevation_masl, type="b",
     xlim=c(-0.5,4.5), yaxt = "n", ylab="",
     xlab="Surface roughness (SD, cm)",
     cex.axis=1.6, cex.lab=1.6,
     pch=16, cex=1.6)
axis(2, at = tick_positions, labels = FALSE,
     cex.axis=1.3)
axis(side = 4, 
     at = tick_positions,
     labels = FALSE) 

segments((Data_elev_4359$R_sd_mean + Data_elev_4359$R_sd_se), 4400, 
         (Data_elev_4359$R_sd_mean - Data_elev_4359$R_sd_se), 4400)
segments((Data_elev_4648$R_sd_mean + Data_elev_4648$R_sd_se), 4648, 
         (Data_elev_4648$R_sd_mean - Data_elev_4648$R_sd_se), 4648)
segments((Data_elev_4660$R_sd_mean + Data_elev_4660$R_sd_se), 4660, 
         (Data_elev_4660$R_sd_mean - Data_elev_4660$R_sd_se), 4660)
segments((Data_elev_4690$R_sd_mean + Data_elev_4690$R_sd_se), 4690, 
         (Data_elev_4690$R_sd_mean - Data_elev_4690$R_sd_se), 4690)
segments((Data_elev_4698$R_sd_mean + Data_elev_4698$R_sd_se), 4698, 
         (Data_elev_4698$R_sd_mean - Data_elev_4698$R_sd_se), 4698)

par(xpd = TRUE)  
text(-0.5, 4717, "B", cex = 1.5)



# =========================  Soil stability ================================
par(mgp=c(2.5,1,0), 
    mar=c(5,0,2,4))
par(family = "serif")

plot(ss_mean,Elevation_masl, type="b", yaxt = "n", ylab="",
     xlim=c(1,6), xlab="Surface stability index",
     cex.axis=1.6, cex.lab=1.6,
     pch=16, cex=1.6)
axis(2, at = tick_positions, labels = FALSE)

# Define manual tick positions (elevations)
tick_positions_r <- c(4400, 4648, 4660, 4690, 4698)

# Define the labels you want at those positions (years)
tick_labels_r <- c(4359, 4648, 4660, 4690, 4698)


axis(side = 4,
     at = tick_positions_r,
     labels = tick_labels_r,
     cex.axis=1.5,
     las=2) 


segments((Data_elev_4698$ss_mean + Data_elev_4698$ss_se), 4698, 
         (Data_elev_4698$ss_mean - Data_elev_4698$ss_se), 4698)
segments(Data_elev_4660$ss_mean, 4656, Data_elev_4359$ss_mean, 4404)



# Right axis with elevation data


mtext("Elevation (m.a.s.l.)", side=4, line=1.5, 
      cex=1.1)

par(xpd = TRUE)  
text(1, 4717, "C", cex = 1.5, cex.lab=1.3)


dev.off()