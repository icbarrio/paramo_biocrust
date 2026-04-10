#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 PÁRAMO BIOCRUST
#       Alejandro Salazar (alejandro@lbhi.is)
#          Isabel C Barrio (isabel@lbhi.is)
#                   Apr 1, 2026
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# this script loads the datasets and runs the analyses presented 
# in the paper by Salazar et al. 

set.seed(132)

# libraries ---- 
# packages to be used
library(tidyverse)
library(lme4)       # to run the models
library(lmerTest)   # to estimate p-values
library(ggplot2)    # to plot
library(gridExtra)
library(extrafont)
library(patchwork)
library(showtext)   # times new roman font


# load datasets ----
# datasets are in separate csv files for cover, surface roughness and stability
# the datasets for include information on:
#   "Elevation": elevation of the sampling site (m a.s.l.)   
#   "Year_edge_glacier": year when the glacier edge was present at the sampling site 
#   "side": position of the sampling grid within the transect (3 possible values:
#       c: center, e: east, w: west; corresponds to grids 1, 2 and 3)
#   "Cover_type": in the cover dataset, type of cover (5 possible values: 
#       Bare_ground, Biocrust, Lichen, Moss, Vasc_plant)
#   "Grid": sampling grid (3 per transect)
#   "Point": for roughness measurements, 9 sampling points were measured per grid
#   "Cover": percent cover of a particular cover type in a grid
#   "Stability": stability index (1-6)
#   "Roughness"

Paramo_biocrust_cover_stats <- read.csv(
    file = "data/Paramo_biocrust_cover_data_GitHub.csv", sep = ",", header = T) %>% 
    rename("Elevation" = "Elevation_long") %>%  # for consistency with the other datasets
    # create a unique identifier per grid
    mutate(GridID = paste0(Year_edge_glacier, "_", Grid))

Paramo_biocrust_roughness_stats<- read.csv(
  file = "data/Paramo_biocrust_roughness_data_GitHub.csv", sep=",", header=T) %>% 
    # create a unique identifier per grid
    mutate(GridID = paste0(Year_edge_glacier, "_", Grid)) 

Paramo_biocrust_stability_stats<- read.csv(
    file = "data/Paramo_biocrust_stability_data_GitHub.csv", sep=",", header=T) %>% 
    # create a unique identifier per grid
    mutate(GridID = paste0(Year_edge_glacier, "_", Grid)) 



# A. Cover ----
## general model ----
# for all models (also when separating by cover type below) we will build 
# models including a non-linear or linear term for elevation 
# and compare which one performs better

# cover model including non-linear term
Cover_VP_B_BG_nonlin <- lmer(Cover ~ poly(Elevation, 2) * Cover_type + (1|GridID),
                   data = Paramo_biocrust_cover_stats)
  summary(Cover_VP_B_BG_nonlin)
    AIC(Cover_VP_B_BG_nonlin)

# cover model without non-linear term
Cover_VP_B_BG <- lmer(Cover ~ Elevation * Cover_type + (1|GridID), 
                      data = Paramo_biocrust_cover_stats)
  summary(Cover_VP_B_BG)
    AIC(Cover_VP_B_BG)

    AIC(Cover_VP_B_BG_nonlin, Cover_VP_B_BG) # compare the two models, first one is better :)

    
## models by cover type ----
# here we do not need to include Grid as a random effect because there is only
# one observation of each cover type per grid (so we use LM)
    
### a) vascular plants ----
# Vasc plant non-linear
Cover_Vasc_plant_cover_nonlin <- lm(Cover ~ poly(Elevation, 2),
                                data = subset(Paramo_biocrust_cover_stats, 
                                              Cover_type == "Vasc_plant"))
  summary(Cover_Vasc_plant_cover_nonlin)
    AIC(Cover_Vasc_plant_cover_nonlin) # 112

# Vasc plant cover linear
Cover_Vasc_plant_cover_lin <- lm(Cover ~ Elevation,
                             data = subset(Paramo_biocrust_cover_stats, 
                                              Cover_type == "Vasc_plant"))
  summary(Cover_Vasc_plant_cover_lin)
    AIC(Cover_Vasc_plant_cover_lin) # 130

       AIC(Cover_Vasc_plant_cover_nonlin, Cover_Vasc_plant_cover_lin) 
       # compare the two models, second one is a little better :)


       
### b) biocrust ----
# Biocrust cover non-linear
Cover_Bioc_cover_nonlin <- lm(Cover ~ poly(Elevation, 2),
                                data = subset(Paramo_biocrust_cover_stats, 
                                              Cover_type == "Biocrust"))
  summary(Cover_Bioc_cover_nonlin)
    AIC(Cover_Bioc_cover_nonlin) # 123

# Biocrust cover linear
Cover_Bioc_cover_lin <- lm(Cover ~ Elevation,
                             data = subset(Paramo_biocrust_cover_stats, 
                                           Cover_type == "Biocrust"))
  summary(Cover_Bioc_cover_lin)
    AIC(Cover_Bioc_cover_lin) # 144

       AIC(Cover_Bioc_cover_nonlin, Cover_Bioc_cover_lin) 
       # compare the two models, first one is slightly better :)

       
### c) moss ----    
# Moss non-linear
Cover_Moss_cover_nonlin <- lm(Cover ~ poly(Elevation, 2),
                                      data = subset(Paramo_biocrust_cover_stats, 
                                                    Cover_type == "Moss"))
  summary(Cover_Moss_cover_nonlin)
    AIC(Cover_Moss_cover_nonlin) # 124

# Moss cover linear
Cover_Moss_cover_lin <- lm(Cover ~ Elevation,
                                   data = subset(Paramo_biocrust_cover_stats, 
                                                 Cover_type == "Moss"))
  summary(Cover_Moss_cover_lin)
    AIC(Cover_Moss_cover_lin) # 144

       AIC(Cover_Moss_cover_nonlin, Cover_Moss_cover_lin) 
       # the two models perform similarly. Elevation is not significant

       
### d) lichen ----
# Lichen non-linear
Cover_Lichen_cover_nonlin <- lm(Cover ~ poly(Elevation, 2),
                                data = subset(Paramo_biocrust_cover_stats, 
                                              Cover_type == "Lichen"))
  summary(Cover_Lichen_cover_nonlin)
    AIC(Cover_Lichen_cover_nonlin) # 119

# Lichen cover linear
Cover_Lichen_cover_lin <- lm(Cover ~ Elevation,
                             data = subset(Paramo_biocrust_cover_stats, 
                                           Cover_type == "Lichen"))
  summary(Cover_Lichen_cover_lin)
    AIC(Cover_Lichen_cover_lin) # 139

       AIC(Cover_Lichen_cover_nonlin, Cover_Lichen_cover_lin) 
       # the two models perform similarly. Elevation is not significant

       
### e) bare ground ----
# Bare_ground non-linear
Cover_Bare_ground_cover_nonlin <- lm(Cover ~ poly(Elevation, 2),
                                  data = subset(Paramo_biocrust_cover_stats, 
                                                Cover_type == "Bare_ground"))
  summary(Cover_Bare_ground_cover_nonlin)
    AIC(Cover_Bare_ground_cover_nonlin) # 122

# Bare_ground cover linear
Cover_Bare_ground_cover_lin <- lm(Cover ~ Elevation,
                               data = subset(Paramo_biocrust_cover_stats, 
                                                Cover_type == "Bare_ground"))
  summary(Cover_Bare_ground_cover_lin)
    AIC(Cover_Bare_ground_cover_lin) # 150

       AIC(Cover_Bare_ground_cover_nonlin, Cover_Bare_ground_cover_lin) 
       # compare the two models, first one is better :)



# B. Roughness ----
# here we include the random effect (GridID) as we have nine measurements per Grid
# roughness model including non-linear term
Roughness_nonlin <- lmer(Roughness ~ poly(Elevation, 2) + (1|GridID),
                             data = Paramo_biocrust_roughness_stats)
  summary(Roughness_nonlin)
    AIC(Roughness_nonlin)

# roughness model without non-linear term
Roughness_lin <- lmer(Roughness ~ Elevation + (1|GridID),
                         data = Paramo_biocrust_roughness_stats)
  summary(Roughness_lin)
    AIC(Roughness_lin)

       AIC(Roughness_nonlin, Roughness_lin) 
       # compare the two models, first one is better :)



# C. Stability ----
# here we include the random effect (GridID) as we have three measurements per Grid
# stability model including non-linear term
Stability_nonlin <- lmer(Stability ~ poly(Elevation, 2) + (1|GridID),
                         data = Paramo_biocrust_stability_stats)
  summary(Stability_nonlin)
    AIC(Stability_nonlin)

# stability model without non-linear term
Stability_non <- lmer(Stability ~ Elevation + (1|GridID),
                         data = Paramo_biocrust_stability_stats)
  summary(Stability_non)
    AIC(Stability_non)

       AIC(Stability_nonlin, Stability_non) 
       # compare the two models, first one is better :)



# Figures -----

## Figure 3 ----
# figure in main text
# load dataset for plotting

Paramo_biocrust_plot <- read.csv(
    file = "data/tundra_paramo_biocrust_data_plot.csv", sep = ",", header = T)

# NOTICE: lowest elevation is 4359 m.a.s.l. Plotted as 4400 m.a.s.l., 
# so the two highest y axis values appear in the plot 
Data_elev_4359 <- subset(Paramo_biocrust_plot, Elevation_masl == 4400)
Data_elev_4648 <- subset(Paramo_biocrust_plot, Elevation_masl == 4648)
Data_elev_4660 <- subset(Paramo_biocrust_plot, Elevation_masl == 4660)
Data_elev_4690 <- subset(Paramo_biocrust_plot, Elevation_masl == 4690)
Data_elev_4698 <- subset(Paramo_biocrust_plot, Elevation_masl == 4698)

tiff('figures/Fig_Paramo_biocrust.tiff', units = "in", width = 8, 
     height = 7, res = 600)

par(mfrow=c(1,3))


### a) cover ----
par(mgp=c(2.5,1,0), 
    mar=c(5,4,2,0))
par(family = "serif")

#Vascular plants
plot(Paramo_biocrust_plot$Cover_vasc_plant_mean, 
     Paramo_biocrust_plot$Elevation_masl, 
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
# Define the labels you want at those positions (years)
tick_labels <- c(1850, 1960, 1990, 2010, 2020)

# Add right axis with manual ticks
axis(2, at = tick_positions, labels = tick_labels, 
     cex.axis=1.5, las=2)

mtext("Year glacier edge", side=2, line=1.5, 
      cex=1.1)

# add segments for each elevation
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

# bare ground
points(Paramo_biocrust_plot$Cover_rock_bare_ground_mean, 
       Paramo_biocrust_plot$Elevation_masl, type="b", col="grey",
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


# moss cover
points(Paramo_biocrust_plot$Cover_moss_mean , 
       Paramo_biocrust_plot$Elevation_masl, type="b", col="palegreen2",
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

# lichens
points(Paramo_biocrust_plot$Cover_lichen_mean, 
       Paramo_biocrust_plot$Elevation_masl, type="b", col="gold",
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

# biocrust
points(Paramo_biocrust_plot$Cover_biocrust_mean, 
       Paramo_biocrust_plot$Elevation_masl, type="b", col=("darkgoldenrod4"),
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



### b) roughness ----
par(mgp=c(2.5,1,0), 
    mar=c(5,2,2,2))
par(family = "serif")

plot(Paramo_biocrust_plot$R_sd_mean, 
     Paramo_biocrust_plot$Elevation_masl, type="b",
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



### c) soil stability ----
par(mgp=c(2.5,1,0), 
    mar=c(5,0,2,4))
par(family = "serif")

plot(Paramo_biocrust_plot$ss_mean,
     Paramo_biocrust_plot$Elevation_masl, type="b", yaxt = "n", ylab="",
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


## Figure Sx ----
# figure in supplementary materials
# load dataset for plotting

Paramo_biocrust_types<- read.csv(
      file = "data/Paramo_biocrust_types_GitHub.csv", sep = ";", header = T) %>% 
      # make sure that the factors are ordered correctly
      mutate(Biocrust_type = factor(Biocrust_type,
              levels = c("Bryophyte", "Lichen", "Cyanobacteria")))

# define custom colors for figure (different types of biocrust)
custom_colors <- c(
  "Cyanobacteria" = "#D2B48C",
  "Lichen" = "#A67C52",
  "Bryophyte" = "#5C4033")

# Max value for top plot
y_max_top <- max(as.numeric(Paramo_biocrust_types$No_hits), na.rm = TRUE)

# Common theme
custom_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    
    text = element_text(family = "serif", size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    
    plot.margin = margin(5, 10, 5, 10)
  )

# plot first the top plot (nr of hits)
No_hits_plot <- ggplot(Paramo_biocrust_types,
  aes(fill = Biocrust_type,
      y = as.numeric(No_hits),
      x = Year_glacier_edge)) +
  geom_bar(position = "dodge", stat = "identity") +
  # add rectangle to show where the biocrust belt is
  annotate("rect",
           xmin = 1952, xmax = 2017,
           ymin = 0, ymax = y_max_top,
           fill = NA, color = "black", linetype = "dashed") +
  # add text to indicate the biocrust belt
  annotate("text",
           x = 1954,
           y = 6.8,
           label = "Biocrust belt",
           hjust = 0, vjust = 1,
           size = 5,
           family = "serif") +
  scale_fill_manual(values = custom_colors, name = "Biocrust type") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "No. of hits") +
  custom_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1))

# and then the bottom plot (percentage of each type of biocrust)
Relative_No_hits_plot <- ggplot(Paramo_biocrust_types,
  aes(fill = Biocrust_type,
      y = as.numeric(No_hits),
      x = Year_glacier_edge)) +
  geom_bar(position = "fill", stat = "identity") +
  # add rectangle to show where the biocrust belt is
  annotate("rect",
           xmin = 1952, xmax = 2017,
           ymin = 0, ymax = 1,
           fill = NA, color = "black", linetype = "dashed") +
  scale_fill_manual(values = custom_colors, name = "Biocrust type") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Year glacier edge",
       y = "Relative no. of hits") +
  custom_theme +
  theme(legend.position = "none")

tiff('figures/Fig_Paramo_biocrust_types.tiff', units = "in", width = 8, 
     height = 7, res = 600)

# combine top and bottom panels
grid.arrange(No_hits_plot,
  Relative_No_hits_plot,
  nrow = 2,
  heights = c(1, 0.9))

No_hits_plot / Relative_No_hits_plot +
  plot_layout(heights = c(1, 0.9))

dev.off()


