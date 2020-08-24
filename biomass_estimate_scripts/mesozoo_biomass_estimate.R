#################################################################################
## Script for Sheldon Revisited group biomass estimates
## Authors: Ian Hatton, Ryan Heneghan, Yinon Bar-On, Eric Galbraith
##
##################### MESOZOOPLANKTON ###########################################
##
## THIS SCRIPT TAKES RAW DATA FOR MESOZOOPLANKTON AND GIANT RHIZARIA 
## AND CALCULATES A GLOBAL BIOMASS ESTIMATE, AND AN ESTIMATE OF UNCERTAINTY.
## IT ALSO RETURNS THE INTERPOLATED BIOMASS OF MESOZOOPLANKTON 
## (RHIZARIA + OTHER MESOZOO) IN 1X1 GRID SQUARES FOR THE ENTIRE GLOBAL OCEAN
#################################################################################

rm(list=ls())
# To plot figures of main effects and partial residuals
library(ggplot2) 
library(egg) 
library(visreg)

# For area of 360x180 global grid squares
library(raster)

# Set working directory to scripts_data_figures folder
#setwd("~/Desktop/Papers/Sheldon_Revisited/scripts_data_figures/")

# Source prediction data file for global biomass/abundance interpolation
pred_dat <- read.csv('./data/prediction_data.csv')

# Source area file of 180x360 global grid (in m2) for global biomass/abundance calculation
areas <- as.vector(raster::area(raster()))*1000*1000 # Areas of 1x1 grid squares, in m2

# Source partial residual/main effect plot function
source("./processing_scripts/residual_plot_function.R")

# Source interpolation function for global abundance/biomass
source("./processing_scripts/interpolate_global_abund_biomass.R")

##########################################################
## STEP 1: FIT STATISTICAL MODEL TO BIOMASS OBSERVATIONS
##########################################################
## Non rhizaria mesozoo
zoo_dat <- read.csv("./data/zooplankton/mesozoo_biom_data.csv")

gm1 <- lm(logBiom ~ poly(logchlo, 3) + poly(logdepth, 3) + SST + BiomassMethod, data = zoo_dat)
summary(gm1)

# PLOT MODEL PARTIAL RESIDUALS
ylabb = c(expression(paste("log"[10], "(Biomass, mgC m"^-3, ")"), sep  = ""))

plot_list = list()
plot_list[[1]] = resid_plot(gm1, xvar = "SST", ylab = ylabb, xlab = expression(paste("Temperature (", degree, "C)", sep = "")), subtitle = "a)", new_xticks = NA)
plot_list[[2]] = resid_plot(gm1, xvar = "logchlo", ylab = ylabb, xlab = expression(paste("log"[10], "(Chlorophyll, mg m"^-3, ")")), subtitle = "b)", new_xticks = NA)
plot_list[[3]] = resid_plot(gm1, xvar = "logdepth", ylab = ylabb, xlab = expression(paste("log"[10], "(Depth, m)", sep = "")), subtitle = "c)", new_xticks = NA)
plot_list[[4]] = resid_plot(gm1, xvar = "BiomassMethod", ylab = ylabb, xlab = 'Measurement Method', subtitle = "d)", new_xticks = c("Ash-free Dry", "Carbon", "CHN Carbon", "Displacement", "Dry", "Settled", "Wet"))

ggsave(filename = "./figures/suppfig4_residuals_mesozoo_biomass.png", plot = ggarrange(plots = plot_list, nrow = 2), width = 13.3, height = 12)

## Giant rhizaria
# Giant rhizaria 0-200m samples
zoo_dat2 <- read.csv("./data/zooplankton/rhizaria_biomtop200_data.csv")

gm1_rhiz1 <- lm(logBiom ~ poly(logchlo, 2) + poly(SST, 2), data = zoo_dat2)
summary(gm1_rhiz1)

# Giant rhizaria 0-500m samples
zoo_dat3 <- read.csv("./data/zooplankton/rhizaria_biomtop500_data.csv")

gm1_rhiz2 <- lm(logBiom ~ poly(logchlo, 2) + poly(SST, 2), data = zoo_dat3)
summary(gm1_rhiz2)

# Plot partial residuals and main effect plots for rhizaria 0-200 and 0-500 models
ylabb1 = c(expression(paste("log"[10], "(Biomass 0-200m, mgC m"^-3, ")"), sep  = ""))
ylabb2 = c(expression(paste("log"[10], "(Biomass 0-500m, mgC m"^-3, ")"), sep  = ""))

plot_list = list()
plot_list[[1]] = resid_plot(gm1_rhiz1, xvar = "SST", ylab = ylabb1, xlab = expression(paste("Temperature (", degree, "C)", sep = "")), subtitle = "a)", new_xticks = NA)
plot_list[[2]] = resid_plot(gm1_rhiz1, xvar = "logchlo", ylab = ylabb1, xlab = expression(paste("log"[10], "(Chlorophyll, mg m"^-3, ")")), subtitle = "b)", new_xticks = NA)
plot_list[[3]] = resid_plot(gm1_rhiz2, xvar = "SST", ylab = ylabb2, xlab = expression(paste("Temperature (", degree, "C)", sep = "")), subtitle = "c)", new_xticks = NA)
plot_list[[4]] = resid_plot(gm1_rhiz2, xvar = "logchlo", ylab = ylabb2, xlab = expression(paste("log"[10], "(Chlorophyll, mg m"^-3, ")")), subtitle = "d)", new_xticks = NA)


ggsave(filename = "./figures/suppfig5_residuals_rhizaria_biomass.png", plot = ggarrange(plots = plot_list, nrow = 2), width = 13.3, height = 12)


##########################################################
## STEP 2: PREDICT GLOBAL MESOZOOPLANKTON BIOMASS AND 
##         CALCULATE UNCERTAINTY
##########################################################

pred_dat <- read.csv('./data/prediction_data.csv')

## Predict global non-rhizaria mesozoo biomass
pred_miscmeso <- pred_func(gm1, xvars = c("SST", "logchlo", "logdepth"), pred_dat = pred_dat, grid_areas = areas, other_xvar_name = "BiomassMethod", other_xvar_val = "Wet", over_depth = TRUE, int = 'confidence')

## Predict global rhizaria 0-200m mesozoo biomass
pred_rhiz0to200 <- pred_func(gm1_rhiz1, xvars = c("SST", "logchlo"), pred_dat = pred_dat, grid_areas = areas, other_xvar_name = "NA", other_xvar_val = "NA", over_depth = FALSE, int = 'confidence')

## Predict global rhizaria 0-500m mesozoo biomass
pred_rhiz0to500 <- pred_func(gm1_rhiz2, xvars = c("SST", "logchlo"), pred_dat = pred_dat, grid_areas = areas, other_xvar_name = "NA", other_xvar_val = "NA", over_depth = FALSE, int = 'confidence')

## Save total and top 200m biomass estimate, across body size range, with 95% uncertainty bounds too
tot_mesozoo_biom <- sum(pred_miscmeso$total, na.rm = TRUE)/1e18 # convert from mg to Pg wet weight
tot_mesozoo_biom_upper200 <- sum(pred_miscmeso$total_top200, na.rm = TRUE)/1e18 # convert from mg to Pg wet weight

tot_rhiz_biom <- sum(pred_rhiz0to500$total, na.rm = TRUE)/1e17 # convert from mg C to Pg wet weight
tot_rhiz_biom_upper200 <- sum(pred_rhiz0to200$total, na.rm = TRUE)/1e17 # convert from mg C to Pg wet weight

# Save interpolated mesozoo biomass to csv
mesozoo_predictions <- data.frame('lon' = pred_dat$Long, 
                                   'lat' = pred_dat$Lat, 
                                   'mesozoo_biom_gm2_alldepths' = (pred_miscmeso$density/1000 + pred_rhiz0to500$density/100),
                                   'mesozoo_biom_gm2_top200m' = (pred_miscmeso$density_top200/1000 + pred_rhiz0to200$density/100))

write.csv(mesozoo_predictions, './global_map_data/mesozoo_predictions.csv', row.names = FALSE)

## Calculate uncertainty of estimate using bootstrap sampling
pred_dat <- read.csv('./data/prediction_data.csv')
n <- 1000 # 1000 times takes about 10 minutes, you can just run the next line to get the result, then skip the loop
biom_uncert <- 7.7

boot_save <- matrix(0, nrow = n, ncol = 2)
pb <- txtProgressBar(min = 0, max = n, style = 3)
for(i in 1:n){
  setTxtProgressBar(pb, i)
  boot_choose1 <- sample(1:dim(zoo_dat)[1], dim(zoo_dat)[1], replace = TRUE)
  boot_choose2 <- sample(1:dim(zoo_dat3)[1], dim(zoo_dat3)[1], replace = TRUE)
  boot_dat1 <- zoo_dat[boot_choose1,]
  boot_dat2 <- zoo_dat3[boot_choose2,]
  gm1 <- lm(logBiom ~ poly(logchlo, 3) + poly(logdepth, 3) + SST + BiomassMethod, data = boot_dat1)
  gm2 <- lm(logBiom ~ poly(logchlo, 2) + poly(SST, 2), data = boot_dat2)
  
  pred1 <- pred_func(gm1, xvars = c("SST", "logchlo", "logdepth"), pred_dat = pred_dat, grid_areas = areas, other_xvar_name = "BiomassMethod", other_xvar_val = "Wet", over_depth = TRUE, int = 'confidence')
  pred2 <- pred_func(gm2, xvars = c("SST", "logchlo"), pred_dat = pred_dat, grid_areas = areas, other_xvar_name = "NA", other_xvar_val = "NA", over_depth = FALSE, int = 'confidence')
  boot_save[i,1] <- sum(pred1$total, na.rm = TRUE)/1e18
  boot_save[i,2] <- sum(pred2$total, na.rm = TRUE)/1e17
}

biom_uncert1 <- 1.96*(10^(0.5*(sd(zoo_dat$logBiom) + sd(log10(boot_save[,1])))))
biom_uncert2 <- 1.96*(10^(0.5*(sd(zoo_dat3$logBiom) + sd(log10(boot_save[,2])))))

biom_uncert <- (tot_mesozoo_biom*biom_uncert1 + tot_rhiz_biom*biom_uncert2)/(tot_mesozoo_biom + tot_rhiz_biom)

##########################################################
## STEP 3: SAVE GLOBAL MACROZOOPLANKTON BIOMASS AND 
##         UNCERTAINTY ESTIMATE
##########################################################

## Save uncertainty estimate
group_standard_errors <- read.csv("./summary_output/group_standard_errors.csv")
group_standard_errors[group_standard_errors$Group == "Mesozooplankton", "log10_Standard_Error"] <- biom_uncert
write.csv(group_standard_errors, file = "./summary_output/group_standard_errors.csv", row.names = FALSE)


## Save mesozoo biomass estimate, with uncertainty, aggregated and broken into size bins
mesozoo_sizes <- c(-5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5) # Size bins mesozoo biomass falls into (-5.4 to 0.6 log10 g)
propmesozoo_biom <- c(0.083,0.167,0.167,0.167,0.167,0.167,0.083)  # Proportion of total mesozoo biomass in each size bin (size range of nanozoo is 10^-14 to 10^-11.5g, split evenly across log10 size bins)

tot_mesozoo_biom_range <- c(tot_mesozoo_biom + tot_rhiz_biom, (tot_mesozoo_biom + tot_rhiz_biom)/biom_uncert, (tot_mesozoo_biom + tot_rhiz_biom)*biom_uncert)
tot_mesozoo_biom_upper200_range <- c(tot_mesozoo_biom_upper200 + tot_rhiz_biom_upper200, (tot_mesozoo_biom_upper200 + tot_rhiz_biom_upper200)/biom_uncert, (tot_mesozoo_biom_upper200 + tot_rhiz_biom_upper200)*biom_uncert)

summary_biomass_table <- read.csv("./summary_output/summary_biomass_table.csv")
summary_biomass_table[7, c(2:4)] <- tot_mesozoo_biom_range
write.csv(summary_biomass_table, file = "./summary_output/summary_biomass_table.csv", row.names = FALSE)

summary_biomass_table_top200 <- read.csv("./summary_output/summary_biomass_table_top200.csv")
summary_biomass_table_top200[7, c(2:4)] <- tot_mesozoo_biom_upper200_range
write.csv(summary_biomass_table_top200, file = "./summary_output/summary_biomass_table_top200.csv", row.names = FALSE)


tot_mesozoo_biom_bybin <- t(tot_mesozoo_biom_range %*% t(propmesozoo_biom))
tot_mesozoo_biom_upper200_bybin <- t(tot_mesozoo_biom_upper200_range %*% t(propmesozoo_biom))

summary_biomass_table_long <- read.csv("./summary_output/summary_biomass_table_long.csv")
summary_biomass_table_long[c(which(summary_biomass_table_long$Group == "Mesozooplankton" & 
                                     summary_biomass_table_long$log10_Size_Midpoint_g %in% mesozoo_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_mesozoo_biom_bybin
write.csv(summary_biomass_table_long, file = "./summary_output/summary_biomass_table_long.csv", row.names = FALSE)

summary_biomass_top200_table_long <- read.csv("./summary_output/summary_biomass_top200_table_long.csv")
summary_biomass_top200_table_long[c(which(summary_biomass_top200_table_long$Group == "Mesozooplankton" & 
                                            summary_biomass_top200_table_long$log10_Size_Midpoint_g %in% mesozoo_sizes)),
                                  c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_mesozoo_biom_upper200_bybin
write.csv(summary_biomass_top200_table_long, file = "./summary_output/summary_biomass_top200_table_long.csv", row.names = FALSE)


