#########################################################################
## Script for Sheldon Revisited group biomass estimates
## Authors: Ian Hatton, Ryan Heneghan, Yinon Bar-On, Eric Galbraith
##
##################### MACROZOOPLANKTON #########################
##
## THIS SCRIPT TAKES RAW DATA FOR MACROZOOPLANKTON AND CALCULATES A
## GLOBAL BIOMASS ESTIMATE, AND AN ESTIMATE OF UNCERTAINTY.
## IT ALSO RETURNS THE INTERPOLATED BIOMASS OF MACROZOOPLANKTON IN 1X1 GRID SQUARES
## FOR THE ENTIRE GLOBAL OCEAN
#########################################################################

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

## Macrozoo biom statistical model
zoo_dat <- read.csv("./data/zooplankton/macrozoo_biom_data.csv")
gm1 <- lm(logBiom ~ poly(logchlo, 3) + poly(logdepth, 3) + SST, data = zoo_dat)
summary(gm1)

# Plot model partial residuals and main effects of macrozooplankton biomass model
ylabb = c(expression(paste("log"[10], "(Biomass, mgC m"^-3, ")"), sep  = ""))

plot_list = list()
plot_list[[1]] = resid_plot(gm1, xvar = "SST", ylab = ylabb, xlab = expression(paste("Temperature (", degree, "C)", sep = "")), subtitle = "a)", new_xticks = NA)
plot_list[[2]] = resid_plot(gm1, xvar = "logchlo", ylab = ylabb, xlab = expression(paste("log"[10], "(Chlorophyll, mg m"^-3, ")")), subtitle = "b)", new_xticks = NA)
plot_list[[3]] = resid_plot(gm1, xvar = "logdepth", ylab = ylabb, xlab = expression(paste("log"[10], "(Depth, m)", sep = "")), subtitle = "c)", new_xticks = NA)

ggsave(filename = "./figures/suppfig6_residuals_macrozoo_biomass.png", plot = ggarrange(plots = plot_list, nrow = 1), width = 20, height = 6)


##########################################################
## STEP 2: PREDICT GLOBAL MACROZOOPLANKTON BIOMASS AND 
##         CALCULATE UNCERTAINTY
##########################################################

## Predict global biomass
pred_dat <- read.csv('./data/prediction_data.csv')
pred <- pred_func(gm1, xvars = c("SST", "logchlo", "logdepth"), pred_dat = pred_dat, grid_areas = areas, other_xvar_name = NA, other_xvar_val = NA, over_depth = TRUE, int = 'confidence')

## Save total and top 200m biomass estimate, across body size range, with 95% uncertainty bounds too
tot_macrozoo_biom <- sum(pred$total, na.rm = TRUE)/1e17 # convert from mg C to Pg wet weight
tot_macrozoo_biom_upper200 <- sum(pred$total_top200, na.rm = TRUE)/1e17 # convert from mg C to Pg wet weight

# Save interpolated macrozoo biomass to csv
macrozoo_predictions <- data.frame('lon' = pred_dat$Long, 
                                   'lat' = pred_dat$Lat, 
                                   'macrozoo_biom_gm2_alldepths' = pred$density/100,
                                   'macrozoo_biom_gm2_top200m' = pred$density_top200/100)

write.csv(macrozoo_predictions, './global_map_data/macrozoo_predictions.csv', row.names = FALSE)


## Calculate uncertainty of estimate using bootstrap sampling
pred_dat <- read.csv('./data/prediction_data.csv')
n <- 1000 # 1000 times takes about 10 minutes, you can just run the next line to get the result, then skip the loop
biom_uncert <- 11.3897

boot_save <- rep(0, n)
pb <- txtProgressBar(min = 0, max = n, style = 3)
for(i in 1:n){
  setTxtProgressBar(pb, i)
  boot_choose <- sample(1:dim(zoo_dat)[1], dim(zoo_dat)[1], replace = TRUE)
  boot_dat <- zoo_dat[boot_choose,]
  gm1 <- lm(logBiom ~ poly(logchlo, 3) + poly(logdepth, 3) + SST, data = boot_dat) 
  pred <- pred_func(gm1, xvars = c("SST", "logchlo", "logdepth"), pred_dat = pred_dat, grid_areas = areas, other_xvar_name = NA, other_xvar_val = NA, over_depth = TRUE, int = 'confidence')
  boot_save[i] <- sum(pred$total, na.rm = TRUE)
}

biom_uncert <- 1.96*(10^(0.5*(sd(zoo_dat$logBiom) + sd(log10(boot_save)))))


##########################################################
## STEP 3: SAVE GLOBAL MACROZOOPLANKTON BIOMASS AND 
##         UNCERTAINTY ESTIMATE
##########################################################

## Save uncertainty estimate
group_standard_errors <- read.csv("./summary_output/group_standard_errors.csv")
group_standard_errors[group_standard_errors$Group == "Macrozooplankton", "log10_Standard_Error"] <- biom_uncert
write.csv(group_standard_errors, file = "./summary_output/group_standard_errors.csv", row.names = FALSE)


## Save macrozoo biomass estimate, with uncertainty, aggregated and broken into size bins
macrozoo_sizes <- c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5) # Size bins macrozoo biomass falls into (-2.4 to 2.7 log10 g)
propmacrozoo_biom <- c(0.075,0.2,0.2,0.2,0.2,0.125)  # Proportion of total macrozoo biomass in each size bin (split evenly across log10 size bins)

tot_macrozoo_biom <- c(tot_macrozoo_biom, tot_macrozoo_biom/biom_uncert, tot_macrozoo_biom*biom_uncert)
tot_macrozoo_biom_upper200 <- c(tot_macrozoo_biom_upper200, tot_macrozoo_biom_upper200/biom_uncert, tot_macrozoo_biom_upper200*biom_uncert)

summary_biomass_table <- read.csv("./summary_output/summary_biomass_table.csv")
summary_biomass_table[8, c(2:4)] <- tot_macrozoo_biom
write.csv(summary_biomass_table, file = "./summary_output/summary_biomass_table.csv", row.names = FALSE)

summary_biomass_table_top200 <- read.csv("./summary_output/summary_biomass_table_top200.csv")
summary_biomass_table_top200[8, c(2:4)] <- tot_macrozoo_biom_upper200
write.csv(summary_biomass_table_top200, file = "./summary_output/summary_biomass_table_top200.csv", row.names = FALSE)


tot_macrozoo_biom_bybin <- t(tot_macrozoo_biom %*% t(propmacrozoo_biom))
tot_macrozoo_biom_upper200_bybin <- t(tot_macrozoo_biom_upper200 %*% t(propmacrozoo_biom))

summary_biomass_table_long <- read.csv("./summary_output/summary_biomass_table_long.csv")
summary_biomass_table_long[c(which(summary_biomass_table_long$Group == "Macrozooplankton" & 
                                     summary_biomass_table_long$log10_Size_Midpoint_g %in% macrozoo_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_macrozoo_biom_bybin
write.csv(summary_biomass_table_long, file = "./summary_output/summary_biomass_table_long.csv", row.names = FALSE)

summary_biomass_top200_table_long <- read.csv("./summary_output/summary_biomass_top200_table_long.csv")
summary_biomass_top200_table_long[c(which(summary_biomass_top200_table_long$Group == "Macrozooplankton" & 
                                            summary_biomass_top200_table_long$log10_Size_Midpoint_g %in% macrozoo_sizes)),
                                  c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_macrozoo_biom_upper200_bybin 
write.csv(summary_biomass_top200_table_long, file = "./summary_output/summary_biomass_top200_table_long.csv", row.names = FALSE)





