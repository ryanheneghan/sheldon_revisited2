#########################################################################
## Script for Sheldon Revisited group biomass estimates
## Authors: Ian Hatton, Ryan Heneghan, Yinon Bar-On, Eric Galbraith
##
############################ PHYTOPLANKTON ##############################
##
## THIS SCRIPT CALCULATES A GLOBAL BIOMASS ESTIMATE OF PHYTOPLANKTON
## (PICO, NANO AND MICRO) AND AN ESTIMATE OF UNCERTAINTY.
## IT ALSO RETURNS THE INTERPOLATED BIOMASS OF PHYTOPLANKTON IN 
## 1X1 GRID SQUARES FOR THE ENTIRE GLOBAL OCEAN
#########################################################################

rm(list=ls())

# Set working directory to scripts_data_figures folder
#setwd("~/Desktop/Papers/Sheldon_Revisited/scripts_data_figures/")

# To open netcdfs of chlorophyll, euphotic zone depth and mixed layer depth
library(ncdf4)

# To obtain grid of surface area of 180x360 global grid
library(raster)

# Source area file of 360x180 global grid (in m2) for global biomass/abundance calculation
areas <- as.matrix(raster::area(raster()))*1000*1000 # Areas of 1x1 grid squares, in m2
areas <- t(areas)

##########################################################
### STEP 1: WORK OUT CONTRIBUTION OF PICO, NANO 
###         AND MICROPHYTOPLANKTON TO TOTAL CHLOROPHYLL, 
###         ACROSS 1X1 GRID OF GLOBAL OCEAN
##########################################################

# Import chlorophyll
chlo_file <- nc_open(list.files(path = "./data/", pattern = '*chlo_modis*', full.names = TRUE))
chlo <- ncvar_get(chlo_file, "chlor_a")

# Equations from Brewin et al., 2010, doi:10.1016/j.ecolmodel.2010.02.014
pico <- (0.13*(1-exp(-0.8/0.13*chlo)))/chlo
nano <- (0.77*(1-exp(-0.94/0.77*chlo)))/chlo - pico
micro <- (chlo - 0.77*(1-exp(-0.94/0.77*chlo)))/chlo


##########################################################
### STEP 2: CALCULATE TOTAL CHLOROPHYLL IN EUPHOTIC ZONE
##########################################################

# We use the equations in table 4 of Uitz et al., 2006 doi:10.1029/2005JC003207
# Import euphotic zone depth
euph_file <- nc_open(list.files(path = "./data/", pattern = '*zeu_lee_*', full.names = TRUE))
euphs <- ncvar_get(euph_file, 'Zeu_lee')

# Import mixed layer depth
mld_file <- nc_open(list.files(path = "./data/", pattern = 'Argo_*', full.names = TRUE))
mlds <- ncvar_get(mld_file, "mld_dt_mean")
mlds <- apply(mlds[,c(181:360,1:180),], c(2,3), mean, na.rm = TRUE) # Calculate annual average mixed layer depth

# Areas where mld is not available - mostly arctic waters, make them stratified (mld > euph zone depth). I did this to be conservative
# as mixed waters have higher chlorophyll than stratified typically
mlds[which(is.na(euphs))] <- 10000
mlds[which(is.na(mlds))] <- 20000
mlds[c(mlds == 20000)] <- 1.2*euphs[c(mlds == 20000)]
mlds[c(mlds == 10000)] <- NA

# Which grid squares are stratified, which are mixed
mix_regime <- euphs/mlds
mix_regime[mix_regime >= 1] <- 'strat'
mix_regime[mix_regime < 1] <- 'mixed'
mix_regime[which(is.na(mix_regime) & !is.na(chlo) == TRUE)] <- 'mixed'

# Calculate total chlorophyll a in euphotic zone
chlo_tot <- matrix(NA, nrow = dim(chlo)[1], ncol = dim(chlo)[2])
chlo_tot[which(mix_regime == 'strat' & chlo <= 1)] <- 36.1*chlo[which(mix_regime == 'strat' & chlo <= 1)]^0.357
chlo_tot[which(mix_regime == 'strat' & chlo > 1)] <- 37.7*chlo[which(mix_regime == 'strat' & chlo > 1)]^0.615
chlo_tot[which(mix_regime == 'mixed')] <- 42.1*chlo[which(mix_regime == 'mixed')]^0.538

# Uncertainty of chlo surface to chlo total conversion
chlo_tot_error <- 1.4 # 1.4 fold 95% prediction interval for total chlorophyll from Uitz 2006. Calculated by digitising figure 3

##########################################################
### STEP 3: CONVERT TOTAL CHLOROPHYLL IN EUPHOTIC ZONE TO
###         CARBON, THEN WET WEIGHT BIOMASS OF PICO,
###         NANO AND MICROPHYTOPLANKTON
##########################################################

chlo_ave <- chlo_tot/euphs # Average chlorophyll a concentration in euphotic zone
# Use chlo2carb relationship from Maranon et al., 2014 - doi:10.1371/journal.pone.0099312
carb_ave <- ((10^1.79)*(chlo_ave)^0.89) # Average phyto biomass (carbon) in euphotic zone (mg C m-3)
carb_tot <- carb_ave*euphs # Total carbon in euphotic zone (mg C m-2)

# Uncertainty of chlo2 carb conversion, coupled with chlo surface to chlo tot conversion
carb_tot_error <- 2.3 # 2.3 fold 95% prediction interval for carb to chlorophyll relationship from Maranon et al., 2014

# Break up total carbon into pico, nano and micro size classes, calculate total in each 1x1 degree grid square
# and sum to get global total, multiply by 10 to convert from carbon to wet weight, then divide by 1e18 to go from mg to Pg
pico_ave <- sum(pico*carb_tot*areas, na.rm = TRUE)*10/1e18
nano_ave <- sum(nano*carb_tot*areas, na.rm = TRUE)*10/1e18
micro_ave <- sum(micro*carb_tot*areas, na.rm = TRUE)*10/1e18
total_error <- chlo_tot_error*carb_tot_error

##########################################################
### STEP 4: CALCULATE GLOBAL BIOMASS OF PICO, NANO AND 
###         MICRO PHYTOPLANKTON, ACROSS INDIVIDUAL CELL
##          SIZE RANGES
##########################################################

# To calculate size bins, convert from ESD to wet weight, assuming 1 g wet weight = 1cm3 biovolume
# Assuming an equal distribution of biomass in log-equal size bins, for example, 
# pico size range is -14.4 to -11.4, total size range is 3 orders of magnitude,
# so 13.3% in -14.4 to -14, 33.3% in -14 to -13, 33.3% in -13 to -12, 20% in -12--11.4 (since minimum size is -14, all biomass for smaller pico goes into -14 to -13 bin)
# nano size range is -11.4 to -8.4, total size range is 3 orders of magnitude,
# so 13.3% in -11.4 to -11, 33.3% in -11 to -10, 33.3% in -10 to -9, 20% in -9--8.4
# micro size range is -8.4 to -5.4, total size range is 3 orders of magnitude,
# so 13.3% in -8.4 to -8, 33.3% in -8 to -7, 33.3% in -7 to -6, 20% in -6--5.4

# Calculate and save total phyto biomass, with uncertainty range
pico_biom <- c(pico_ave, pico_ave/total_error, pico_ave*total_error)
nano_biom <- c(nano_ave, nano_ave/total_error, nano_ave*total_error)
micro_biom <- c(micro_ave, micro_ave/total_error, micro_ave*total_error)

summary_biomass_table <- read.csv("./summary_output/summary_biomass_table.csv")
summary_biomass_table[2, c(2:4)] <- pico_biom
summary_biomass_table[3, c(2:4)] <- nano_biom
summary_biomass_table[4, c(2:4)] <- micro_biom
write.csv(summary_biomass_table, file = "./summary_output/summary_biomass_table.csv", row.names = FALSE)

summary_biomass_table_top200 <- read.csv("./summary_output/summary_biomass_table_top200.csv")
summary_biomass_table_top200[2, c(2:4)] <- pico_biom
summary_biomass_table_top200[3, c(2:4)] <- nano_biom
summary_biomass_table_top200[4, c(2:4)] <- micro_biom
write.csv(summary_biomass_table_top200, file = "./summary_output/summary_biomass_table_top200.csv", row.names = FALSE)

# Plankton size bins
pico_sizes <- c(-13.5,-12.5,-11.5)
nano_sizes <- c(-11.5,-10.5,-9.5,-8.5)
micro_sizes <- c(-8.5,-7.5,-6.5,-5.5)

# Proportion of each group in size bins
pico_props <- c(0.466,0.333,0.20) 
nano_props <- c(0.133,0.333,0.333,0.2)
micro_props <- c(0.133,0.333,0.333,0.2)

# Biomass for each phyto group in size bins
tot_pico_biom_bybin <- t(pico_biom %*% t(pico_props))
tot_nano_biom_bybin <- t(nano_biom %*% t(nano_props))
tot_micro_biom_bybin <- t(micro_biom %*% t(micro_props))

## Save uncertainty estimate
group_standard_errors <- read.csv("./summary_output/group_standard_errors.csv")
total_error <- chlo_tot_error*carb_tot_error
group_standard_errors[group_standard_errors$Group %in% c("Picophytoplankton","Nanophytoplankton","Microphytoplankton"), "log10_Standard_Error"] <- total_error
write.csv(group_standard_errors, file = "./summary_output/group_standard_errors.csv", row.names = FALSE)

## Save total and top 200m biomass estimate, across body size range
summary_biomass_table_long <- read.csv("./summary_output/summary_biomass_table_long.csv")
summary_biomass_table_long[c(which(summary_biomass_table_long$Group == "Picophytoplankton" & 
                                     summary_biomass_table_long$log10_Size_Midpoint_g %in% pico_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_pico_biom_bybin
summary_biomass_table_long[c(which(summary_biomass_table_long$Group == "Nanophytoplankton" & 
                                     summary_biomass_table_long$log10_Size_Midpoint_g %in% nano_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_nano_biom_bybin
summary_biomass_table_long[c(which(summary_biomass_table_long$Group == "Microphytoplankton" & 
                                     summary_biomass_table_long$log10_Size_Midpoint_g %in% micro_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_micro_biom_bybin
write.csv(summary_biomass_table_long, file = "./summary_output/summary_biomass_table_long.csv", row.names = FALSE)

summary_biomass_top200_table_long <- read.csv("./summary_output/summary_biomass_top200_table_long.csv")
summary_biomass_top200_table_long[c(which(summary_biomass_top200_table_long$Group == "Picophytoplankton" & 
                                            summary_biomass_top200_table_long$log10_Size_Midpoint_g %in% pico_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_pico_biom_bybin
summary_biomass_top200_table_long[c(which(summary_biomass_top200_table_long$Group == "Nanophytoplankton" & 
                                            summary_biomass_top200_table_long$log10_Size_Midpoint_g %in% nano_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_nano_biom_bybin
summary_biomass_top200_table_long[c(which(summary_biomass_top200_table_long$Group == "Microphytoplankton" & 
                                            summary_biomass_top200_table_long$log10_Size_Midpoint_g %in% micro_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_micro_biom_bybin
write.csv(summary_biomass_top200_table_long, file = "./summary_output/summary_biomass_top200_table_long.csv", row.names = FALSE)


### save phytoplankton predictions to csv
pred_dat <- read.csv('./data/prediction_data.csv')

### CALCULATE GLOBAL BIOMASS MAPS, FOR PICO, NANO, MICRO AND TOTAL CONVERT FROM MG C TO G WET WEIGHT
### THEN SHIFT MATRIX TO MATCH LAT LONG FROM PREDICTION FILE
pico_global_density <- (pico*carb_tot)/100
pico_global_density <- (pico_global_density[c(181:360,1:180),order(ncol(pico_global_density):1)])

nano_global_density <- (nano*carb_tot)/100
nano_global_density <- (nano_global_density[c(181:360,1:180),order(ncol(nano_global_density):1)])

micro_global_density <- (micro*carb_tot)/100
micro_global_density <- (micro_global_density[c(181:360,1:180),order(ncol(micro_global_density):1)])

global_density <- pico_global_density + nano_global_density + micro_global_density

phyto_predictions <- data.frame('lon' = pred_dat$Long, 
                                'lat' = pred_dat$Lat, 
                                'pico_phyto_biom_gm2_total' = as.vector(pico_global_density),
                                'nano_phyto_biom_gm2_total' = as.vector(nano_global_density),
                                'micro_phyto_biom_gm2_total' = as.vector(micro_global_density),
                                'all_phyto_biom_gm2_total' = as.vector(global_density))


write.csv(phyto_predictions, './global_map_data/phyto_predictions.csv', row.names = FALSE)
