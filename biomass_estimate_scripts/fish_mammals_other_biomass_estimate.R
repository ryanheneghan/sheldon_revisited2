#################################################################################
## Script for Sheldon Revisited group biomass estimates
## Authors: Ian Hatton, Ryan Heneghan, Yinon Bar-On, Eric Galbraith
##
## MESOPELAGICS AND NON-MESOPELAGIC FISH, CEPHALOPODS, INVEREBRATES AND MAMMALS #
##
## THIS SCRIPT TAKES RAW DATA FOR MESOPELAGICS, NON-MESOPELAGICS FISH AND INVEREBRATES
##  AND MAMMALS AND CALCULATES GLOBAL BIOMASS ESTIMATES, AND ESTIMATES OF UNCERTAINTY.
## IT ALSO RETURNS THE INTERPOLATED BIOMASS OF NON-MESOPELAGIC FISH 
## IN 1X1 GRID SQUARES FOR THE ENTIRE GLOBAL OCEAN
#################################################################################


rm(list=ls())

# To open fish netcdfs
library(ncdf4)

# For area of 360x180 global grid squares
library(raster)

# Set working directory to scripts_data_figures folder
#setwd("~/Desktop/Papers/Sheldon_Revisited/scripts_data_figures/")

# Source area file of 360x180 global grid (in m2) for global biomass/abundance calculation
area_global1 <- t(as.matrix(raster::area(raster()))*1000*1000) # Area of 360x180 global grid squares in m^2

##########################################################
## MESOPELAGIC FISH BIOMASS
##########################################################
 # Global mesopelagic fish biomass range from Proud et al., 2019 (doi:10.1093/icesjms/fsy037)
meso_fish_range <- c(1.8, 16) # Pg wet weight
meso_fish_estimate <- 10^(mean(log10(meso_fish_range)))
meso_fish_uncert <- meso_fish_range[2]/meso_fish_estimate
meso_fish_biom_range <- c(meso_fish_estimate, meso_fish_estimate/meso_fish_uncert, meso_fish_estimate*meso_fish_uncert)

## Save uncertainty estimate
group_standard_errors <- read.csv("./summary_output/group_standard_errors.csv")
group_standard_errors[group_standard_errors$Group == "Mesopelagics", "log10_Standard_Error"] <- meso_fish_uncert
write.csv(group_standard_errors, file = "./summary_output/group_standard_errors.csv", row.names = FALSE)

summary_biomass_table <- read.csv("./summary_output/summary_biomass_table.csv")
summary_biomass_table[9, c(2:4)] <- meso_fish_biom_range
write.csv(summary_biomass_table, file = "./summary_output/summary_biomass_table.csv", row.names = FALSE)



# Mesopelagic fish size range of 0.01gm - 500gm from Eduardo et al., 2018; Escobar-Flores et al., 2020; Lopez-Perez et al., 2020
meso_fish_size_range <- c(-1.5, -0.5, 0.5, 1.5, 2.5)
meso_fish_size_props <- c(0.22,0.22,0.22,0.22,0.12)

tot_meso_fish_biom_bybin <- t(meso_fish_biom_range %*% t(meso_fish_size_props))

summary_biomass_table_long <- read.csv("./summary_output/summary_biomass_table_long.csv")
summary_biomass_table_long[c(which(summary_biomass_table_long$Group == "Mesopelagics" & 
                                     summary_biomass_table_long$log10_Size_Midpoint_g %in%  meso_fish_size_range)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_meso_fish_biom_bybin
write.csv(summary_biomass_table_long, file = "./summary_output/summary_biomass_table_long.csv", row.names = FALSE)


##################################################################
## EPIPELAGIC, DEMERSAL FISH, CEPHALOPDS AND INVERTEBRATE BIOMASS
##################################################################

### For each marine model, extract average map of global fish biomass from 1850-1960
# Process CESM
boats_cesm_tcb_1850_1860 <- apply(ncvar_get(nc_open(list.files(path = "./data/fish_mammals_other/",pattern = glob2rx("boats*cesm*hist*tcb*1850*"), recursive = TRUE, full.names = TRUE)), 'tcb')[,,1:120], c(1,2), mean, na.rm = TRUE)
boats_cesm_tcb_1850_1860 <- boats_cesm_tcb_1850_1860[c(181:360,1:180),]

### Process FEISTY
# Process CESM
feisty_cesm_tcb_1850_1860 <- apply(ncvar_get(nc_open(list.files(path =  "./data/fish_mammals_other/",pattern = glob2rx("feisty*cesm*hist*tcb*1850*"), recursive = TRUE, full.names = TRUE)), 'tcb')[,,1:10], c(1,2), mean, na.rm = TRUE)
feisty_cesm_b10_1850_1860 <- apply(ncvar_get(nc_open(list.files(path =  "./data/fish_mammals_other/",pattern = glob2rx("feisty*cesm*hist*b10*1850*"), recursive = TRUE, full.names = TRUE)), 'b10cm')[,,1:10], c(1,2), mean, na.rm = TRUE)

####### CREATE ARRAYS OF ALL AVERAGE MAPS FOR 1950-1960, 2020-2030 AND 2090-100

hist_1850_1860 <- array(c(boats_cesm_tcb_1850_1860,feisty_cesm_tcb_1850_1860, feisty_cesm_b10_1850_1860), dim = c(360,180,3))

hist_1850_1860[which(hist_1850_1860 == 0)] <- NA ## All land values are NA

# 1850_1860 global biomass average (Pg wet weight)
ave_hist_1850_1860_all <- apply(sweep(hist_1850_1860, c(1,2), area_global1, '*'), c(3), sum, na.rm = TRUE)*10/1e15

# < 10cm biomass estimate
small_fish <- ave_hist_1850_1860_all[2] - ave_hist_1850_1860_all[3]
small_fish_pred <- as.vector(hist_1850_1860[,,2]) - as.vector(hist_1850_1860[,,3])

# > 10cm fish biomass estimate
large_fish <- 10^((log10(ave_hist_1850_1860_all[1]) + log10(ave_hist_1850_1860_all[3]))/2)
large_fish_pred <- as.vector(10^((log10(hist_1850_1860[,,1])+log10(hist_1850_1860[,,3]))/2))

# Total fish biomass
total_fish <- small_fish + large_fish
total_fish_pred <- small_fish_pred + large_fish_pred

# Fish uncertainty, from Jennings and Collingridge, 2015 (doi:10.1371/journal.pone.0133794)
fish_uncert <- 5 

## Save uncertainty estimate
group_standard_errors <- read.csv("./summary_output/group_standard_errors.csv")
group_standard_errors[group_standard_errors$Group == "Fish", "log10_Standard_Error"] <- fish_uncert
write.csv(group_standard_errors, file = "./summary_output/group_standard_errors.csv", row.names = FALSE)

  
### Save fish biomass, broken into size bins
# Small fish
small_fish_sizes <- c(-2.5,-1.5,-0.5,0.5)
small_fish_props <- c(0.25,0.25,0.25,0.25)

small_fish_bioms <- c(small_fish, small_fish/fish_uncert, small_fish*fish_uncert)
tot_small_fish_biom_bybin <- t(small_fish_bioms %*% t(small_fish_props))

# Large fish
large_fish_sizes <- c(1.5, 2.5, 3.5, 4.5)
large_fish_props <- c(0.25, 0.25, 0.25, 0.25)

large_fish_bioms <- c(large_fish, large_fish/fish_uncert, large_fish*fish_uncert)
tot_large_fish_biom_bybin <- t(large_fish_bioms %*% t(large_fish_props))

summary_biomass_table_long <- read.csv("./summary_output/summary_biomass_table_long.csv")
summary_biomass_table_long[c(which(summary_biomass_table_long$Group == "Fish" & 
                                     summary_biomass_table_long$log10_Size_Midpoint_g %in%  small_fish_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_small_fish_biom_bybin
summary_biomass_table_long[c(which(summary_biomass_table_long$Group == "Fish" & 
                                     summary_biomass_table_long$log10_Size_Midpoint_g %in%  large_fish_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_large_fish_biom_bybin
write.csv(summary_biomass_table_long, file = "./summary_output/summary_biomass_table_long.csv", row.names = FALSE)


summary_biomass_top200_table_long <- read.csv("./summary_output/summary_biomass_top200_table_long.csv")
summary_biomass_top200_table_long[c(which(summary_biomass_top200_table_long$Group == "Fish" & 
                                            summary_biomass_top200_table_long$log10_Size_Midpoint_g %in%  small_fish_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_small_fish_biom_bybin
summary_biomass_top200_table_long[c(which(summary_biomass_top200_table_long$Group == "Fish" & 
                                            summary_biomass_top200_table_long$log10_Size_Midpoint_g %in%  large_fish_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_large_fish_biom_bybin
write.csv(summary_biomass_top200_table_long, file = "./summary_output/summary_biomass_top200_table_long.csv", row.names = FALSE)


### Save predictions of non-meso fish, demersals, cephalopod biomass
dat <- expand.grid("lon" = -179.5:179.5, "lat" = 89.5:-89.5)
dat$fish_biom_gm2_1mg_to_1e5g <- total_fish_pred*10

write.csv(dat, "./global_map_data/fish_predictions.csv", row.names = FALSE)


##########################################################
## MARINE MAMMALS AND VERY LARGE FISH (> 100KG)
##########################################################

marmam <-  read.csv("./data/fish_mammals_other/marine_mammal_biom_data.csv")
mammal_uncert <- 7
gmean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

xb <- marmam$Body_mass_g
yb <- as.numeric(marmam$Pop_pre.whaling_reg)
ybm <- yb * xb /1e15
z <- 10^(4:9)
Z1 <- data.frame("mass"=rep(0,length(z)-1),"y_var"=rep(0,length(z)-1), "ybio"=rep(0,length(z)-1))
for (i in 1:(length(z)-1)){
  Z1[i,1] <- gmean(c(z[i],z[i+1]))
  Z1[i,2] <- sum(yb[which(xb>z[i] & xb<=z[i+1])])
  Z1[i,3] <- sum(ybm[which(xb>z[i] & xb<=z[i+1])])
}
mams <- data.frame("Cat"=rep(5,5),"Group"=rep("Mammals",5),"Body_mass_class_g"=Z1[,1],"Abundance"=Z1[,2],"Mean_Gt"=Z1[,3],"Lo95_Gt"=Z1[,3]/mammal_uncert,"Up95_Gt"=Z1[,3]*mammal_uncert,"Size_Lo"=10^(4:8),"Size_Hi"=10^(5:9))

## Save uncertainty estimate
group_standard_errors <- read.csv("./summary_output/group_standard_errors.csv")
  group_standard_errors[group_standard_errors$Group == "Mammals", "log10_Standard_Error"] <- mammal_uncert
write.csv(group_standard_errors, file = "./summary_output/group_standard_errors.csv", row.names = FALSE)

tot_mam_bioms <- colSums(mams[,c("Mean_Gt", "Lo95_Gt", "Up95_Gt")])

summary_biomass_table <- read.csv("./summary_output/summary_biomass_table.csv")
summary_biomass_table[11, c(2:4)] <- tot_mam_bioms
write.csv(summary_biomass_table, file = "./summary_output/summary_biomass_table.csv", row.names = FALSE)

summary_biomass_table_top200 <- read.csv("./summary_output/summary_biomass_table_top200.csv")
summary_biomass_table_top200 [11, c(2:4)] <- tot_mam_bioms
write.csv(summary_biomass_table_top200 , file = "./summary_output/summary_biomass_table_top200.csv", row.names = FALSE)

log_mam_sizebins <- log10(mams$Body_mass_class_g)

summary_biomass_table_long <- read.csv("./summary_output/summary_biomass_table_long.csv")
summary_biomass_table_long[c(which(summary_biomass_table_long$Group == "Mammals" & 
                                     summary_biomass_table_long$log10_Size_Midpoint_g %in%  log_mam_sizebins)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- mams[,c("Mean_Gt", "Lo95_Gt", "Up95_Gt")]
write.csv(summary_biomass_table_long, file = "./summary_output/summary_biomass_table_long.csv", row.names = FALSE)

summary_biomass_top200_table_long <- read.csv("./summary_output/summary_biomass_top200_table_long.csv")
summary_biomass_top200_table_long[c(which(summary_biomass_top200_table_long$Group == "Mammals" & 
                                            summary_biomass_top200_table_long$log10_Size_Midpoint_g %in%  log_mam_sizebins)),
                                  c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- mams[,c("Mean_Gt", "Lo95_Gt", "Up95_Gt")]
write.csv(summary_biomass_top200_table_long, file = "./summary_output/summary_biomass_top200_table_long.csv", row.names = FALSE)

##################################################
##### VERY LARGE FISH (100kg - 1 tonne, 1 tonne - 10 tonne)
##################################################
# 0.57 Pg for 100kg - 1tonne and 0.48 Pg for 1 tonne to 10 tonnes,
# comes from Jennings and Collingridge 2015 model for total consumer biomass in those size bins.
# Subtract mammal biomass to get fish biomass in those bins

fish_100kg_1tonne <- 0.57 - mams[2,"Mean_Gt"] 
fish_1tonne_10tonne <- 0.48 - mams[3,"Mean_Gt"] # 0.48 Pg comes from Jennings and Collingridge 2015 model

largest_fish_sizes <- c(5.5, 6.5)
largest_fish_bioms <- matrix(c(fish_100kg_1tonne, fish_100kg_1tonne/fish_uncert, fish_100kg_1tonne*fish_uncert,
                               fish_1tonne_10tonne, fish_1tonne_10tonne/fish_uncert, fish_1tonne_10tonne*fish_uncert), byrow = TRUE,
                             nrow = 2, ncol = 3)

summary_biomass_table_long <- read.csv("./summary_output/summary_biomass_table_long.csv")
summary_biomass_table_long[c(which(summary_biomass_table_long$Group == "Fish" & 
                                     summary_biomass_table_long$log10_Size_Midpoint_g %in%  largest_fish_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- largest_fish_bioms
write.csv(summary_biomass_table_long, file = "./summary_output/summary_biomass_table_long.csv", row.names = FALSE)

summary_biomass_top200_table_long <- read.csv("./summary_output/summary_biomass_top200_table_long.csv")
summary_biomass_top200_table_long[c(which(summary_biomass_top200_table_long$Group == "Fish" & 
                                            summary_biomass_top200_table_long$log10_Size_Midpoint_g %in%  largest_fish_sizes)),
                                  c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- largest_fish_bioms
write.csv(summary_biomass_top200_table_long, file = "./summary_output/summary_biomass_top200_table_long.csv", row.names = FALSE)

tot_fish_bioms <- small_fish_bioms + large_fish_bioms + colSums(largest_fish_bioms)

summary_biomass_table_top200 <- read.csv("./summary_output/summary_biomass_table_top200.csv")
summary_biomass_table_top200[10, c(2:4)] <-  tot_fish_bioms 
write.csv(summary_biomass_table_top200 , file = "./summary_output/summary_biomass_table_top200.csv", row.names = FALSE)

summary_biomass_table <- read.csv("./summary_output/summary_biomass_table.csv")
summary_biomass_table[10, c(2:4)] <- tot_fish_bioms 
write.csv(summary_biomass_table, file = "./summary_output/summary_biomass_table.csv", row.names = FALSE)


