#########################################################################
## Script for Sheldon Revisited group biomass estimates
## Authors: Ian Hatton, Ryan Heneghan, Yinon Bar-On and Eric Galbraith
##
############################ BACTERIA ###################################
##
## THIS SCRIPT TAKES RAW DATA FOR BACTERIA AND CALCULATES A
## GLOBAL BIOMASS ESTIMATE, AND AN ESTIMATE OF UNCERTAINTY.
## IT ALSO RETURNS THE INTERPOLATED BIOMASS OF BACTERIA IN 1X1 GRID SQUARES
## FOR THE ENTIRE GLOBAL OCEAN
##
## CLIMATE CHANGE PROJECTIONS FOR BACTERIA ARE OBTAINED AT THE BOTTOM OF
## THIS SCRIPT
#########################################################################
rm(list=ls())
# To plot figures of main effects and partial residuals
library(ggplot2) 
library(egg) 
library(visreg)
library(raster)

# To fit lognormal distribution to individual bacteria biomass data
library(MASS)

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
## STEP 1: GLOBAL BACTERIA ABUNDANCE
##########################################################

## Fit statistical model to bacteria abundance data
bac_abunds <- read.csv('./data/bacteria/bacteria_abundance_data.csv') # Import bacteria abundance data

gm1 <- lm(logabund ~ poly(logdepth,4) + logchlo + SST, data = bac_abunds) # Statistical model linking bacteria abundance with environmental variables
summary(gm1)

# Plot main effects with partial residuals
ylabb = c(expression(paste("log"[10], "(Abundance, ml"^-1, ")"), sep  = ""))

plot_list = list()
plot_list[[1]] = resid_plot(gm1, xvar = "SST", ylab = ylabb, xlab = expression(paste("Temperature (", degree, "C)", sep = "")), subtitle = "a)", new_xticks = NA)
plot_list[[2]] = resid_plot(gm1, xvar = "logchlo", ylab = ylabb, xlab = expression(paste("log"[10], "(Chlorophyll, mg m"^-3, ")")), subtitle = "b)", new_xticks = NA)
plot_list[[3]] = resid_plot(gm1, xvar = "logdepth", ylab = ylabb, xlab = expression(paste("log"[10], "(Depth, m)", sep = "")), subtitle = "c)", new_xticks = NA)

ggsave(filename = "./figures/suppfig1_residuals_bacteria_abundance.png", plot = ggarrange(plots = plot_list, nrow = 1), width = 20, height = 6)

# Calculate global bacteria abundance across entire water column, and top 200m
glob_bac_abund1 <-pred_func(gm1, xvars = c("SST", "logchlo", "logdepth"), pred_dat = pred_dat, grid_areas = areas, other_xvar_name = NA, other_xvar_val = NA, over_depth = TRUE, int = 'confidence')
glob_bac_abund <- sum(glob_bac_abund1$total, na.rm = TRUE)*1e6 # Multiply by 1e6 to go from ml to m^3
glob_bac_abund_top2001 <- pred_func(gm1, xvars = c("SST", "logchlo", "logdepth"), pred_dat = pred_dat, grid_areas = areas, other_xvar_name = NA, other_xvar_val = NA, over_depth = TRUE, int = 'confidence')
glob_bac_abund_top200 <- sum(glob_bac_abund_top2001$total_top200, na.rm = TRUE)*1e6

## Calculate uncertainty of global bacteria abundance estimate, using bootstrap sampling and prediction method n times
n <- 1000 # 1000 times takes about 10 minutes, you can just run the next line to get the result, then skip the loop
abund_uncert <- 3.391625

boot_save <- rep(0, n)
pb <- txtProgressBar(min = 0, max = n, style = 3)
for(i in 1:n){
  setTxtProgressBar(pb, i)
  boot_choose <- sample(1:dim(bac_abunds)[1], dim(bac_abunds)[1], replace = TRUE)
  boot_dat <- bac_abunds[boot_choose,]
  gm1 <- lm(logabund ~ poly(logdepth,4) + logchlo + SST, data = boot_dat) # Statistical model linking bacteria abundance with environmental variables
  bac_abund_preds <- pred_func(gm1, xvars = c("SST", "logchlo", "logdepth"), pred_dat = pred_dat, grid_areas = areas, other_xvar_name = NA, other_xvar_val = NA, over_depth = TRUE, int = 'confidence')
  boot_save[i] <- sum(bac_abund_preds$total, na.rm = TRUE)
}

abund_uncert <- 1.96*(10^(0.5*(sd(bac_abunds$logabund) + sd(log10(boot_save)))))

##########################################################
## STEP 2: GLOBAL BACTERIA BIOMASS, ESTIMATE AND SAVE
##         GLOBAL MAP, AGGREGATE AND SPLIT BIOMASS ESTIMATES
##########################################################

## Import individual bacteria biomass data
bac_bioms <- read.csv("./data/bacteria/bacteria_ind_bioms.csv")

## Fit lognormal distribution and plot histogram of all in situ individual bacteria biomass data together
fit_params1 <- fitdistr(10^(bac_bioms$log10_Vol_um3),"lognormal")

bac_biom_hist <- ggplot(bac_bioms, aes(log10_Vol_um3)) + 
  geom_histogram(aes(y=..density..),binwidth = .1, col="black", fill = "darkgray",size=0.5,position="identity")+
  stat_function(fun = dnorm, args = list(mean = log10(exp(fit_params1$estimate[1])), sd = log10(exp(fit_params1$estimate[2]))),
                col = "darkred", size = 1)+
  theme_bw()+theme(axis.text = element_text(size = 24),
                   plot.subtitle = element_text(size = 24),
                   axis.title = element_text(size = 24),
                   panel.border = element_rect(size = 1),
                   axis.ticks = element_line(size = 1),
                   legend.position = "none")+
  ylab('Density') + 
  xlab(expression(paste('log'[10],'(Body Size ', mu, 'g)', sep = '')))

ggsave(filename = "./figures/suppfig2_density_bacteria_abundance.png", plot =bac_biom_hist, width = 8, height =7)

## Calculate total bacteria biomass and 95% uncertainty bounds
tot_bac <- c(glob_bac_abund, glob_bac_abund/abund_uncert, glob_bac_abund*abund_uncert) # 95% interval and mean global bacteria abundance estimate
tot_bac_top200 <- c(glob_bac_abund_top200, glob_bac_abund_top200/abund_uncert, glob_bac_abund_top200*abund_uncert) # 95% interval and mean global bacteria abundance estimate in top 200m
mean_distn <- c(log10(exp(fit_params1$estimate[1])),log10(exp(fit_params1$estimate[1] - 1.96*fit_params1$sd[1])),
                log10(exp(fit_params1$estimate[1] + 1.96*fit_params1$sd[1]))) # 95% confidence interval bounds and average of the mean of the lognormal size distribution of individual bacteria biomass

## Predict global bacteria abundance, in 360x180 grid squares
bacteria_predictions <- data.frame('lon' = pred_dat$Long, 'lat' = pred_dat$Lat, 'bacteria_abund_m2_alldepths' = glob_bac_abund1$density*1e6,
                                   'bacteria_biom_gm2_top200m' = glob_bac_abund_top2001$density_top200*10^(mean_distn[1])/1e6, 
                                   'bacteria_biom_gm2_alldepths' = glob_bac_abund1$density*10^(mean_distn[1])/1e6)

## Save bacteria predictions to csv
write.csv(bacteria_predictions, './global_map_data/bacteria_predictions.csv', row.names = FALSE)

## Save uncertainty estimate
group_standard_errors <- read.csv("./summary_output/group_standard_errors.csv")
tot_bac_uncert <- abund_uncert*exp(1.96*fit_params1$sd[1])
group_standard_errors[group_standard_errors$Group == "Hetero_Bacteria", "log10_Standard_Error"] <- tot_bac_uncert
write.csv(group_standard_errors, file = "./summary_output/group_standard_errors.csv", row.names = FALSE)

## Save total and top 200m biomass estimate, across body size range
tot_bac_biom <- tot_bac*10^(mean_distn)/1e27 # convert ug to Pg
tot_bac_biom_upper200 <- tot_bac_top200*10^(mean_distn)/1e27 # convert ug to Pg

summary_biomass_table <- read.csv("./summary_output/summary_biomass_table.csv")
summary_biomass_table[1, c(2:4)] <- tot_bac_biom
write.csv(summary_biomass_table, file = "./summary_output/summary_biomass_table.csv", row.names = FALSE)

summary_biomass_table_top200 <- read.csv("./summary_output/summary_biomass_table_top200.csv")
summary_biomass_table_top200[1, c(2:4)] <- tot_bac_biom_upper200
write.csv(summary_biomass_table_top200, file = "./summary_output/summary_biomass_table_top200.csv", row.names = FALSE)

bac_sizes <- c(-13.5, -12.5, -11.5) # Size bins bacteria biomass falls into (log10 g)
prop_biom <- c(0.4, 0.4, 0.2) # Proportion of total bacteria biomass in each size bin (size range of bacteria is 10^-14 to 10^-11.5g, split evenly across log10 size bins)

tot_bac_biom_bybin <- t(tot_bac_biom %*% t(prop_biom))
tot_bac_biom_upper200_bybin <- t(tot_bac_biom_upper200 %*% t(prop_biom))

summary_biomass_table_long <- read.csv("./summary_output/summary_biomass_table_long.csv")
summary_biomass_table_long[c(which(summary_biomass_table_long$Group == "Hetero_Bacteria" & 
                                     summary_biomass_table_long$log10_Size_Midpoint_g %in% bac_sizes)),
                           c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_bac_biom_bybin
write.csv(summary_biomass_table_long, file = "./summary_output/summary_biomass_table_long.csv", row.names = FALSE)

summary_biomass_top200_table_long <- read.csv("./summary_output/summary_biomass_top200_table_long.csv")
summary_biomass_top200_table_long[c(which(summary_biomass_top200_table_long$Group == "Hetero_Bacteria" & 
                                            summary_biomass_top200_table_long$log10_Size_Midpoint_g %in% bac_sizes)),
                                  c("Biomass_Pg_wet_weight_estimate",  "Biomass_Pg_wet_weight_95CI_lower", "Biomass_Pg_wet_weight_95CI_upper")] <- tot_bac_biom_upper200_bybin
write.csv(summary_biomass_top200_table_long, file = "./summary_output/summary_biomass_top200_table_long.csv", row.names = FALSE)





#########################################################################
################## CLIMATE CHANGE BACTERIA PROJECTIONS ###################
##
##### PREDICT CHANGE IN GLOBAL BACTERIA ABUNDANCE IN 1950-1960, 2020-2030 AND 
##### 2090-2100 USING IPSL AND GFDL ESM MODEL PROJECTIONS
#########################################################################

pred_dat1 <- read.csv('./data/prediction_data.csv')
ipsl_dat <- read.csv("./data/climate_change_forcings/ipsl_clim_forcings.csv")
gfdl_dat <- read.csv("./data/climate_change_forcings/gfdl_clim_forcings.csv")

num_store <- matrix(0, nrow = 3, ncol = 3)
rownames(num_store) <- c("IPSL", "GFDL","Mean")
colnames(num_store) <- c("1950-1960", "2020-2030", "2090-2100")

for(i in 1:2){
  if(i == 1){ #IPSL ESTIMATES
    for(j in 1:3){
      pred_dat <- pred_dat1
      pred_dat$SST <- ipsl_dat[,((j+1)*2)]
      pred_dat$logchlo <- log10(ipsl_dat[,((j+1)*2-1)])
      ppred <- pred_func(gm1, xvars = c("SST", "logchlo", "logdepth"), pred_dat = pred_dat, grid_areas = areas, other_xvar_name = NA, other_xvar_val = NA, over_depth = TRUE, int = 'confidence')
      num_store[i,j] <- sum(ppred$total, na.rm = TRUE)
    }
  }
  
  if(i == 2){ # GFDL ESTIMATES
    for(j in 1:3){
      pred_dat <- pred_dat1
      
      pred_dat$SST <- gfdl_dat[,((j+1)*2)]
      pred_dat$logchlo <- log10(gfdl_dat[,((j+1)*2-1)])
      ppred <-pred_func(gm1, xvars = c("SST", "logchlo", "logdepth"), pred_dat = pred_dat, grid_areas = areas, other_xvar_name = NA, other_xvar_val = NA, over_depth = TRUE, int = 'confidence')
      num_store[i,j] <- sum(ppred$total, na.rm = TRUE)
    }
  }
}

num_store[1,] <- num_store[1,]/num_store[1,1]
num_store[2,] <- num_store[2,]/num_store[2,1]
num_store[3,] <- colMeans(num_store[1:2,])
num_store

write.csv(num_store, "./summary_output/bacteria_climchange_projections", row.names = FALSE)


