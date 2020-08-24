# sheldon_revisited
Authors: Ian Hatton, Ryan Heneghan, Yinon Bar-On and Eric Galbraith

**Input data to calculate global estimates and uncertainty can be found in the ‘./data’ folder.

**Scripts to estimate biomass for each group can be found in the './biomass_estimate_scripts/' folder

Step 1: Run “./processing_scripts/create_summary_output_dataframes.R” to create the empty data frames to save biomass estimate and uncertainty outputs

Step 2: Run the scripts for each group in './biomass_estimate_scripts/ to calculate biomass estimates and uncertainty: “bacteria_biomass_estimate.R”, “phytoplankton_biomass_estimate.R”,
“nano-microzoo_biomass_estimate.R”, “mesozoo_biomass_estimate.R”, “macrozoo_biomass_estimate.R”, “fish_mammal_other_biomass_estimate.R”

Step 3: Run slope sensitivity test to assess the sensitivity of the global biomass spectrum slope to uncertainty in biomass and distribution of each group's biomass within their respective size range. Test is run through the script './processing_scripts/slope_sensitivity_tests.R'

Other: Climate change projections for bacteria can be obtained by running the code at the bottom of the "bacteria_biomass_estimate.R" script.

***Figures are saved in ./figures, global map data for each group in ‘./global_map_data’, processing scripts for plots and global map data in ‘./scripts’, summary output files in ‘./summary_output’
