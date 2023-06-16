This Repo is a compilation of ongoing work by Sitara Baboolal on Combined Sewage Overflow Design Optimization.

 This repository contains codes and data that are a work in progress.

Keywords:
Extreme Precipitation, Flood Risk, Stormwater Infrastructure Design, 
Climate Change, Decision-Making under Deep Uncertainty

Question? Sitara Baboolal sitara.baboolal08@gmail.com
    
Project Summary
--------------
* CSO Optimization -
* PYSWMM
* PipeDesign-master - work done by Sharma et al 2021 used as a base case for the CSO Optimization study
*
*
* "dataset" - Directory containing all data needed for simulation and to generate figures used in Sharma et al. 2021. ** (Directory needs seperating into input and output data to be less confusing) **
  * "cordex" - Directory containing dynamically downscalled Climate data from 6 Regional Climate Models
  * "maca" - Directory containing statstically downscalled Climate data from 3 Global Climate Models
* "projections" - Directory containing code to generate climate projections provided in "dataset"
* "SourceCodes" - Directory containing functions used for climate projections and SOF calculation for Pipe diameter
  * "Prior2SourceMu.R" - non-Stationary
  * "Prior2SourceStat.R" - Stationalry
  * "batchmeans.R"
  * "failureprob_diameter.R" - Reliability function for pipe diameter calculation
  
* "Figure" - Directory containing code to generate visualizations from simulation results presented in Sharma et al. 2021. 
  * "design_pipe.R"
  * "precip_estimates.R"
  * "precip_projections.R" - 
* "SampleFigures" - Directory containing sample figures presented in Sharma et al. 2021
  


