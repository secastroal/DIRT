# Time-Varying Dynamic Partial Credit Model (TV-DPCM)

## Description
This project aims to develop polytomous item response theory models suitable to analyze intensive longitudinal psychological data. In particular, we extended the partial credit model in order to analyze multivariate time series data (i.e., N = 1). The proposed model combines the partial credit model (PCM) and the time-varying autoregressive model (TV-AR as presented in [Bringmann et al., 2017](https://psycnet.apa.org/doiLanding?doi=10.1037%2Fmet0000085)). In a nutshell, the TV-DPCM is a measurement model for intensive longitudinal data that also allows modeling trend-stationary time series when non-linear trends are present. A preprint introducing the TV-DPCM is available [here](https://psyarxiv.com/udnbt/).

Additionally, as the model was implemented within the Bayesian framework in Stan, we are also developing posterior predictive model checking methods for the TV-DPCM. These can be found in the file [PPMC.R](R/PPMC.R).

## Folders and Files Description
This repository contains the main files needed to fit the TV-DPCM, as well as, code needed to run the simulation study. Output files from the simulation studies (preliminary and final) are also available. The most important files and folders are: 
* R: This folder contains R files with custom functions.
  + IRT_models.R: Functions to simulate standard IRT model such as the 1PL, 2PL, 3PL, PCM, and GRM.
  + IRT_plots.R: Functions to plot the IICs, IIFs, and TIF of the TV-DPCM.
  + PPMC.R: Functions to compute posterior predictive model checking methods of the TV-DPCM.
  + gen.TVDPCM.R: Function to simulate data based on the TV-DPCM.
  + tvdpcm2stan.R: Function to turn data into a list as required by Stan.
* Simulation: This folder contains the output files from the simulation studies. These files are read and analyzed in [Results_simulation_tv_dpcm.R](Results_simulation_tv_dpcm.R).
* Stan: This folder contains diverse stan models.
  + [tv_dpcm_int_v5.1.stan](Stan/tv_dpcm_int_v5.1.stan): This is the stan TV-DPCM model used in the simulation and in the empirical example.
* [Results_simulation_tv_dpcm.R](Results_simulation_tv_dpcm.R): R script to summarize and plot the output from the simulations of the TV-DPCM.
* [Sim_TV_DPCM.R](Sim_TV_DPCM.R): R script of the simulation study to test the performance of the TV-DPCM. This simulation was run on [peregrine](https://www.rug.nl/society-business/centre-for-information-technology/research/services/hpc/facilities/peregrine-hpc-cluster), the High Performance Computing cluster available at the University of Groningen.
