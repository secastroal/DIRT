# Dynamic Item Response Theory (DIRT)

## Description
This project aims to develop polytomous item response theory models suitable to analyze intensive longitudinal psychological data. In particular, we are extending the partial credit model in order to analyze multivariate time series data (i.e., N = 1).

## Folders and Files Description
This repository contains the main files needed to fit the dynamic item response theory models that we have developed and tested. Output files from the simulation studies (preliminary and final) are also available. The most important files and folders are: 
* R: This folder contains R files with custom functions.
  + IRT_models.R: Functions to simulate standard IRT model such as the 1PL, 2PL, 3PL, PCM, and GRM.
  + PPC.R: Functions to compute posterior predictive model checking methods of the AR-GRM and the AR-PCM models.
* Simulation: This folder contains the output files from the simulation studies. These files are read and analyzed in [Results_simulation.R](Results_simulation.R).
* Stan: This folder contains diverse stan models.
  + [ar_irt_pcm_na.stan](Stan/ar_irt_pcm_na.stan): This is the most updated AR-PCM model. This stan model can handle missing data.
  + [ar_irt_na.stan](Stan/ar_irt_na.stan): This is the most updated AR-GRM model. This stan model can handle missing data.
* [AR-IRT-PCM-NA.R](AR-IRT-PCM-NA.R): R script to simulate data based on the AR-PCM model and to fit the model as it is in [ar_irt_pcm_na.stan](Stan/ar_irt_pcm_na.stan). 
* [AR-IRT-NA.R](AR-IRT-NA.R): R script to simulate data based on the AR-GRM model and to fit the model as it is in [ar_irt_na.stan](Stan/ar_irt_na.stan).
* [Results_simulation.R](Results_simulation.R): R script to summarize and plot the output from the simulations of the AR-PCM and the AR-GRM.
* [Sim_AR_PCM.R](Sim_AR_PCM.R): R script of the simulation study to test the performance of the AR-PCM across 72 conditions.
* [Sim_AR_GRM.R](Sim_AR_GRM.R): R script of the simulation study to test the performance of the AR-GRM across 72 conditions.
