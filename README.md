# Dynamic Item Response Theory (DIRT)

## Description
This project aims to develop polytomous item response theory models suitable to analyze intensive longitudinal psychological data. In particular, we are extending the partial credit model and the graded response model in order to analyze multivariate time series data (i.e., N = 1).

## Folders and Files Description
This repository contains the main files needed to fit the dynamic item response theory models that we have developed and tested. The most important files and folders are: 
* R: This folder contains R files with custom functions.
+ IRT_models.R: Functions to simulate standard IRT model such as the 1PL, 2PL, 3PL, PCM, and GRM.
+ PPC.R: Functions to compute posterior predictive model checking methods of the AR-GRM and the AR-PCM models.
* Stan: This folder contains the diverse stan models.
+ ar_irt_pcm_na.stan: This is the most updated AR-PCM model. This stan model can handle missing data.
+ ar_irt_na.stan: This is the most updated AR-GRM model. This stan model can handle missing data.
* AR-IRT-PCM-NA.R: This file simulates data based on the AR-PCM model and fits the model as it is in [ar_irt_pcm_na.stan](Stan/ar_irt_pcm_na.stan). 
* AR-IRT-NA.R: This file simulates data based on the AR-GRM model and fits the model as it is in [ar_irt_na.stan](Stan/ar_irt_na.stan).
