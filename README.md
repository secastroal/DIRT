# Time-Varying Dynamic Partial Credit Model (TV-DPCM)

## Description
In this project, we developed an item response theory model suitable to analyze intensive longitudinal data in psychology. In particular, we extended the partial credit model in order to analyze multivariate time series data (i.e., N = 1). The proposed model combines the partial credit model (PCM) and the time-varying autoregressive model (TV-AR as presented in [Bringmann et al., 2017](https://psycnet.apa.org/doiLanding?doi=10.1037%2Fmet0000085)). In a nutshell, the TV-DPCM is a measurement model for intensive longitudinal data that also allows modeling trend-stationary time series when non-linear trends are present. We refer to this model as the time-varying dynamic partial credit model (TV-DPCM). A detailed description of the TV-DPCM can be found in [Castro-Alvarez et al. (2023a)](https://www.tandfonline.com/doi/full/10.1080/00273171.2023.2214787).

Additionally, as the model was implemented within the Bayesian framework in Stan, we also developed different test statistics and discrepancy measures based on the posterior predictive model checking method for the TV-DPCM. A preprint introducing these fit statistics can be found in [Castro-Alvarez et al. (2023b)](https://osf.io/preprints/psyarxiv/kufnh/).

## Folders and Files Description
This repository contains the main files needed to fit the TV-DPCM, as well as, code needed to run the simulation studies performed for each paper. Output files from the simulation studies (preliminary and final) are also available. The most important files and folders are: 
* R: This folder contains R files with custom functions.
  + IRT_models.R: Functions to simulate standard IRT model such as the 1PL, 2PL, 3PL, PCM, and GRM.
  + IRT_plots.R: Functions to plot the IICs, IIFs, and TIF of the TV-DPCM.
  + PPMC.R: Functions to compute posterior predictive model checking methods for the TV-DPCM.
  + gen.TVDPCM.R: Function to simulate data based on the TV-DPCM.
  + tvdpcm2stan.R: Function to turn data into a list as required by Stan.
* Simulation: This folder contains the output files from the simulation studies. These files are read and analyzed in [Results_simulation_tv_dpcm.R](Results_simulation_tv_dpcm.R) and [Results_ppmc.R](Results_ppmc.R).
* Stan: This folder contains diverse stan models.
  + [tv_dpcm_int_v5.1.stan](Stan/tv_dpcm_int_v5.1.stan): This is the stan TV-DPCM model used in the simulations and in the empirical examples from the two papers.
* [ManuscriptPlots.R](ManuscriptPlots.R): R script used to create some of the figures included in [Castro-Alvarez et al. (2023a)](https://www.tandfonline.com/doi/full/10.1080/00273171.2023.2214787).
* [PPMCManuscriptPlots.R](PPMCManuscriptPlots.R): R script used to create some of the figures included in [Castro-Alvarez et al. (2023b)](https://osf.io/preprints/psyarxiv/kufnh/).
* [Results_ppmc.R](Results_ppmc.R): R script to summarize and plot the output from the simulation study presented in [Castro-Alvarez et al. (2023b)](https://osf.io/preprints/psyarxiv/kufnh/).
* [Results_simulation_tv_dpcm.R](Results_simulation_tv_dpcm.R): R script to summarize and plot the output from the simulation study presented in [Castro-Alvarez et al. (2023a)](https://www.tandfonline.com/doi/full/10.1080/00273171.2023.2214787).
* [Sim_PPMC.R](Sim_PPMC.R): R script to run the simulation study, which was used to assess the performance of the different measures proposed within the PPMC method for the TV-DPCM as presented in [Castro-Alvarez et al. (2023b)](https://osf.io/preprints/psyarxiv/kufnh/).  
* [Sim_TV_DPCM.R](Sim_TV_DPCM.R): R script to run the simulation study, which was used to test the performance of the TV-DPCM as presented in [Castro-Alvarez et al. (2023a)](https://www.tandfonline.com/doi/full/10.1080/00273171.2023.2214787). 
* Rmarkdown: This folder contains Rmarkdown files used to fit the model to empirical data, to run simulation tests, and to summarize the results of preliminary simulations.
  + [PGdata_analysis_SE.Rmd](Rmarkdown/PGdata_analysis_SE.Rmd): Includes code to run the analysis of the empirical data as presented in [Castro-Alvarez et al. (2023a)](https://www.tandfonline.com/doi/full/10.1080/00273171.2023.2214787). 
  + [PGdata_PPMC.Rmd](PGdata_PPMC.Rmd): Includes code to run the analysis of the empirical data as presented in [Castro-Alvarez et al. (2023b)](https://osf.io/preprints/psyarxiv/kufnh/).
  
## Notes
The data used for the empirical examples in both papers was retrieved and is openly available at [Kossakowski et al. (2017)](http://doi.org/10.5334/jopd.29).

The simulation studies performed within this research were run on [peregrine](https://www.rug.nl/society-business/centre-for-information-technology/research/services/hpc/facilities/peregrine-hpc-cluster), the High Performance Computing cluster available at the University of Groningen.

