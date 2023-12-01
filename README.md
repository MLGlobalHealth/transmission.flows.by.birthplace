# Sources of HIV infections among MSM with a migration background:  a viral phylogenetic case study in Amsterdam, the Netherlands

This repository includes code and partial data for the analyses in Blenkinsop, Sofocleous, Kostaki et al. aRxiv (2023)

## Data
The data folder contains the following input files:
* Anonymised transmission pairs including covariates:
    * age of transmitter and recipient at estimated time of transmission of the incident case
    * time elapsed
    * patristic distance
* Data from Belgian study for fitting the evolutionary clock meta-analysis model

## Code
The scripts folder contains the following scripts to run the analysis:
1) run-stan.R - run stan model using cmdstan on anonymised data from formulated pairs
2) post-processing.R - write script to produce convergence diagnostics and generate figures and tables
