# Sources of HIV infections among MSM with a migration background:  a viral phylogenetic case study in Amsterdam, the Netherlands

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![aRxiv link](https://img.shields.io/badge/aRxiv-link%20to%20paper-blue)](https://doi.org/10.48550/arXiv.2401.08308)

This repository includes code and partial data for the analyses in Blenkinsop, Sofocleous, Kostaki et al. *aRxiv* (2024)

- [Data](#floppydisk-data)
- [Code](#computer-code)
- [License](#pagefacingup-license)
- [Warranty](#shield-warranty)
- [Citation](#books-citation)
- [Acknowledgments](#acknowledgments)

## :floppy_disk: Data
The data folder contains the following input files:
* Data from Belgian study for fitting the evolutionary clock meta-analysis model
* Phylogenetic trees of major subtypes and circulating recombinant forms among Amsterdam MSM, with taxa labels for Amsterdam sequences comprising risk group (MSM/non-MSM), estimated age group at infection, estimated infection date and 95% credible intervals, world region of birth. U corresponds to unknown.
* Anonymised possible transmission pairs formulated from phylogenies and participant clinical data and metadata including covariates:
    * time elapsed
    * patristic distance
    * age of source and recipient at estimated time of transmission of the incident case
    * world region of birth of recipient
    * world region of birth of putative source
    * estimated year of possible transmission event

## :computer: Code
The scripts folder contains the following scripts to run the analysis:
1) run-stan.R - run stan model using cmdstan on anonymised data from formulated pairs
2) post-processing.R - write script to produce convergence diagnostics and generate figures and tables

**Note** the results from these scripts do not incorporate uncertainty in phylogenies and time since infection estimates, as in the manuscript which requires more detailed participant data, and uses pairs formed from central alignments and posterior median time since infection estimates only.

## :page_facing_up: License
The code in this repository is licensed under [CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg).

## :shield: Warranty
Imperial makes no representation or warranty about the accuracy or completeness of the data nor that the results will not constitute in infringement of third-party rights. Imperial accepts no liability or responsibility for any use which may be made of any results, for the results, nor for any reliance which may be placed on any such work or results.

## :books: Citation
Please cite this work as:

Blenkinsop, A., Pantazis, N., Kostaki, E.G., Sofocleous, L., van Sighem, A., Bezemer, D., van de Laar, T., van der Valk, M., Reiss, P., de Bree, G. and Ratmann, O., 2024. Sources of HIV infections among MSM with a migration background: a viral phylogenetic case study in Amsterdam, the Netherlands. *arXiv preprint arXiv:2401.08308*.


## Acknowledgements
This study received funding as part of the H-TEAM initiative from Aidsfonds (project number P29701), the Engineering and Physical Sciences Research Council (EP/X038440/1) and the Imperial College Department of Mathematics Cecilia Tanner Research Funding Scheme. The H-TEAM initiative is being supported by Aidsfonds (grant number: 2013169, P29701, P60803), Stichting Amsterdam Dinner Foundation, Bristol-Myers Squibb International Corp. (study number: AI424-541), Gilead Sciences Europe Ltd (grant number: PA-HIV-PREP-16-0024), Gilead Sciences (protocol numbers: CO-NL-276-4222, CO-US-276-1712, CO-NL-985-6195), and M.A.C AIDS Fund. The ATHENA database is maintained by Stichting HIV Monitoring and supported by a grant from the Dutch Ministry of Health, Welfare and Sport through the Centre for Infectious Disease Control of the National Institute for Public Health and the Environment. We thank all contributors to ATHENA and SHM, program staff and ATHENA participants; the Imperial College Research Computing Service (https://doi.org/10.14469/hpc/2232); and Zulip for sponsoring team communications through the Zulip Cloud Standard chat app.
