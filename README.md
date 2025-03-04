# Diversification in the Peninsular Indian Plate
## Authors
- [Pragyadeep Roy](https://sites.google.com/view/jahnavijoshi/team/phd-students#h.7rw6e76q98yf) (Orcid - [0009-0002-2215-2171](https://orcid.org/0009-0002-2215-2171), pragyadeep@csirccmb.org)
- [Jahnavi Joshi](https://sites.google.com/view/jahnavijoshi/team/jahnavi-joshi#h.6x6n3sbsag2p) (Orcid - [0000-0002-6015-4138](https://orcid.org/0000-0002-6015-4138), jahnavi@csirccmb.org)

## Relevant Manuscript Title:
Idiosyncrasies unveiled: examining the pace, patterns and predictors of biotic diversification in peninsular India

## Summary
Hello,

Welcome to the this project. It deals with understanding the pace, pattern and predictors of in-situ diversification in the Peninsular Indian Plate (PIP). The data that we have used in the project includes dated tree files and paleoclimate data. The manuscript relevant to this work is up on [BiorXiv](https://www.biorxiv.org/content/10.1101/2023.11.15.567174v3) and is under review. 

## Data used in the study
All the data would be present in the "data" directory.

### Phylogenetic trees
The "all_trees" folder as the name suggests contains all the dated phylogenetic trees of the clades of our interest. These files were used in the analyses that involved using the homogeneous birth-death models of package RPANDA and to run the CoMET (CPP on Mass Extinction Time) analysis. Few of the clades used in our analysis have smaller clade size and it has been seen that diversification estimates show more uncertainty with smaller clades. Hence, in addition to using the homogenous birth-death model from RPANDA and the CoMET analysis, we analysed our data using a hetergenous birth-death model - ClaDS (Cladogenetic Diverisification rate Shift) model. And, this model was implemented on larger trees which contained those small clades of our interest and those "super trees" were kept in the "all_trees_or_super_trees" folder. For clades with fairly large clade size we have used the same phylogenetic tree as is present in the "all_trees" folder.

### Paleoclimatic data
Paleoclimate data includes global paleotemperature (for the cenozoic era), Himalayan elevation, Expansion of C4-plants through time. The "paleoclimate_data" folder doesn't however include the temperature data. It is available in the RPANDA package
```{r}
library(RPANDA)
data(InfTemp)
head(InfTemp)
```

Outside these folders the file - "Final_summary_updated.csv" summary information about the identity of each of the lineages and the corresponding inference about its diversification. More details about the analyses can be found in the supplementary information of the study. The "SUPER_TREE.tre" file is a dated phylogenetic tree connecting all the clades used in the study. It was generated in [www.timetree.org](www.timetree.org) and also using the information in the clades and respective studies and was used in one of the analysis (PGLS) in the study.

## Codes
All the relevant codes are present in the "codes" folder. Here is the list of all the codes.

Although the code names are self-explanatory, here's a brief account of the actions that the codes would perform:

"1_spline_interpolation_PR_JJ.R" - for smoothing the paleoclimatic variables used in the analyses

"2_COMET_analysis_and_CRABS_PR_JJ.R" - for performing the CoMET analysis, to extract information about diverisification rates (speciation and extinction) through time and times of strong episodic rate shifts, and to assess congruence classes using the CRABS analysis

"3_rpanda_code_and_RTT_comparisons_PR_JJ.R" - for estimating diversification rates, patterns and drivers using homogenous birth-death models from the package RPANDA and to compare the estimated rates through time with that of CoMET analysis

"4_run_ClaDS_in_Julia_PR_JJ.jl" - to peform the ClaDS analysis in Julia and to save the results in Rdata format

"5_export_ClaDS_data_and_summarise_in_R_PR_JJ.R" - to load the Rdata result files relevant to the ClaDS analysis and to extract rates through time for clades of interest

"6_pulled_diversification_rates_finding_optimal_age_gridding.R" - to fit the homogenous pulled diversification rate model on the data to assess general trends in pulled rates

"7_make_plots_PR_JJ.R" - to make plots submitted with this manuscript

"8_phyloANOVA_and_phyloRegression.R" - to perform standard and phylogenetic analysis of variance

"9_simulate_trees_and_perform_model_fitting.R" - to simulate trees using RPANDA model parameters, check the propensity of the recovery of the starting model parameters and calculate the credibility of empirical estimates

## Graphs
![Alt text](https://github.com/pragyr/Diversification_Analysis_Peninsular_India/blob/main/1_Conceptual%20figure2-1-1.png)
*A conceptual figure illustrating the predictions for the tempo and mode of diversification of endemic PI biota and the data used in the study*
