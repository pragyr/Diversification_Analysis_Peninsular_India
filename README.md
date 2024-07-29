# Diversification in the Peninsular Indian Plate
## Authors
- Pragyadeep Roy
- Jahnavi Joshi

Hello,
Welcome to the this project. It deals with understanding the pace, pattern and predictors of in-situ diversification in the Peninsular Indian Plate (PIP). The data that we have used in the project includes dated tree files and paleoclimate data.

## Data used in the study
All the data would be present in the "data" directory.
Within the "data" directory, there are three folders.
```{r}
list.dirs("./data")
```

The "all_trees" folder as the name suggests contains all the dated phylogenetic trees of the clades of our interest. These files were used the analyses that involved using the homogeneous birth-death models of package RPANDA and to run the CoMET (CPP on Mass Extinction Time) analysis. Few of the clades used in our analysis have smaller clade size and it has been seen that diversification estimates show more uncertainty with smaller clades. Hence, in addition to using the homogenous birth-death model from RPANDA and the CoMET analysis, we analysed our data using a hetergenous birth-death model - ClaDS (Cladogenetic Diverisification rate Shift) model. And, this model was implemented on larger trees which contained those small clades of our interest and those "super trees" were kept in the "all_trees_or_super_trees" folder. For clades with fairly large clade size we have used the same phylogenetic tree as is present in the "all_trees" folder.

Paleoclimate data includes global paleotemperature (for the cenozoic era), Himalayan elevation, Expansion of C4-plants through time. The "paleoclimate_data" folder doesn't however include the temperature data. It is available in the RPANDA package
```{r}
library(RPANDA)
data(InfTemp)
head(InfTemp)
```

Outside these folders the file - "Final_summary_updated.csv" summary information about the identity of each of the lineages and the corresponding inference about its diversification. More details about the analyses can be found in the supplementary information of the study. The "SUPER_TREE.tre" file is a dated phylogenetic tree connecting all the clades used in the study. It was generated in [www.timetree.org](www.timetree.org) and also using the information in the clades and respective studies and was used in one of the analysis (PGLS) in the study.

### Codes
All the relevant codes are present in the "codes" folder. Here is the list of all the codes.
```{r}
list.files("./codes")
```
Although the code names are self-explanatory, here's a brief account of the actions that the codes would perform:

"1_spline_interpolation_PR_JJ.R" - for smoothing the paleoclimatic variables used in the analyses

"2_COMET_analysis_and_CRABS_PR_JJ.R" - for performing the CoMET analysis, to extract information about diverisification rates (speciation and extinction) through time and times of strong episodic rate shifts, and to assess congruence classes using the CRABS analysis

"3_rpanda_code_and_RTT_comparisons_PR_JJ.R" - for estimating diversification rates, patterns and drivers using homogenous birth-death models from the package RPANDA and to compare the estimated rates through time with that of CoMET analysis

"4_run_ClaDS_in_Julia_PR_JJ.jl" - to peform the ClaDS analysis in Julia and to save the results in Rdata format

"5_export_ClaDS_data_and_summarise_in_R_PR_JJ.R" - to load the Rdata result files relevant to the ClaDS analysis and to extract rates through time for clades of interest

"6_pulled_diversification_rates_finding_optimal_age_gridding.R" - to fit the homogenous pulled diversification rate model on the data to assess general trends in pulled rates

"7_make_plots_PR_JJ.R" - to make plots submitted with this manuscript
