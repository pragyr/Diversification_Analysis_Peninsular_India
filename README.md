# Diversification in the Peninsular Indian Plate
## Authors
- [Pragyadeep Roy](https://sites.google.com/view/jahnavijoshi/team/phd-students#h.7rw6e76q98yf) (Orcid - [0009-0002-2215-2171](https://orcid.org/0009-0002-2215-2171), pragyadeep@csirccmb.org)
- [Jahnavi Joshi](https://sites.google.com/view/jahnavijoshi/team/jahnavi-joshi#h.6x6n3sbsag2p) (Orcid - [0000-0002-6015-4138](https://orcid.org/0000-0002-6015-4138), jahnavi@csirccmb.org)

## Versions used
R version used: 4.4.1\
Julia version used: 1.11.4

## Relevant Manuscript:
Title: Idiosyncrasies unveiled: examining the pace, patterns and predictors of biotic diversification in peninsular India\
DOI: [10.1111/ele.70160](https://onlinelibrary.wiley.com/doi/epdf/10.1111/ele.70160)
## Summary of the project
Hello,


Welcome to the this project. This project deals with understanding the pace, pattern and predictors of in-situ diversification in the Peninsular Indian Plate (PIP). Here we provide the list of packages, the data set and the computer scripts used in the study.

The Peninsular Indian Plate, one of the oldest regions of diversification in tropical Asia, harbours highly diverse and endemic biota. However, our understanding of the diversification dynamics of its biota within a quantitative framework remains limited. To address this, we used time-calibrated molecular phylogenies and birth-death models to examine the tempo, mode, and drivers of diversification across 33 well-studied endemic clades (~770 species). Among peninsular Indian clades, angiosperms diversified the fastest, while invertebrates diversified the slowest. Younger clades of Asian origin diversified more rapidly than the older, relictual Gondwanan clades. Evolutionary relatedness explained the disparities in diversification rates across taxonomic groups and biogeographic origins. A gradual accumulation of diversity was supported in 17 clades, suggesting that the historical stability of their habitat was an important driver. Miocene intensification of monsoons and aridification and fluctuations in paleotemperature explained diversification patterns in the remaining 16 clades. Our results highlight the role of regional biogeography, geoclimatic processes, and phylogenetic history in governing diversification dynamics in the tropics.

## Packages used in the study
Following is the list of packages used in the study. Let's load them.
```{r}
library(ape)
library(permute)
library(lattice)
library(vegan)
library(nlme)
library(picante)
library(npreg)
library(RPANDA)
library(coda)
library(deSolve)
library(TESS)
library(tibble)
library(pracma)
library(dplyr)
library(ggplot2)
library(CRABS)
library(iterators)
library(parallel)
library(foreach)
library(doParallel)
library(cluster)
library(maps)
library(phytools)
library(DDD)
library(dispRity)
library(phangorn)
library(Rcpp)
library(castor)
library(see)
library(ggpubr)
library(rr2)
library(MASS)
library(mvtnorm)
library(caper)
library(geiger)
library(grid)
library(Cairo)
library(qpdf)
library(deeptime)
library(ggdist)
library(gridExtra)
library(tidyr)
library(kmed)
library(FactoMineR)
library(factoextra)
library(mgcv)
library(ecole) # to install, execute this command: remotes::install_github("phytomosaic/ecole"). If you don't have the remotes package then install it using: install.packages("remotes")
library(ggstatsplot)
library(ggsci)
library(ggtree)
library(RRPP)
library(survival)
library(multcomp)
library(FSA)
library(ggbreak)
```
To install all the packages sequentially, replace library( with install.packages(\" and ) with \") from the above mentioned list.

## Data used in the study
The data that we have used in the project includes dated tree files and paleoclimate data. All the data would be present in the "data" directory.
Within the "data" directory, there are three folders and five files.
```{r}
list_dirs_first <- function(path = ".", ...) {
  files <- list.files(path, full.names = TRUE, ...)
  info <- file.info(files)
  dirs <- files[info$isdir]
  non_dirs <- files[!info$isdir]
  c(dirs, non_dirs)
}

# List folders first and then the files
list_dirs_first("data")
```
1. The **"all_trees"** folder contains the dated phylogenetic trees of the clades of our interest. These files were used the analyses that involved using the homogeneous birth-death models of package *RPANDA* and to run the CoMET (CPP on Mass Extinction Time) analysis. Few of the clades used in our analysis have smaller clade size and it has been seen that diversification estimates show more uncertainty with smaller clades. Hence, in addition to using the homogenous birth-death model from *RPANDA* and the CoMET analysis, we analysed our data using a hetergenous birth-death model - ClaDS (Cladogenetic Diverisification rate Shift) model. And, this model was implemented on larger trees (present in the **"trees_for_ClaDS_analyses"** folder) which contained those small clades of our interest. These dated trees were collected from multiple studies which have been cited in the in the supporting information (Section 1.4 and Table S1) of our study.

2. The **"paleoclimate_data"** folder contains the data for two environmental variables - Himalayan elevation and the expansion of C4-plants through time. This folder doesn't however include the temperature data. It is available in the *RPANDA* package with the name - *InfTemp*.
    - Data and column descriptions:
      - **InfTemp** - Data available with the RPANDA package. Columns:
        - *Age*: Age in million years ago (Mya)
        - *Temperature*: Temperature in $^\circ$C
      - **C_dat2_Potwar.csv** - Pedogenic Carbonate content through time. The data is relevant to the Potwar Plateau in Pakistan but serves as a proxy for intensified aridification across the Indian subcontinent. Columns:
        - *Age*: Age in million years ago (Mya)
        - *Delta.13_C*: Values of the carbon isotopic signature - $\delta$^13^C indicating the prevalent vegatation type. Higher values indicate high prevalence of C4 plants. Lower values indicate high prevalence of C3 plants.
      - **Him_Oro_data.csv** - Reconstructed Elevations of the Himalayan-Orogen. Columns:
        - *Age*: Age in million years ago (Mya)
        - *Elevation*: Average elevation in metres (m)
```{r}
# read/load the data
  data(InfTemp)
  C_exp <- read.csv("data/paleoclimate_data/C_dat2_Potwar.csv")
  Him_oro <- read.csv("data/paleoclimate_data/Him_Oro_data.csv")
  
  # print the column names and the few starting rows of each matrix
  head(InfTemp)
  head(C_exp)
  head(Him_oro)
```
3. The **"trees_for_ClaDS_analyses"** folder contains either the same phylogenies or the larger phylogenies to which a focal clade belongs. For clades with fairly large clade size we have used the same phylogenetic tree as is present in the "all_trees" folder. The larger phylogenies are majorly for the depauperate clades used in our study and these trees were specifically used in the ClaDS analysis.

4. The **"clads_tip_rates_with_attributes.csv"** file contains the speciation, extinction and net-diversification rates for each of the species (tip-rates) belonging to the 33 clades used in the study, estimated using the ClaDS model. In addition to the rates, the file also contains columns for the taxonomic, biogeographic and the habitat affinities of each species. Information from this file is used in generating Figure 2 (**Code 7**, mentioned below) and also to perform phylogenetic ANOVA (**Code 8**, mentioned below).
    - Columns and descriptions:
      - **Lineage**: Clade/Lineage a species/tip belongs to
      - **Speciation_Rate**: Speciation rate (estimated using ClaDS) relevant to a species (tip rate)
      - **Extinction_Rate**: Extinction rate (estimated using ClaDS) relevant to a species (tip rate)
      - **Net.diversification_Rate**: Net-diversification rate (estimated using ClaDS) relevant to a species (tip rate)
      - **Taxonomic.group"**: Taxonomic affinity of a species (Angiosperm/Invertebrate/Herpetofauna)
      - **Biogeographic.origin**: Biogeographic origin of a clade that a species belongs to (Asian/Gondwanan)
      - **Habitat**: Habitat preference of a species (Wet forests - "Wet"/Wet+Dry forests - "Wet+Dry"/Freshwater - "FW")
      - **tip_names**: Species name
```{r}
colnames(read.csv("data/clads_tip_rates_with_attributes.csv"))
```

5. The **"Final_summary_updated.csv"** file contains the summary of the identity of each of the lineages and the corresponding inference about its diversification. The cleaner version of this table has been presented with the relevant metadata as Table S2 in the Supporting Information of this study.
    - Columns and descriptions:
      - **Taxonomic.group"**: Taxonomic affinity of a clade/lineage (Angiosperm/Invertebrate/Herpetofauna)
      - **Lineage**: Clade/Lineage name
      - **Stem.age..Mya.**: The stem age of a clade, i.e., the age when a clade diverged from its sister lineage
      - **Crown.age..Mya.**: The crown age of a clade, i.e., the age of the first divergence in a lineage. This age is used as the clade age in the PGLS regression
      - **Stem.Clade.age..Mya.**: The phylogenetic pulse, i.e. the difference between the stem and crown age of a clade
      - **No..of.species**: The total species diversity of a clade
      - **Endemic.species**: The number of species among the total number of species within a clade that are endemic to the peninsular Indian plate
      - **Speciation.rate..RPANDA.**: The speciation rate of a clade estimated by fitting the constant rate birth-death model in the *RPANDA* framework
      - **Extinction.rate..RPANDA.**: The extinction rate of a clade estimated by fitting the constant rate birth-death model in the *RPANDA* framework
      - **Net.diversification.rate..RPANDA.**: The net-diversification rate of a clade estimated by fitting the constant rate birth-death model in the *RPANDA* framework
      - **Speciation.rate..CoMET.**: The mean speciation rate of a clade across time estimated by performing the CoMET analysis
      - **Extinction.rate..CoMET.**: The mean extinction rate of a clade across time estimated by performing the CoMET analysis
      - **Net.diversification.rate..CoMET.**: The mean net-diversification rate of a clade across time estimated by performing the CoMET analysis
      - **Speciation.rate..ClaDS.**: The mean speciation rate of a clade across branches and tips estimated by fitting the ClaDS model
      - **Extinction.rate..ClaDS.**: The mean extinction rate of a clade across branches and tips estimated by fitting the ClaDS model
      - **Net.diversification.rate..ClaDS.**: The mean net-diversification rate of a clade across branches and tips estimated by fitting the ClaDS model
      - **Pulled.Diversification.Rate**: The pulled diversification rate estimated by fitting the time-independent birth-death model with pulled diversification rates
      - **Biogeographic.origin**: Biogeographic origin of the clade (Asian/Gondwanan)
      - **Habitat**: Habitat preference of the clade (Wet forests - "Wet"/Wet+Dry forests - "Wet+Dry"/Freshwater - "FW")
      - **Speciation.rate.shift.times..TESS.**: The time(s) of strong episodic diversification rate shift(s) for a clade (if inferred) estimated by the CoMET analysis
      - **Results.without.CoMET**: The diversification scenario for a clade inferred by fitting the birth-death models of the *RPANDA* package
      - **Final.Results**: The diversification scenario for a clade inferred by combining the results of both the CoMET as well as the *RPANDA* models
      - **Num_Rate_Shifts**: The total number of distinct episodic rate shifts inferred for a clade by the SES (Symmetrical Episodic rate-Shift) model in the RPANDA framework and the CoMET analysis
      - **Rate_Shift_Time.CoMET.**: The age of the first rate shift (if inferred) within a clade estimated by the CoMET analysis
      - **Peak_or_Dip_Time.SES.**: The age of a strong peak or dip in diversification rates (if inferred) for a clade estimated by the SES model in the RPANDA framework
      - **Rate_Shift_Time_2.CoMET.**: The age of the second rate shift (if inferred in addition to an earlier rate shift) within a clade estimated by the CoMET analysis
      - **Driver**: The combination of strong drivers of diversification for a clade among - temperature (TDD), diversity (DDD), himalayan orogeny (Him.oro)
      - **DDD**: Presence or absence of strong diversity-dependent diversification (DDD) for a clade
      - **TDD**: Presence or absence of strong paleotemperature-dependent diversification (TDD) for a clade
      - **Him.oro**: Presence or absence of strong Himalayan-orogeny-dependent diversification for a clade
      - **C4.Exp**: Presence or absence of strong C4-Plant-Expansion-dependent diversification for a clade
      - **Other**: If none of the earlier mentioned drivers have had any substantial impact on the diversification of a clade
      - **Results..PDR.**: Inferences based on the pulled diversification rates through time for a clade
      - **CRABS**: Inferences about temporal diversification trends of a clade based on the CRABS test
```{r}
colnames(read.csv("data/Final_summary_updated.csv"))
```
  
6. The **"Scenarios_and_criteria.csv"** file contains information about the criteria for choosing different diversification scenarios which have been tested in the study. This file is read by one of the scripts mentioned below *(3_rpanda_code_and_RTT_comparisons_PR_JJ.R)* while inferring the diversification scenarios. Values in the cells for the parameters are the conditions, the combinations of which define a diversification scenario. Graphical illustrations for each of these scenarios has been presented in Figure S1.
    - Columns and descriptions:
      - **Scenario**: Diversification scenario. SC1 - Gradual Accumulation, SC2 - Saturated Accumulation, SC3 - Waxing and waning and SC4 - Exponential accumulation
      - **Sub-scenario**: Sub-scenarios as has been numbered and illustrated in Fig S1
      - **Scenarios**: Diversification scenario combined with the sub-scenario
      - **alpha**: The degree of decay or gain in speciation rates, i.e., the slope in the time-varying equation for speciation rates
      - **beta**: The degree of decay or gain in extinction rates, i.e., the slope in the time-varying equation for extinction rates
      - **lambda0-mu0**: The difference between speciation and extinction rates, i.e., the net-diversification rate
      - **additional remarks**: Additional conditions apart from the conditions for the other parameters
```{r}
head(read.csv("data/Scenarios_and_criteria.csv"))
```

7. The **"SUPER_TREE.nwk"** file contains a dated phylogenetic tree connecting all the clades used in the study in the newick format. It was generated in [www.timetree.org](www.timetree.org) and also using the information in the clades and respective studies and was used in some statistical analyses (Phylogenetic ANOVA and Phylogenetic Generalised Least-squares Regression - **Code 8**).

8. The **"SUPER_TREE_with_all_tips_ultrametric.nwk"** file contains a dated phylogenetic tree that connects all the clades and their respective phylogenies. This tree was built by integrating information from the "SUPER_TREE.nwk" file with the phylogenies of individual clades. This tree was also used in the same statistical analyses.


### Codes
All the relevant codes are present in the "codes" folder. Here is the list of all the codes.
```{r}
list.files("./codes")
```
The number that each code starts with indicates its position in the sequence of execution. Although the code names are self-explanatory, here's a brief account of the actions that the codes would perform:

1. The **"1_spline_interpolation_PR_JJ.R"** script facilitates smoothing the timeseries data for the paleoclimatic variables used in the analyses - Paleotemperature, Himalayan elevation and Expansion of C4-Plants through time. These smoothed curves can then be used in **Code 3** for fitting environmental-variable-dependent birth-death models.

2. The **"2_COMET_analysis_and_CRABS_PR_JJ.R"** script performs the CoMET (CPP on Mass Extinction Times) analysis and then extracts information about diverisification rates (speciation and extinction) through time and times of strong episodic rate shifts. Additionally, it assesses congruence classes of the inferred birth and death rates through time using the CRABS (Congruent Rate Analyses in Birth–death Scenarios) analysis.

3. The **"3_rpanda_model_and_scenario_selection_PR_JJ.R"** script performs likelihood-based homogeneous birth-death model-fitting in the *RPANDA* framework. Following this, diversification tempo (mean rates), temporal patterns (scenarios), and their drivers (diversity-dependence, environmental variables, etc) are estimated and inferred.

4. The **"4_run_ClaDS_in_Julia_PR_JJ.jl"** script peforms the ClaDS analysis in *Julia* and saves the results in .Rdata format. The first line of **Code 5** directly runs this code from R. So, one does not really need to explicitly run this code on the terminal. Just make sure that Julia is installed on your system.

5. The **"5_run_ClaDS_from_R_import_results_and_summarise_PR_JJ.R"** script runs **Code 4** from R, loads the .Rdata result files relevant to the ClaDS analyses and extracts mean rates and rates through time for the clades of interest belonging to the trees present in the **"trees_for_ClaDS_analyses"** folder.

6. The **"6_pulled_diversification_rates_finding_optimal_age_gridding.R"** script fits the homogeneous Pulled Diversification Rate (PDR) model on the data to assess the congruence classes in homogeneous birth-death models for each of the trees present in the **"all_trees"** folder. It tests if the age of a clade can be split into multiple time grids/bins (we checked for 1 to 5) so that diversification rates can vary across and within those time bins or remain constant throughout.

7. The **"7_make_plots_PR_JJ.R"** script is used to make plots submitted with the manuscript and the supporting information of this study.

8. The **"8_phyloANOVA_and_phyloRegression.R"** script performs phylogenetic ANOVA and regression using data from the  **"Final_summary_updated.csv"** and **"clads_tip_rates_with_attributes.csv"** files utilising the phylogenetic information from the **"SUPER_TREE.nwk"** and **"SUPER_TREE_with_all_tips_ultrametric.nwk"** tree files.

9. The **"9_simulate_trees_and_perform_model_fitting.R"** script simulates trees using the best model parameters from the RPANDA estimates (for temporal patterns) and performs posterior predictive tests to assess the credibility of the parameter estimates. The trees are simulated using the *tess_sim_age()* function from the *TESS* package and then the same analyses, as in **Codes 2 and 3** are performed. The estimates relevant to the simulated trees are then compared with the ones corresponding to the relevant empirical trees.

## Graphs
![Alt text](https://github.com/pragyr/Diversification_Analysis_Peninsular_India/blob/main/1_Conceptual%20figure2-1-1.png)
*A conceptual figure illustrating the predictions for the tempo and mode of diversification of endemic PI biota and the data used in the study*
