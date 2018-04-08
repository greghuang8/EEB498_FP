# EEB498_FP

This project contains all necessary code and files to produce the personal project (excluding work-study gibberish) for EEB498.
Code developed for this project is in beta_diversity.R, with some minor edits in watt_bipart_script.R. 

## File names and contents/purpose

Here are the descriptions for each of the files in this project. 

### CSV files

**Fortin_TWC_LU_2016**

All relevant data stored in this file. This includes lat-long location, host and parasite information, date, and zoning. 

**HH_Master**

Shows all Host-Host interactions and its frequencies. A host-to-host interaction is counted when two different hosts (different in order) in the same zone were recorded to be infected by the same kind of parasite. 

**HP_combined_counts**

Shows all Host-parasite interactions.

**ID_PopDen**

Shows the calculated population density (pdSTD) for each of the location IDs

**LCBDs_699**

Shows the LCBD values for the 699 unique locations found in the TWC data set after they have been calculated by beta_diversity.R. 

**pd_included**

Somewhat similar to ID_PopDen, but one thing to note is that htis csv file also contains zoning information which is used for the main analysis.

**twc_cellid_xy**

This .csv file shows the x-y coordinates of each of the Cell IDs. This csv file is used to plot out the beta-diversity results in R. 

**TWC_PopDen_GME_March2015**

Includes population density data for all locations.

### R scripts

**beta_diversity**

Main focus of this study. This R script calculates the beta-diversities (as well as species richness and simpson's and shannons diversities) of hosts, parasites, host-parasites, and hosts-hosts in the GTA, using the updated urbanization gradient. Plots of the Local contribution to beta diversity (LCBD) values are also created here. 

**pop_density**

A short script determining the bin sizes of the population density heat map. The bins become the guide to dividing up the data into zones.

**watts_bipart_script**

Script initially developed for Dr. Watt's PhD thesis; here I make couple changes to accomodate the updated zones.


### PDF and TXT files

**Project_Proposal_Greg_498**

Pretty self explanatory. 

**correlation_outputs**

Results for correlation analysis what compares pdSTD to each of the sub-categories for the initial zoning configuration.
