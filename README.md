# EEB498_FP

## Introduction
This project contains all necessary code and files to produce the personal project (excluding work-study gibberish) for EEB498.

## File names and contents/purpose

Here are the descriptions for each of the files in this project. 

### CSV files

**Fortin_TWC_LU_2016**

All relevant data stored in this file. This includes lat-long location, host and parasite information, date, and zoning. 

**ID_PopDen**

Shows the calculated population density (pdSTD) for each of the location IDs

**pd_included**

Somewhat similar to ID_PopDen, but one thing to note is that htis csv file also contains zoning information which is used for the main analysis.

**TWC_PopDen_GME_March2015**

Includes population density data for all locations

### R scripts

**pop_density**

A short script determining the bin sizes of the population density heat map. The bins become the guide to dividing up the data into zones.

**watts_bipart_script**

Script initially developed for Dr. Watt's PhD thesis; here I make couple changes to accomodate the updated zones.


### PDF and TXT files

**Project_Proposal_Greg_498**

Pretty self explanatory. 

**correlation_outputs**

Results for correlation analysis what compares pdSTD to each of the sub-categories for the initial zoning configuration.
