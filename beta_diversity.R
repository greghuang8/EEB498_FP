# beta_diversity.R
#
# Calculate the beta-diversity for each zones 
# After the zones have been reclassified into four separate sections by pop.
# density (prepped by zone_reclassify.R and excel), it is now time to run 
# beta-diversity calculations in each of those zones. The code will first
# run beta-diversity in the general cases (e.g. updated_zone1_all.csv), 
# then potentially for the more refined, yearly zone csv files.
# The beta diversities will encompass the hosts and parasites. 
# 
# Version: 1.4
# Author: Greg Huang
# Last update: March 20, 2018
#
# Versions: 
#         1.1  Quick write up of the code
#         1.2  Clean up; Apply code to actual data
#         1.3  Added H-P beta-diversity analysis 
#         1.4  Plot edits
# ==============================================================================

#### Install required packages ####
if(!require(vegan, quietly = TRUE)){
  install.packages("vegan")
  library(vegan)
}

# package adespatial contains beta.div,
# the tool I'll use to get the beta diversity calculations
if(!require(adespatial, quietly = TRUE)){
  install.packages("adespatial")
  library(adespatial)
}

# For rbind with dplyr
library(dplyr)


#### Load in .csv files ####
order_counts <- read.csv("zones_orders_counts_key.csv")
colnames(order_counts)[1] <- "Bats"

infection_counts <- read.csv("zones_infections_counts_key.csv")
colnames(infection_counts)[1] <- "Bacteria"

HP_zone1 <- read.csv("zone1_hp_interactions.csv")
HP_zone2 <- read.csv("zone2_hp_interactions.csv")
HP_zone3 <- read.csv("zone3_hp_interactions.csv")
HP_zone4 <- read.csv("zone4_hp_interactions.csv")
HP_counts <- bind_rows(HP_zone1, HP_zone2, HP_zone3, HP_zone4)
HP_counts[is.na(HP_counts)] <- 0

#### Beta diversity calculation ####

# LCBD - Local Contribution to Beta Diversity
# SCBD - Species Contribution to Beta Diversity

# Calculate beta-diversities for orders across zones
order_counts_hel <- decostand(order_counts, "hel") # Hellinger transform
(beta_r1_order <- beta.div(order_counts))  
beta_r1_order<-(beta_r1_order = beta.div(order_counts))

# Calculate beta-diversities for infections across zones
infection_counts_hel <- decostand(infection_counts, "hel")
(beta_r1_infections <- beta.div(infection_counts))  
beta_r1_infections<-(beta_r1_infections = beta.div(infection_counts))

# Calculate beta-diversities for Host-parasite interactions 
# across zones
HP_counts_hel <- decostand(HP_counts, "hel")
(beta_r1_HP <- beta.div(HP_counts))
beta_r1_HP <- (beta_r1_HP = beta.div(HP_counts))

# Save results into data frames
# order
order_SCBD <- data.frame(beta_r1_order$SCBD)
order_LCBD <- data.frame(beta_r1_order$LCBD)

# infections
infections_SCBD <- data.frame(beta_r1_infections$SCBD)
infections_LCBD <- data.frame(beta_r1_infections$LCBD)

# H-P interactions (order - infections)
HP_SCBD <- data.frame(beta_r1_HP$SCBD)
HP_LCBD <- data.frame(beta_r1_HP$LCBD)


#### Check for significance ####

# count the locations of p.LCBD values that has significance (p) <= 0.05
(order_r1_sig_LCBD <- which(beta_r1_order$p.LCBD <= 0.05))
# not good: it's all pretty high, return null
(infection_r1_sig_LCBD <- which(beta_r1_infections$p.LCBD <= 0.05))
(HP_r1_sig_LCBD <- which(beta_r1_HP$p.LCBD <= 0.05))


#### Maps of LCBD values and richness per quadrat ####
# Leave out for now: twc.xy unavailable 
# plot(twc.xy, asp=1, type="n", xlab="Longitude (m)", ylab="Latitude (m)", main="LCBD indices, microphytophagous mites", xlim=c(0,2.5), ylim=c(0,10))
# points(twc.xy,pch=21, col="white", bg="steelblue2", cex=18*sqrt(beta.res1$LCBD))
# points(twc.xy[res1.signif.LCBD,],pch=21, col="white", bg="red", cex=18*sqrt(beta.res1$LCBD[res1.signif.LCBD]))
# text(twc.xy, labels=1:70, pos=4, cex=0.8, offset=0.7)


#### INAKI - SCBD Plots ####

#plot order
plot(sort(beta_r1_order$SCBD,decreasing=T),type="n", 
     main = "SCBD of Host (order)", ylab = "SCBD")
text(sort(beta_r1_order$SCBD,decreasing=T),
     labels=names(sort(beta_r1_order$SCBD,decreasing=T)),cex=0.5)

#plot infections
plot(sort(beta_r1_infections$SCBD,decreasing=T),type="n",
     main = "SCBD of Parasite (infection)", ylab = "SCBD")
text(sort(beta_r1_infections$SCBD,decreasing=T),
     labels=names(sort(beta_r1_infections$SCBD,decreasing=T)),cex=0.5)

#plot HP
plot(sort(beta_r1_HP$SCBD,decreasing=T),type="n", 
     main = "SCBD of H-P interaction (H-P)", ylab = "SCBD")
text(sort(beta_r1_HP$SCBD,decreasing=T),
     labels=names(sort(beta_r1_HP$SCBD,decreasing=T)),cex=0.5)
