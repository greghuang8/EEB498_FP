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
# Version: 1.6
# Author: Greg Huang
# Last update: March 21, 2018
#
# Versions: 
#         1.1  Quick write up of the code
#         1.2  Clean up; Apply code to actual data
#         1.3  Added H-P beta-diversity analysis 
#         1.4  Plot edits
#         1.5  Add code for H-H beta diversity analysis,
#              also added code for SCBD significance testing
#         1.6  Add code for LCBD plotting. Script nearing completion.
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

# For ddply with plyr
library(plyr)


#### Load in .csv files ####
order_counts <- read.csv("zones_orders_counts_key.csv")
colnames(order_counts)[1] <- "Bats"

infection_counts <- read.csv("zones_infections_counts_key.csv")
colnames(infection_counts)[1] <- "Bacteria"

#prep for Host-parasite interaction beta diversity
#load in individual host-parasite interactions in each zone
HP_zone1 <- read.csv("zone1_hp_interactions.csv")
HP_zone2 <- read.csv("zone2_hp_interactions.csv")
HP_zone3 <- read.csv("zone3_hp_interactions.csv")
HP_zone4 <- read.csv("zone4_hp_interactions.csv")

#combine all these interactions
HP_counts <- bind_rows(HP_zone1, HP_zone2, HP_zone3, HP_zone4)
HP_counts[is.na(HP_counts)] <- 0
write.csv(HP_counts, file = "HP_combined_counts.csv")

#prep for Host-Host interaction beta diversity
#HH_all contains all host-host pairs that share the same parasite
HH_all <- read.csv("HH_Master.csv", header = FALSE)

#HH_unique takes only the unique entries and aggregate the duplicates
HH_unique <- ddply(HH_all,"V1",numcolwise(sum))

#Transpose the rows and columns
HH_transpose <- data.frame(t(HH_unique[-1]))

#Give new column and row names
colnames(HH_transpose) <- HH_unique[,1]
rownames(HH_transpose) <- c(1:4)

# For mapping: twc.xy contains xy coordinates of each point
twc.xy <- read.csv("twc_cellid_xy.csv")

# The LCBD values are mapped to the 699 unique locations
mapped.lcbds <- read.csv("LCBDs_699.csv")

#### Beta diversity calculation and Output ####

# SCBD - Species Contribution to Beta Diversity
# LCBD - Local Contribution to Beta Diversity

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

# Calculate beta-diversities for Host-Host interactions 
HH_transpose_hel <- decostand(HH_transpose, "hel")
(beta_r1_HH <- beta.div(HH_transpose))
beta_r1_HH <- (beta_r1_HH = beta.div(HH_transpose))

# Save results into data frames and txt outputs
# order
order_SCBD <- data.frame(beta_r1_order$SCBD)
order_LCBD <- data.frame(beta_r1_order$LCBD)
write.csv(beta_r1_order$SCBD, file="beta_r1_order_SCBD.csv")
write.csv(beta_r1_order$LCBD, file="beta_r1_order_LCBD.csv")

# infections
infections_SCBD <- data.frame(beta_r1_infections$SCBD)
infections_LCBD <- data.frame(beta_r1_infections$LCBD)
write.csv(beta_r1_infections$SCBD, file="beta_r1_infections_SCBD.csv")
write.csv(beta_r1_infections$LCBD, file="beta_r1_infections_LCBD.csv")

# H-P interactions (order - infections)
HP_SCBD <- data.frame(beta_r1_HP$SCBD)
HP_LCBD <- data.frame(beta_r1_HP$LCBD)
write.csv(beta_r1_HP$SCBD, file="beta_r1_HP_SCBD.csv")
write.csv(beta_r1_HP$LCBD, file="beta_r1_HP_LCBD.csv")

# H-H interactions (order - order)
HH_SCBD <- data.frame(beta_r1_HH$SCBD)
HH_LCBD <- data.frame(beta_r1_HH$LCBD)
write.csv(beta_r1_HH$SCBD, file="beta_r1_HH_SCBD.csv")
write.csv(beta_r1_HH$LCBD, file="beta_r1_HH_LCBD.csv")

#### Check for SCBD significance ####

# Hosts
order_SCBD_mean <- mean(order_SCBD$beta_r1_order.SCBD)
order_sig_SCBD <- which(order_SCBD$beta_r1_order.SCBD >= order_SCBD_mean)
# int[1:3] 4 7 9
# three significant SCBD values for our hosts

# Parasites
infection_SCBD_mean <- mean(infections_SCBD$beta_r1_infections.SCBD)
infection_sig_SCBD <- which(infections_SCBD$beta_r1_infections.SCBD >= 
                              infection_SCBD_mean)
# 2 significant SCBD values for parasites

# Host-Parasite Interactions
HP_SCBD_mean <- mean(HP_SCBD$beta_r1_HP.SCBD)
HP_sig_SCBD <- which(HP_SCBD$beta_r1_HP.SCBD >= HP_SCBD_mean)

# Host-Host Interactions
HH_SCBD_mean <- mean(HH_SCBD$beta_r1_HH.SCBD)
HH_sig_SCBD <- which(HH_SCBD$beta_r1_HH.SCBD >= HH_SCBD_mean)

#### Check for LCBD significance ####
# count the locations of p.LCBD values that has significance (p) <= 0.05
(order_r1_sig_LCBD <- which(beta_r1_order$p.LCBD <= 0.05))
(infection_r1_sig_LCBD <- which(beta_r1_infections$p.LCBD <= 0.05))
(HP_r1_HP_LCBD <- which(beta_r1_HP$p.LCBD <= 0.05))
(HP_r1_HH_LCBD <- which(beta_r1_HH$p.LCBD <= 0.05))
# not too good, but we'll plot those

#### Maps of LCBD values and richness per quadrat ####

plot(twc.xy, asp=1, type="n",
     xlab="Longitude (m)", ylab="Latitude (m)",
     main="LCBD indices, for hosts",
     xlim=c(-81.5,-78), ylim=c(42.5,44.5))
points(twc.xy,pch=15, col="red", bg="red",
       cex=1*sqrt(mapped.lcbds$orderLCBD))

plot(twc.xy, asp=1, type="n",
     xlab="Longitude (m)", ylab="Latitude (m)",
     main="LCBD indices, for infections",
     xlim=c(-81.5,-78), ylim=c(42.5,44.5))
points(twc.xy,pch=15, col="steelblue2", bg="steelblue2",
       cex=1*sqrt(mapped.lcbds$InfectionLCBD))

plot(twc.xy, asp=1, type="n",
     xlab="Longitude (m)", ylab="Latitude (m)",
     main="LCBD indices, for H-P interactions",
     xlim=c(-81.5,-78), ylim=c(42.5,44.5))
points(twc.xy,pch=15, col="green", bg="steelblue2",
       cex=1*sqrt(mapped.lcbds$H.PLCBD))

plot(twc.xy, asp=1, type="n",
     xlab="Longitude (m)", ylab="Latitude (m)",
     main="LCBD indices, for H-H interactions",
     xlim=c(-81.5,-78), ylim=c(42.5,44.5))
points(twc.xy,pch=15, col="orange", bg="steelblue2",
       cex=1*sqrt(mapped.lcbds$HHLCBD))

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

#plot HH
plot(sort(beta_r1_HH$SCBD,decreasing=T),type="n", 
     main = "SCBD of H-H interaction", ylab = "SCBD")
text(sort(beta_r1_HH$SCBD,decreasing=T),
     labels=names(sort(beta_r1_HH$SCBD,decreasing=T)),cex=0.5)

#### END ####