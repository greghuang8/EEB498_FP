# zone_reclassify.R
#
# Reclassifying the zones
# The zones will be based on standardized population density (pdSTD).
# pdSTD is logged for 700 locations. 
# Since the main file with all entries have duplicates (one location
# may have more than one H-P interaction logged), the purpose of this
# script is to map the 700 pdSTD values to the 3000+ entries, so the 
# new zones can be recalculated.
#
# Version:  1.2
# Author: Greg Huang
# Last update: March 16, 2018
#
# Versions: 
#         1.1  Quick write up of the code
#         1.2  Edit the comments
# =============================================================================
####Prep stage####
#load in the original file with 3060 entries
original <- read.csv("Fortin_TWC_LU_2016.csv")
#load in the file with the 699 pdSTD entries (all unique)
popDen <- read.csv("ID_PopDen.csv")

#now we have a essentially a directory - 
#for each unique ID there is an associated pdSTD
id_popDen <- data.frame(popDen$CELL_ID,popDen$pdSTD)
colnames(id_popDen) <- c("CELL_ID","pdSTD")

####Mapping to Original####
result <- merge(original,id_popDen)

####Export####
write.csv(result, "all_entries.csv")
#Re-write the 4 zones in excel for submission. Finished. 

####End####

