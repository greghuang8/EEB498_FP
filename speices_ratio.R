# species_ratio.R
#
# Script for species percentages
# Large data set: TWC_Species_LU_Feb2015, 26100 entries 
# Against small data set: all_entries_updated_zones, 3060 entries
#
# Version: 1.0
# Author: Greg Huang
# Last update: May 15th, 2018
#
# Version history:
#                 1.0 Initial script creation 
#
# ==============================================================================

#### Required packages ####
install.packages("stringr")
library(stringr)

#### Load in data sets ####
big <- read.csv("TWC_Species_LU_Feb2015.csv", header = T, 
                stringsAsFactors = FALSE)
small <- read.csv("all_entries_updated_zones.csv", header = T, 
                  stringsAsFactors = FALSE)


#### Determine which ratios to look for ####
unique(small$ClassOrder)
# [1] "Small Carnivores"        "Pigeons and Doves"       "Large Carnivores"       
# [4] "Small Rodents and Moles" "Waterfowl and Waders"    "Bats"                   
# [7] "Songbirds"               "Seabirds and Shorebirds" "Large Rodents"          
# [10] "Lagomorphs"              "Falcons and Hawks"       "Owls" 

unique(big$Group_)
#[1] "Rodents and Moles"        "Songbirds"               "Small Carnivores and Opossums"
#[4] "Seabirds and Shorebirds"  "Large-bodied Mammals"    "Waterfowl and Waders"         
#[7] "Reptiles"                 "Raptor"                  "Winged Mammals"               
#[10] "Exotics"                 NA                        "Amphibians"                   
#[13] "Fish"                    "Gamebirds"    

# replace NA with string "NA"
big[is.na(big)] <- "NA"

#### Group-wise ratio calculations ####

# Small carnivores
a <- sum(str_count(small$ClassOrder, "Small Carnivores")) #595
b <- sum(str_count(big$Group_, "Small Carnivores and Opossums")) #2589
small_carnivore_ratio <- a/b  #0.23
  
# Small rodents and moles
a <- sum(str_count(small$ClassOrder, "Small Rodents and Moles")) #844
b <- sum(str_count(big$Group_, "Rodents and Moles")) #4403
small_rodent_ratio <- a/b     #0.19

# Lagomorphs
a <- sum(str_count(small$ClassOrder, "Lagomorphs")) #5
b <- sum(str_count(big$Group_, "Rodents and Moles")) #4403
lagomorph_ratio <- a/b        #0.001

# Songbirds 
a <- sum(str_count(small$ClassOrder, "Songbirds")) #83
b <- sum(str_count(big$Group_, "Songbirds")) #11999
songbirds_ratio <- a/b       #0.007

# Pigeons 
a <- sum(str_count(small$ClassOrder, "Pigeons and Doves")) #1236 (1/3)
b <- sum(str_count(big$Group_, "Songbirds")) #11999
pigeons_ratio <- a/b         #0.10

# Waterfowl and Waders
a <- sum(str_count(small$ClassOrder, "Waterfowl and Waders")) #54
b <- sum(str_count(big$Group_, "Waterfowl and Waders")) #4249
waterfowl_ratio <- a/b       #0.013

# Seabirds and Shorebirds
a <- sum(str_count(small$ClassOrder, "Seabirds and Shorebirds")) #31
b <- sum(str_count(big$Group_, "Seabirds and Shorebirds")) #1126
seabirds_ratio <- a/b        #0.028

# Falcons and Hawks
a <- sum(str_count(small$ClassOrder, "Falcons and Hawks"))  #29
b <- sum(str_count(big$Group_, "Raptor")) #850
falcons_ratio <- a/b          #0.034

#Large Carnivores !!!!!!!!!
a <- sum(str_count(small$ClassOrder, "Large Carnivores")) #118
b <- sum(str_count(big$Group_, "Large-bodied Mammals")) #99


# Bats
a <- sum(str_count(small$ClassOrder, "Bats")) #55
b <- sum(str_count(big$Group_, "Winged Mammals")) #349
bats_ratio <- a/b             #0.157

# Large Rodents
a <- sum(str_count(small$ClassOrder, "Large Rodents")) #5
b <- sum(str_count(big$Group_, "Rodents and Moles")) #4403
large_rodents_ratio <- a/b    #0.001

# Owls
a <- sum(str_count(small$ClassOrder, "Owls")) #5
b <- sum(str_count(big$Group_, "Raptor")) #850
owls_ratio <- a/b             #0.006


#### Species-wise ratio calculations ####

#uniq_spec <- as.data.frame(unique(big$Species))
#simple_spec <- gsub(",.*", "", uniq_spec$`unique(big$Species)`)
#uniq_simple_spec <- as.data.frame(unique(simple_spec))

# Step 1: Merge the two data sets 
big_w_spec <- read.csv("Big_ID_spec_only.csv", stringsAsFactors = FALSE)
small_IDs <- read.csv("Small_ID_only.csv", stringsAsFactors = FALSE)

matched <- big_w_spec$BIG_ID %in% small_IDs$Small_ID
table(matched)["TRUE"]
new <- as.data.frame(big_w_spec[matched,])
write.csv(new, "small_set_with_species.csv")

list_of_species <- unique(new$Species)
#Identification/ratio of the 45 different species

match2 <- small_IDs$Small_ID %in% big_w_spec$BIG_ID
table(match2)["TRUE"]
new2 <- as.data.frame(big_w_spec[match2,])

#Test
a <- sum(str_count(new$Species, "Coyote")) #37
b <- sum(str_count(big_w_spec$Species, "Coyote")) #74

export_df <- data.frame("Test1", 5, 10, 0.5)
names(export_df) <- c("Species", "Count in subset", "Count in original", "ratio")

# For each of the unique items in list_of_species, calculate the ratios
for (i in list_of_species){
  a <- sum(str_count(new$Species, i))
  b <- sum(str_count(big_w_spec$Species, i))
  ratio <- a/b
  df <- data.frame(i, a, b, ratio)
  names(df) <- c("Species", "Count in subset", "Count in original", "ratio")
  export_df <- rbind(export_df,df) 
}

write.csv(export_df, "species_ratios.csv")
unique(big$Species)
unique(big_w_spec$Species)
test <- as.data.frame(table(big_w_spec$Species))
match3 <- !(big_w_spec$Species %in% export_df$Species)
new3 <- as.data.frame(big_w_spec[match3,])
not_included <- as.data.frame(table(new3$Species))

write.csv(not_included, "unselected_species_counts.csv")
barplot()