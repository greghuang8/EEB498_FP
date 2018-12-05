#### BIPARTITE RCODE FOR HOST-PARASITE MUTUALISTIC NETWORKS##############
#### A.WATTS 2015 #######################################################
#### Selected functions from bipartite package (C. Dormann et al. 2014)##
#### Edited by Greg Huang################################################
#########################################################################

##############LIBRARIES####################

library(vegan)
library(bipartite)

######READ IN URBAN DATA (MARCH2015)

twc = read.csv("updated_zone4_all.csv", header=TRUE)
str(twc)

###TRANSPOSE: Create an interaction matrix from a prexisting dataframe (community-matrix style, with sites as rows and observations by column)
##Choose the name of the columns that represent two interaction species or species groups (i.e. one column for hosts, one column for parasites, at the species level or the functional level)
#In this example, I use 'Species.x' column and the 'Infection' column in the 'twc.urban.region' datatable, representing host species and parasite species, respectively.
#call this table something new (I use bipart) to be used in the bipartite package.

xtabs(~ ClassOrder + InfectionType, data=twc)

bipart_twc = as.table(xtabs(~ ClassOrder + InfectionType, data=twc))
##I Ususally write this matrix so that I only have to do this step once. 
write.csv(bipart_twc, file="bipart_twc.csv")

bipart_twc

###I needed to 'clean' my data matrix so that it is suitable for bipartite analyses. 
bipart_twc_all = read.csv("bipart_twc.csv", header=T, check.names=FALSE)

bipart_twc_all

rownames(bipart_twc_all) = bipart_twc_all[,1]
bipart_twc_all = bipart_twc_all[,-1]

bipart_twc_all


####BIPARTITE PACKAGE:

###PLOTWEB: Visualize a bipartite interaction matrix (e.g. a foodweb)
plotweb_twc = plotweb(bipart_twc_all,
                      method = "normal", empty = TRUE, labsize = 0.9, ybig = 1.8, y.width.low = 0.1,
                      y.width.high = 0.1, low.spacing = NULL, high.spacing = NULL,
                      arrow="down", col.interaction="gray30", col.high = "ivory3",
                      col.low="grey39", bor.col.interaction ="gray30", bor.col.high="ivory3",
                      bor.col.low="grey39", high.lablength = NULL, low.lablength = NULL,
                      sequence=NULL, low.abun=NULL, low.abun.col="green",
                      bor.low.abun.col ="black", high.abun=NULL, high.abun.col="red",
                      bor.high.abun.col="black", text.rot=90, text.high.col="black",
                      text.low.col="black", adj.high=NULL, adj.low=NULL, plot.axes = FALSE,
                      low.y=1.5, high.y=2.4, add=FALSE, y.lim=NULL, x.lim=NULL, low.plot=TRUE,
                      high.plot=TRUE, high.xoff = 0, low.xoff = 0, high.lab.dis = NULL,
                      low.lab.dis = NULL, abuns.type="additional")

#visweb: Plotting function to visualize a bipartite food web

visweb_twc = visweb(bipart_twc_all, type="nested", prednames=TRUE, preynames=TRUE, labsize=1,
                    plotsize=25, square="interaction", text="no", frame=NULL, textsize=1,
                    textcol="red", pred.lablength=NULL, prey.lablength=NULL, clear=TRUE,
                    xlabel="", ylabel="", boxes=TRUE, circles=FALSE, circle.col="black",
                    circle.min=0.2, circle.max=2, outerbox.border="white",
                    outerbox.col="white", box.border="black", box.col="black", def.col="blue",
                    max.digits=4, NA.col="red")

######networklevel: Analysis of bipartite webs at the level of the entire network
networkmetrics_twc = networklevel(bipart_twc_all, index="ALLBUTDD", level="both", weighted=TRUE,
                                  ISAmethod="Bluethgen", SAmethod = "Bluethgen", extinctmethod = "r",
                                  nrep = 100, CCfun=median, dist="horn", normalise=TRUE, empty.web=TRUE,
                                  logbase="e", intereven="prod", H2_integer=TRUE, fcweighted=TRUE,
                                  fcdist="euclidean", legacy=FALSE)

##I usually export this as a *.csv file so that I can consolidate the metrics in Excel. HL refers to 'higher-level species' and LL refers to 'lower-level species'.
write.csv(networkmetrics_twc, file="networkmetrics_twc.csv")


##### specieslevel: Calculate various indices for network properties at the species level
##Choose 'lower' or 'higher' depending on the level you need 
speciesmetrics_twc = specieslevel(bipart_twc_all, index="ALLBUTD", level="lower", logbase=exp(1), low.abun=NULL,
                                  high.abun=NULL, PDI.normalise=TRUE, PSI.beta=c(1,0), nested.method="NODF",
                                  nested.normalised=TRUE, nested.weighted=TRUE, empty.web=TRUE)

write.csv(speciesmetrics_twc, file="speciesmetrics_twc.csv")

######SPECIALIZATION METRICS
#dfun Calculates standardised specialisation index d' (d prime) for each
##species in the lower trophic level of a bipartite network.
dfun_bipart_twc<-dfun(bipart_twc_all, abuns=NULL)
write.csv(dfun_bipart_twc, file="dfun_bipart_twc.csv")

###computeModule: This function takes a bipartite weighted graph and computes modules by applying M. E. J. Newman?????s
#modularity measure in a bipartite weighted version to it.
modules_bipart_twc = computeModules(bipart_twc_all, deep = FALSE, deleteOriginalFiles = TRUE,
                                    steps = 1000000, tolerance = 1e-10, experimental = FALSE)
###Function "plotModuleWeb"
#This function takes an object of class "moduleWeb" and plots the modules found by function
#computeModules(...) onto the graph.
plotModuleWeb(modules_bipart_twc, plotModules = TRUE,
              rank = FALSE, weighted = TRUE, displayAlabels = TRUE,
              displayBlabels = TRUE, labsize = 1, xlabel = "", ylabel = "",
              square.border = "white", fromDepth = 0, upToDepth = -1)

######wine: Weighted-Interaction Nestedness Estimator
#Calculates the nestedness of a network taking into account the weight of the interactions, according
#to the method proposed by Galeano et al. (2008).

bipart_twc_wine = wine(bipart_twc_all, nreps = 1000)
bipart_twc_wine
plot.wine(bipart_twc_wine)
write.csv(bipart_twc_wine$dij.w, file="bipart_twc_wine_dijw.csv")

##############NULL MODELS#############

#null.distr: Null model based on fitted marginal distribution
#Given a network, this function fits a distribution to the marginal totals and then draws randomly
#from this distribution to yield a new network.

null.distr(N=2, bipart_twc_all, distr="lognormal")

####nullmodel: Generates null models for network analysis
#Methods 1-4 for Quantitative Bipartite Graphs 
#Methods 4-5 for Binary Bipartite graphs
#method: Null model type. Can be given as an integer or name: 1/"r2dtable", 2/"swap.web",
#3/"vaznull", 4/"shuffle.web", 5/"mgen"; allows for partial match of names; methods
#1 to 4 works for quantitative webs, 4 and 5 for binary.

data(bipart_twc_all)
nullmodel(bipart.urban.region, N=2, method=1)
nullmodel(bipart.urban.region>0, N=2, method=4)
obs <- unlist(networklevel(bipart.urban.region, index="generality"))
nulls <- nullmodel(bipart.urban.region, N=100, method=3)
null <- unlist(sapply(nulls, networklevel, index="generality")) #takes a while ...
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="generality: Region")
abline(v=obs, col="red", lwd=2)
praw <- sum(null>obs) / length(null)
ifelse(praw > 0.5, 1-praw, praw) # P-value

#vaznull: Null model with constrained totals and connectance
#Implements Diego Vazquez proposal of a null model for pollination networks. "The algorithm randomized the total number of individual
# interactions observed in the original interaction matrix, F."
vaznull(1000, bipart.urban.region)

###null.t.test: Compares observed pattern to random web using r2dtable
#A little null-model function to check, if the observed values actually are much different to what one
#would expect under random numbers given the observed row and column totals.

null.t.test(bipart.urban.region, index=c("generality", "vulnerability","weighted connectance", "H2", "ISA", "SA"), nrep=2, N=1000)
