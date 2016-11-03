#######################
### Code for example metacommunity simulation and beta-null deviation calculations
### with "Differentiating between niche and neutral assembly in metacommunities using
### null models of beta-diversity"
### Prepared May 14, 2014 
### Authors Caroline Tucker, Lauren Shoemaker, Brett Melbourne
#######################


## Load required source files and libraries
library(reldist)
library(vegan)
library(bipartite)
source("MetacommunityDynamicsFctsOikos.r")
source("PANullDevFctsOikos.r")

## Set metacommunity size, time length
patches <- 25  # Number of patches
species <- 25  # Number of species
time <- 150	   # Length of model run (generations)


## Set parameters for metacommunity dynamics
Rbest <- rep(NA, species)	# intrinsic growth-rate patch type 2 (optimal patch)
d <- rep(0.05, species)		# migration rate

Rbest <- rep(1.45, species) # Species sorting parameterization

a <- 1/600					# Beverton-holt alpha for all species. This is fixed.

R <- matrix(1.095833, nrow=species, ncol=patches) # intrinsic growth-rate patch type 1 (non-optimal patch)
diag(R) <- Rbest

# Data array for storing all results
dat <- array(NA, c(species, patches, time))
dat[,,1] <- rep(150, species*patches)	# starting conditions--150 species in each patch

## Stochastic metacommunity - run for 150 years
stochasticRun <- rstoch(R, d, a, patches, species, dat) 
stoch <- stochasticRun[,,time] # use the last time point for null deviation analyses
stoch <- ifelse(stoch[,] <2, 0, stoch[,]) # threshold number of individuals

dune<-t(stoch)

## OR, deterministic metacommunity - run for 150 years
# deterministicRun <- rdet(R, d, a, patches, species, dat)
# det <- deterministicRun[,,time]
# det <- ifelse(det[,] <2, 0, det[,])

# dune<-t(det)

################

##Calculate beta-diversity for metacommunity
dune_bin <- siteXsp_prep(dune, plot_names_in_col1=FALSE)

beta_comm <- vegdist(dune_bin, "jaccard")  	# Calculate beta-diversity
res_beta_comm <- as.matrix(as.dist(beta_comm))
diag(res_beta_comm) <- NA
beta_div_stoch <- apply(res_beta_comm, 2, FUN="mean", na.rm=TRUE) 	# Calculate patch mean value

## Calculate PA null deviation
nullobj <- null_distributions(dune_bin, test_func = vegdist, method = "jaccard", reps=999)		# generate null distributions
nulldev <- null_deviations(dune_bin, nullobj, vegdist, method="jaccard")		# Calculate null deviations using the expected value method

### Store results, where each value is the PA beta-null deviation value for a pairwise comparison (distance matrix) between patches
res_nulldev <- as.matrix(as.dist(nulldev$null_devs_eval))
diag(res_nulldev) <- NA
PA_null_dev <- apply(res_nulldev, 2, FUN="mean", na.rm=TRUE) #Calculate patch mean value


### Prepare and calculate abundance beta-null deviation metric
## Adjusted from Stegen et al 2012 GEB
bbs.sp.site <- dune
rand <- 999
null.alphas <- matrix(NA, ncol(dune), rand)
null.alpha <- matrix(NA, ncol(dune), rand)
expected_beta <- matrix(NA, 1, rand)
null.gamma <- matrix(NA, 1, rand)
null.alpha.comp <- numeric()
bucket_bray_res <- matrix(NA, patches, rand)

bbs.sp.site = ceiling(bbs.sp.site/max(bbs.sp.site)) 
mean.alpha = sum(bbs.sp.site)/nrow(bbs.sp.site) #mean.alpha
gamma <- ncol(bbs.sp.site) #gamma
obs_beta <- 1-mean.alpha/gamma
obs_beta_all <- 1-rowSums(bbs.sp.site)/gamma

##Generate null patches
for (randomize in 1:rand) {  
	null.dist = dune
	for (species in 1:ncol(null.dist)) {
		tot.abund = sum(null.dist[,species])
		null.dist[,species] = 0
		for (individual in 1:tot.abund) {
			sampled.site = sample(c(1:nrow(bbs.sp.site)), 1)
			null.dist[sampled.site, species] = null.dist[sampled.site, species] + 1
		}
	}
	
	##Calculate null deviation for null patches and store
	null.alphas[,randomize] <- apply(null.dist, 2, function(x){sum(ifelse(x > 0, 1, 0))})
	null.gamma[1, randomize] <- sum(ifelse(rowSums(null.dist)>0, 1, 0))
	expected_beta[1, randomize] <- 1 - mean(null.alphas[,randomize]/null.gamma[,randomize])
	null.alpha <- mean(null.alphas[,randomize])
	null.alpha.comp <- c(null.alpha.comp, null.alpha)
	
	bucket_bray <- as.matrix(vegdist(null.dist, "bray"))
	diag(bucket_bray) <- NA
	bucket_bray_res[,randomize] <- apply(bucket_bray, 2, FUN="mean", na.rm=TRUE)
	
} ## end randomize loop

## Calculate beta-diversity for obs metacommunity
beta_comm_abund <- vegdist(dune, "bray")
res_beta_comm_abund <- as.matrix(as.dist(beta_comm_abund))
diag(res_beta_comm_abund) <- NA
# output beta diversity (Bray)
beta_div_abund_stoch <- apply(res_beta_comm_abund, 2, FUN="mean", na.rm=TRUE)

# output abundance beta-null deviation
abund_null_dev <- beta_div_abund_stoch - mean(bucket_bray_res)


### Outputs:

#beta_div_stoch  - Jaccard beta-diversity for the metacommunity, average value (of all pairwise comparisons) for each patch
#beta_div_abund_stoch - Bray-Curtis beta-diversity for the metacommunity, average value (of all pairwise comparisons) for each patch
#PA_null_dev - presence-absence null deviation values or the metacommunity, average value (of all pairwise comparisons) for each patch
#abund_null_dev - abundance null deviation values or the metacommunity, average value (of all pairwise comparisons) for each patch

