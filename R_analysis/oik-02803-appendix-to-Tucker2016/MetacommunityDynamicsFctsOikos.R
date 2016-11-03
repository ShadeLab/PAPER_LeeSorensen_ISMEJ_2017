#######################
### Code for deterministic and stochastic metacommunity models 
### Helper functions for code "ExampleSimulationNullDevCalc.R"
### Prepared May 14, 2014 
### Authors Lauren Shoemaker, Caroline Tucker, Brett Melbourne
#######################

## Code for the metacommunity and patch dynamics functions described in Tucker et al. 
## Required for Beta-null deviation simulations

## Function definition for Beverton-Holt competition model
bevHoltC <- function(N,a) {
	return( 1 / ( 1 + a*(sum(N))))
}

## Define birth rate (lambda) functions
# Function definition for deterministic birth
b <- function(R,current, patches, species, a) {
	birth <- matrix(NA, nrow=species, ncol=patches)
	for (p in 1:patches) {
		birth[,p] <- R[,p]*current[,p]*bevHoltC(current[,p],a)
	}
	return(birth)
}

# Function definition for stochastic birth
bstoch <- function(R,current, patches, species, a) {
	birth <- matrix(NA, nrow=species, ncol=patches)
	for (p in 1:patches) {
		birth[,p] <- rpois(species,R[,p]*current[,p]*bevHoltC(current[,p],a))
	}
	return(birth)
}


## Define metacommunity simulations. 
# Function definition for deterministic model run with global dispersal
rdet <- function(R,d,a,patches,species,data) {
	
	for (t in 1:(time-1)){
		
		# birth
		
		birth <- b(R, data[,,t],patches, species, a)
		
		# dispersal
		migrate <- d*birth
		stay <- birth - migrate
		
		totalmigrate <- rep (NA, species)
		for (s in 1:species) {
			totalmigrate[s] <- sum(migrate[s,])
		}
		
		# determine where each individual migrates
		migratelocation <- matrix(NA, nrow=species, ncol=patches)
		for (s in 1:species) {
			migratelocation[s,] <- rep(totalmigrate[s]/patches, patches)
		}
		
		# Determine number of individuals in each patch
		data[,,t+1] <- stay + migratelocation	
		
	}
return(data)
}


# Function definition for stochastic model run with global dispersal 
rstoch <- function(R,d,a,patches,species,data) {

	for (t in 1:(time-1)){
		
		# birth	
		birth <- bstoch(R, data[,,t],patches, species, a)
		
		# dispersal
		migrate <- matrix(rbinom(patches*species,birth,d), nrow=species, ncol=patches)
		stay <- birth - migrate
		
		totalmigrate <- rep (NA, species)
		for (s in 1:species) {
			totalmigrate[s] <- sum(migrate[s,])
		}
		
		# determine where each individual migrates
		migratelocation <- matrix(NA, nrow=species, ncol=patches)
		for (s in 1:species) {
			patchLocation <- floor(runif(totalmigrate[s], min=1, max=patches+1))
			for (p in 1:patches) {
				migratelocation[s,p] <- length(which(patchLocation==p))
			}
		}
		
		# Determine number of individuals in each patch
		data[,,t+1] <- stay + migratelocation	
		
	}
return(data)
}


