#######################
### Code for presence/absence Beta-null deviation calculations.
### Helper functions for code "ExampleSimulationNullDevCalc.R"
### Prepared May 14, 2014 
### Authors Brett Melbourne, Caroline Tucker, Lauren Shoemaker 
#######################

## This is a modularized version of Chase et al's raup_crick function. It now calculates the null deviation in a range of ways. The null distributions take the longest, so these can be calculated and stored in an object, and saved to disk. From these stored distributions the null deviations can be calculated using a range of methods.

## Other things that were changed from Chase et al's raup_crick function spXsite became siteXsp

###----siteXsp_prep()----
# Helper function to clean up the siteXsp matrix. Removes names from column 1 if any, removes redundant species columns, and converts to presence-absence.
# siteXsp	matrix for siteXsp, with row names for plots in first column, or optionally no names.


siteXsp_prep <- function(siteXsp, plot_names_in_col1=TRUE){

	# Move plot names in column 1 (if specified as being present) into the 
	# row names of the matrix and drop the column of names
	if(plot_names_in_col1){
		row.names(siteXsp) <- siteXsp[,1]
		siteXsp <- siteXsp[,-1]
	}
	
	# Remove any redundant species columns (i.e. species that were never observed).
	colsums <- colSums(siteXsp)
	siteXsp <- siteXsp[,colsums > 0]

	# Make the siteXsp matrix into a pres/abs. (overwrites intitial siteXsp 
	# matrix):
	siteXsp <- ifelse(siteXsp > 0, 1, 0)

	return(siteXsp)
}



####----null_distributions()----
# Generates null distributions for the chosen test function.

# siteXsp	Matrix for siteXsp. This must be a presence-absence matrix fully prepared for use. Use siteXsp_prep() first.
# occur		Optional vector of occurrences (number of sites each species was observed in) for the global species pool. Must not contain zeros.
#test_func	Function used to compare sites, e.g. number of shared species or Jaccard. The test function can return a matrix or a dist object (which allows any of the distance functions in dist or vegan's vegdist or betadiver to be 			used). Only the lower triangular needs entries for a siteXsite matrix (the rest can be NA). Distance functions are considerably slower due to converting to a distance object.


# Value Returns a list of null distributions. Each null distribution is a numeric vector. There is one distribution for each pairwise combination of alpha diversity.
##


null_distributions <- function(siteXsp, occur = NULL, reps=9999,test_func,...){

	# If occurrence vector is not supplied, create one from siteXsp matrix.
	if (is.null(occur)){
		occur <- colSums(siteXsp)
	}

	# Calculate gamma, the global species richness.
	gamma <- length(occur)

	# Count number of sites
	n_sites <- nrow(siteXsp)

	# Determine how many unique species richness values are in the dataset.
	# This is used to limit the number of null communities that have to be 
	# calculated.
	alpha_levels <- sort(unique(rowSums(siteXsp)))
	
	# Make_null:
	
	# alpha_table is used as a lookup to help identify which null distribution
	# to use for the tests later.  It contains one row for each combination 
	# of alpha richness levels. 
	
	alpha_table <- data.frame(c(NA), c(NA))
	names(alpha_table) <- c("smaller_alpha", "bigger_alpha")
	
	# null_array will hold the actual null distributions.  Each element
	# of the array corresponds to a distribution for one combination 
	# of alpha values.  The alpha_table is used to point to the correct null
	# distribution - the row numbers of alpha_table correspond to the [[x]] 
	# indices of the null_array.  Later the function will find the row of 
	# alpha_table with the right combination of alpha values.  That row 
	# number is used to identify the element of null_array that contains the
	# correct null distribution for that combination of alpha levels. 
	
	null_array <- list()

	# Calculate total loops for monitoring
	total_loops <- sum(1:length(alpha_levels))

	# For each combination of alpha levels:
	col_count <- 1
	for(a1 in 1:length(alpha_levels)){
		for(a2 in a1:length(alpha_levels)){
			
			# Generate random communities for pairs of alpha values and
			# calculate a null distribution for the test function.
			null_distr <- rep(NA,reps)
			for(i in 1:reps){
				
				# two empty null communities of size gamma:
				com1 <- rep(0,gamma)
				com2 <- rep(0,gamma)
				
				# add alpha1 number of species to com1, weighting by 
				# species occurrence frequencies:
				com1[sample(1:gamma, alpha_levels[a1], replace=FALSE, prob=occur)] <- 1
				
				# same for com2:
				com2[sample(1:gamma, alpha_levels[a2], replace=FALSE, prob=occur)] <- 1	

				# Calculate the test function for the pair
				# as.matrix allows dist objects to be used - slower.
				null_pair <- t(cbind(com1,com2))
				null_distr[i] <- as.matrix( test_func(null_pair,...) )[1,2]

			}
			
			# Store null distribution, record values for alpha 1 and 2 in the alpha_table to help 
			# find the correct null distribution later:
			null_array[[col_count]] <- null_distr
			
			alpha_table[col_count, which(names(alpha_table)=="smaller_alpha")] <- alpha_levels[a1]
			alpha_table[col_count, which(names(alpha_table)=="bigger_alpha")] <- alpha_levels[a2]

			# increment the counter for the columns of the alpha table/ elements of the null array
			col_count <- col_count+1
			
		}
	
	}
	
	# create a new column with both alpha levels to match on:
	alpha_table$matching <- paste(alpha_table[,1], alpha_table[,2], sep="_")
	
# In future this could also return the testfunc name here too. Don't know how. Need to extract
# the call to null_deviations. This could then be accessed by null_deviations() without the need
# to specify again in arguments.
	return( list(null_array=null_array,alpha_table=alpha_table,
			 siteXsp=siteXsp, occur=occur) )
}


#----n_shared()----
# One option for test_func in null_distributions() and null_deviations(). 
# Calculates a matrix of the number of shared species for pairs of sites.
# This test_func in combination with pval_1tail_L gives the classic Raup-
# Crick measure. Thanks to betadiver in vegan for the tcrossprod solution.
#
n_shared <- function(siteXsp,asdist=FALSE){
	n_shared <- tcrossprod(siteXsp)
	if (asdist==FALSE) return( n_shared )
	return( as.dist(n_shared) )
}



###----null_deviations----
## Calculates the null deviations by the chosen method. The function null_distributions should be called first using the global siteXsp matrix to assemble the null distributions and store in an object. The classic Raup-Crick measure is given by using the pval_1tail_L method with n_shared as the test_func. The Chase et al. modified Raup-Crick measure is given by using the chase_L method with n_shared() as the test_func. RC and modified RC respectively measure the probability of a beta diversity as large or larger than the observed (actually the probability of the number of shared species as small or smaller than observed), given the null model is true. Since the number of shared species is used, it is the left tail that is wanted to test the hypothesis of high beta diversity. When a dissimilarity or beta diversity is used directly as the test function, the right tail will be needed to test the same hypothesis of high beta diversity.

# siteXsp	 		Matrix for siteXsp. This must be a presence-absence matrix fully prepared for use. Use siteXsp_prep() first. This can be different to the global siteXsp matrix, such as a subset of an experiment, 						because all we are going to do is compare pairs of sites in this matrix with pairs in the null that share the same species richness.
# null_object     The object returned from the call to null_distributions.
# test_func	 	See null_distributions().
# devmethod       Method to calculate the deviation:
	#			chase_L		- Chase et al (left tail)
	#			chase_R		- Chase et al (right tail)
	#			pval_1tail_L	- Classic Raup-Crick style (left tail)
	#			pval_1tail_R	- Classic Raup-Crick style (right tail)
	#			pval_2tail	- Two tailed p-value
	#			lnlik		- ln likelihood
	#			eval		- Expected deviation
	#			eval_rel	- Relative expected deviation
#...		 Parameters passed to test_func. See null_distributions().
##


null_deviations <- function(siteXsp,null_object,test_func,...){

	alpha_table <- null_object$alpha_table
	null_array <- null_object$null_array
	reps <- length(null_array[[1]])
	null_object <- NULL #Clean out the large object

	# Number of sites
	n_sites <- nrow(siteXsp)	

	# Matrix of observed values for the test function. Only the lower 
	# triangular of this is needed (so we delete the rest to match the
	# other returned matrices. as.matrix allows dist objects to be used.
	testfunc_obs <- as.matrix(test_func(siteXsp,...))
	diag(testfunc_obs) <- NA
	testfunc_obs[upper.tri(testfunc_obs)] <- NA

	# Initialize site by site matrices for the results, with the names of the 
	# sites in the row and col names (only lower tri is used).

## Chase L
	null_devs_Chase_L<- matrix(data=NA, nrow=n_sites, ncol=n_sites, 
				dimnames=list(row.names(siteXsp), row.names(siteXsp)))


### pval 1tail L
	null_devs_pval_1tail_L <- matrix(data=NA, nrow=n_sites, ncol=n_sites, 
				dimnames=list(row.names(siteXsp), row.names(siteXsp)))


### eval	
	null_devs_eval <- matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(siteXsp), row.names(siteXsp)))
	enulls <- null_devs_eval


# More to keep
obs_a1_keep <- null_devs_eval
obs_a2_keep <- null_devs_eval
null_index_keep <- null_devs_eval

	# Calculate deviation for each pair of sites. Keep lower triangular only.
	for(i in 1:(n_sites-1)){
		for(j in (i+1):n_sites){
			
			# What was the observed richness of each site?
			obs_a1 <- sum(siteXsp[i,])
			obs_a2 <- sum(siteXsp[j,])
	
			# Place these alphas into an object to match against 
			# alpha_table (sort so smaller alpha is first)
			obs_a_pair <- sort(c(obs_a1, obs_a2))
			
			# Match against the alpha table- row index identifies which 
			# element of the null array contains the correct null 
			# distribution for the observed combination of alpha values:
			null_index <- which(alpha_table$matching==paste(obs_a_pair[1], obs_a_pair[2], sep="_"))

			# What is the expected value of the null distribution?
			expected <- mean(null_array[[null_index]])

			# How many null observations are equal, greater, or less than
			# the observed value?
			num_equal_in_null <- sum( null_array[[null_index]] == testfunc_obs[j,i] )
			num_greater_in_null <- sum( null_array[[null_index]] > testfunc_obs[j,i] )
			num_less_in_null <- sum( null_array[[null_index]] < testfunc_obs[j,i] )

			
			# Calculate different measures of the null deviation
			value_Chase_L <- 1 - ( ( num_greater_in_null + num_equal_in_null / 2 ) / reps )
			
  			# Classic Raup-Crick
			value_pval_1tail_L <- 1 -(num_greater_in_null)/reps

			# eval
			value_eval <- testfunc_obs[j,i] - expected

			# store the results:
			null_devs_Chase_L[j,i] <- value_Chase_L

			null_devs_pval_1tail_L[j,i] <- value_pval_1tail_L

			null_devs_eval[j,i] <- value_eval


obs_a1_keep[j,i] <- obs_a1
obs_a2_keep[j,i] <- obs_a2
null_index_keep[j,i] <- null_index

		}
	}
	
return(list(null_devs_Chase_L=null_devs_Chase_L,null_devs_pval_1tail_L=null_devs_pval_1tail_L,null_devs_eval=null_devs_eval,null_distrbns=null_array,enulls=enulls,testfunc_obs=testfunc_obs,obs_a1=obs_a1_keep,obs_a2=obs_a2_keep,null_index=null_index_keep,alpha_table=alpha_table))
}


