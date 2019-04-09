require(MASS)

##read in file:
d <- read.csv("data/old_data/2013_Neighbor_survey_prepped.csv")


###############################################
###############################################
## set up a few key data subsets here:

## plants from plots designated as lambda plots:
lam <- subset(d, d$density=="lambda")

## drop lambda rows without seed production:
lam<-lam[!is.na(lam$seeds),]

## intermediate subset- plants from plots designated as competition plots ('not lambda')
nl <- subset(d, d$density!="lambda")

## drop rows with no seed production from the not lambda subset:
nl<-nl[!is.na(nl$seeds),]

## set up a df of just plants in competition with seed production (keep this):
comp <- subset(nl, nl$neighbours_number>0)

## pull out rows from 'not lambda' plots without backgrounds:
extra_zeros <- subset(nl, nl$neighbours_number==0)

## add these extra no-background plants into the lambda df:
lam <- rbind(lam, extra_zeros)

##two df to use at this point- 'comp' and 'lam'


################################################
################################################

## log likelihood functions, modeled on those Janneke wrote:

## while weren't interested in the simpler models per se,
## the way that Janneke sets these nested models up uses the parameter estimates
## from the simpler models as the starting parameters of the more complex models
## which should be helpful in getting the estimation for the more complex models to converge properly

## OJO!!! these functions call data objects that aren't passed to them as parameters like 'log_seed'

#model 1 - no effect of density (no competitive effects)
compmodel1<-function(par){
	
	## mu (=lambda later) and sigma parameters for the normal distribution (assuming lognormal error- seed data are logged)
	lambda<-par[1]
	sigma<-par[2]
	
	#this the predictive model- here is just fitting a horizontal line through the data:
	pred<-rep(lambda, times=length(log_seeds)) 
	
	#these are the log likelihoods of the data given the model + parameters
	llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE) 
	
	#return the sum of negative log likelihoods - what optim minimizes
	return(sum(-1*llik)) 
}

#model 2 - competition, but no difference between species
compmodel2<-function(par){
	
	lambda<-par[1] ## same as model 1
	alpha<-par[2]  ## new parameter introduced in model 2
	sigma<-par[3] ## same as model 1
	
	## apply is used to get a single vector of densities rather than a matrix
	## there is only one competitor at a time here- just adding each density to a 
	## column of other zeros:
	
	tot_density<-apply(dens, MARGIN=2, FUN=sum)
	
	## predictive model:
	pred<-lambda/(1+alpha*(tot_density)) 
	
	## log likelihoods of data given the model + parameters:
	llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
	
	## return sum of negative log likelihoods:
	return(sum(-1*llik)) 
	
}


#model 3 - all species have different competitive effects
compmodel3<-function(par){
	lambda<-par[1] #same as model 2
	a_AGHE<-par[2]	## new parameters- use alpha estimate from model 2 as start value for fitting
	a_AGRE<-par[3]
	a_AMME<-par[4]
	a_ANAR<-par[5]
	a_CEME<-par[6]
	a_CLPU<-par[7]
	a_ERBO<-par[8]
	a_ERCI<-par[9]
	a_EUPE<-par[10]
	a_GECA<-par[11]
	a_HECO<-par[12]
	a_LACA<-par[13]
	a_LOPU<-par[14]
	a_LOWR<-par[15]
	a_MEPO<-par[16]
	a_NAAT<-par[17]
	a_PLER<-par[18]
	a_SACA<-par[19]
	a_SIGA<-par[20]
	sigma<-par[21] ## same as model 2
	
	##probably overkill to name all the alphas, but it helps keep things straight
	
	# same form as model 1, but super long given all the species:
	pred<- lambda/ (1+ a_AGHE* dens['AGHE',] + a_AGRE * dens['AGRE',] + a_AMME* dens['AMME',] + a_ANAR* dens['ANAR',] + a_CEME* dens['CEME',] + a_CLPU* dens['CLPU',] + a_ERBO* dens['ERBO',] + a_ERCI* dens['ERCI',] + a_EUPE* dens['EUPE',] + a_GECA* dens['GECA',] + a_HECO* dens['HECO',] + a_LACA* dens['LACA',] + a_LOPU* dens['LOPU',] + a_LOWR* dens['LOWR',] + a_MEPO* dens['MEPO',] + a_NAAT* dens['NAAT',] + a_PLER* dens['PLER',] + a_SACA* dens['SACA',] + a_SIGA* dens['SIGA',] )
	
	
	# likelihood as before:
	llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
	
	# return sum of negative log likelihoods
	return(sum(-1*llik)) #sum of negative log likelihoods
}





################################################







## get a list of target species to work through sequentially:
splist<-unique(comp$target)
## drop CLBO, which has no lambda plants
splist<-paste(splist[-(which(splist=="CLBO")), drop=TRUE])

##alpha order splist:

splist<-splist[order(splist)]



## objects to hold the final parameter estimates from model 3:

alpha_matrix<- matrix(0, nrow=length(splist), ncol=length(splist))
alpha_SE_matrix<- matrix(0, nrow=length(splist), ncol=length(splist))
N_away_0<-matrix(0, nrow=length(splist), ncol=length(splist))


row.names(alpha_matrix)<-splist
colnames(alpha_matrix)<-splist
lambda_est<-NULL
lambda_SE<-NULL
sigma_est<-NULL
convergence_code<-NULL


pdf(file="figures/alphas_and_lambdas_single_model.pdf", width=11, height=8)


## for each species in turn as a target:
for(i in 1:length(splist)){
	
	#set up grid for plotting and add first panel indiciating focal sp:
	layout(matrix(1:21, nrow=3))
	plot(x=0, y=0, type='n')
	text(0, 0, paste(splist[i], "as focal", sep=" "), cex=1, col="red")
	
	
	## subset out the rows that are needed from the competition df
	comp_points<-subset(comp, comp$target==splist[i], drop=TRUE)
	
	## and the correct lambda plants
	
	lam_points<-subset(lam, lam$target==splist[i])
	
	## now need to build up a vector of nonzero seed production
	## and corresponding density vectors (held in a matrix) for background species 
	## to use in model fitting
	
	## start with the lambda seeds- will add on to this for each background:
	
	seeds<-lam_points$seeds
	
	
	##build density matrix (each row will be density of a species)
	dens<-matrix(0, nrow=length(splist), ncol=(nrow(lam_points)+ nrow(comp_points)))
	
	row.names(dens)<-splist
	
	
	##for each background species the target competes against:
	background_list<-unique(comp_points$background)
	
	## use this counter to keep track of which column is next to have data added to
	## set it to begin after the lambda points:
	
	start <- nrow(lam_points) +1
	
	
	for(j in 1:length(background_list)){
		
		## LOOP creates a seeds vector and a corresponding dens ("density") matrix with 
		## background species as rows, and columns corresponding to the seed production
		## vector.
		## beginning columns correspond to lambda plants
		
		#take just the rows pertaining to a specific background sp:
		bg_points<- subset(comp_points, comp_points$background==background_list[j])
		
		## which row of the density matrix corresponds?
		rownum <- which(row.names(dens)==background_list[j])
		
		#column to end with:
		
		end<-start + nrow(bg_points)-1
		
		## drop in density values into matrix:
		
		dens[rownum, start:end]<-bg_points$neighbours_number
		
		## add seed numbers into the seeds vector
		
		seeds<-c(seeds, bg_points$seeds)
		
		
		
				
		##increase start counter so that next species is offset from this one
		start<-end+1
		
		
		
		
	}
	
	## should now have "seeds" vector and corresponding "dens"ity matrix
	## can test ncol(dens)==length(seeds)
	
	## we'll be working with log seeds (for giving a lognormal error structure):
		
	log_seeds<-log(seeds)
	
	
	#model fitting using optim and earlier likelihood functions
	
	#############################
	## model 1, no competition ##
	#############################
	
	###recall parameters are lambda and sigma- initialize these with estimates from the data:
	par1<-c(mean(log_seeds), sd(log_seeds))
	
	##repeat optimization until we get convergence (or we try 25 times)
	for(k in 1:25){
		testcomp1<-optim(par1,compmodel1)
		##update start parameters to final estimate to use in next run in case of nonconvergence
		par1<-testcomp1$par
		if(testcomp1$convergence==0){
			print(paste(splist[i], "model 1 converged on rep", k, sep=" "))
			break
			}
		
		}

		
	#############################
	## model 2, one alpha      ##
	#############################
	
	## pars here are lambda, alpha and sigma- use lambda and sigma from model 1 as starting esimtates
	
	par2<-c(testcomp1$par[1],0.001,testcomp1$par[2])
	
	##as before:
	for(k in 1:25){
		##now using a specific method that permits constrained optimization so that alpha has to be nonzero- this is an issue in some of the fits, especially in model 3. lower parameter has same order as par2
		testcomp2<-optim(par2,compmodel2, method="L-BFGS-B", lower=c(1,0,0.0000000001))
		par2<-testcomp2$par
		if(testcomp2$convergence==0){
			print(paste(splist[i],  "model 2 converged on rep", k, sep=" "))
			break
			}
		}

	#############################
	## model 3, unique alphas  ##
	#############################
	
	## pars here are lambda, 19 alphas and sigma- use lambda and sigma from model 2, and alphas from model 2 as starting esimtate for each species
	
	par3<-c(testcomp2$par[1], rep(testcomp2$par[2], times=19), testcomp2$par[3])
	
	##as before
	for(k in 1:25){
		##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
		testcomp3<-optim(par3,compmodel3, method="L-BFGS-B", lower=c(1, rep(0.001, times=19),0.0000000001), control=list(maxit=1000), hessian=TRUE)
		par3<-testcomp3$par
		if(testcomp3$convergence==0){
			print(paste(splist[i], "model 3 converged on rep", k, sep=" "))
			break
			}
		}

	
	##################################
	### save estimates from model 3 ##
	##################################
	
	lambda_est<-c(lambda_est, testcomp3$par[1])
	lambda_SE<-c(lambda_SE, sqrt(ginv(testcomp3$hessian)[1,1]))
	sigma_est<-c(sigma_est, testcomp3$par[21])
	convergence_code<-c(convergence_code, testcomp3$convergence)
	
	
	## in keeping with Lotka Volterra notation, we'll use alpha1_2 to indicate effect of
	## sp 2 on growth of 1.  Following convention, i refers to rows and j to cols in a matrix
	## so each step of the loop here (for a target sp) corresponds to one row of this matrix:
	
	alphas<-testcomp3$par[2:20]
	
	alpha_matrix[i,]<-alphas
	
	alpha_SE_matrix[i,]<-sqrt(diag(ginv(testcomp3$hessian))[2:20])
	
	
	
	##note that in cases where there is no data for a particular species the alpha estimate for that species ends up as the starting value- we need to be careful of these as they are basically gargbage numbers.  Keeping them in up to now to keep the structure of the data constant, but will set them to NA here:
	
	#identify which species have no data in this fit:
	which(apply(dens, MARGIN=1, FUN=mean)==0)->no_data
	
	pres_abs<-dens>0
	

	
	
	N_away_0[i,]<-apply(pres_abs, MARGIN=1, FUN=sum)
	
	## set their alphas to NA in the matrix:
	
	alpha_matrix[i,no_data]<-NA
	
	
	
	
	
	###############################
	## some diagnostics +  plots ##
	###############################
	
	## print an error to the console if any one of the three models failed to converge:
	
	if(testcomp1$convergence + testcomp2$convergence + testcomp3$convergence !=0){print(paste("at least one model did not converge for", splist[i], sep=" "))}
	
	##################################
	## plot observed vs predicted:
	##################
	
	par<-testcomp3$par
	
	#from model 3 code:
	lambda<-par[1] 
	a_AGHE<-par[2]	
	a_AGRE<-par[3]
	a_AMME<-par[4]
	a_ANAR<-par[5]
	a_CEME<-par[6]
	a_CLPU<-par[7]
	a_ERBO<-par[8]
	a_ERCI<-par[9]
	a_EUPE<-par[10]
	a_GECA<-par[11]
	a_HECO<-par[12]
	a_LACA<-par[13]
	a_LOPU<-par[14]
	a_LOWR<-par[15]
	a_MEPO<-par[16]
	a_NAAT<-par[17]
	a_PLER<-par[18]
	a_SACA<-par[19]
	a_SIGA<-par[20]

	pred<- lambda/ (1+ a_AGHE* dens['AGHE',] + a_AGRE * dens['AGRE',] + a_AMME* dens['AMME',] + a_ANAR* dens['ANAR',] + a_CEME* dens['CEME',] + a_CLPU* dens['CLPU',] + a_ERBO* dens['ERBO',] + a_ERCI* dens['ERCI',] + a_EUPE* dens['EUPE',] + a_GECA* dens['GECA',] + a_HECO* dens['HECO',] + a_LACA* dens['LACA',] + a_LOPU* dens['LOPU',] + a_LOWR* dens['LOWR',] + a_MEPO* dens['MEPO',] + a_NAAT* dens['NAAT',] + a_PLER* dens['PLER',] + a_SACA* dens['SACA',] + a_SIGA* dens['SIGA',] )	
	
	
	plot(seeds, pred, xlim=c(min(c(seeds, pred)), max(c(seeds, pred)) ), ylim=c(min(c(seeds, pred)), max(c(seeds, pred)) ), log='xy', xlab="observed seeds", ylab="predicted seeds", main=splist[i] )
	
	abline(a=0, b=1, lwd=2)
	
	
	#####################
	#### plot each fit
	##########
	
	names(alphas)<-splist
	
	names(apply(dens, MARGIN=1, FUN=mean)!=0)->plotlist
	for(l in 1:length(plotlist)){
		
		## which columns in the density dataframe have nozero values for species l ?
		cols<-which(dens[l,]>0)
		
		
		x<-dens[l,cols]
		y<-seeds[cols]
		
		
		##add lambdas:
		x<-c(x, rep(0,nrow(lam_points)))
		y<-c(y, lam_points$seeds)
		
		x_det<-seq(min(x), max(x), by=(  (max(x)-min(x))/1000  ))
		
		alpha_temp<-alphas[which(names(alphas)==plotlist[l])]
		y_pred<-lambda/(1 + alpha_temp*x_det)
		
		if(length(unique(x))>1){
			plot(y~x, xlab="density", ylab="seeds", main=plotlist[l])
			} else {
				
				plot(x=0, y=0, main=plotlist[l], type='n')
				text(0, 0, "no data")
				
			}
		
		lines(y_pred~x_det, col="red", lwd=2)	
	}

	
	
	
	
	

}


dev.off()

results<-data.frame(splist, lambda_est, lambda_SE, sigma_est, convergence_code)

names(results)<-c("species", "lambda", 'lambda_SE', "sigma", "convergence_code")

write.csv(results, "output/lambda_estimates.csv")
write.csv(alpha_matrix, "output/alpha_estimates_row_is_target.csv")

row.names(N_away_0)<-splist
colnames(N_away_0)<-splist
write.csv(N_away_0, "output/N_away_0_row_is_target.csv")

row.names(alpha_SE_matrix)<-splist
colnames(alpha_SE_matrix)<-splist
write.csv(alpha_SE_matrix, "output/alpha_SE_row_is_target.csv")
