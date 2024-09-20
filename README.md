# MSOM
Simulation and model fitting code to accompany "Sample size considerations for species co-occurrence models".
To run the code smoothly, we recommend opening the project ".RProj" object and run scripts from there. 

R Script(s) descrption:
# Data simulation
 - Co-occurrence simulation_v2.R: Code to simulate and fit data for the 13 different scenarios increasing in model complexity
	- 6 Null model with two species
	- 4 With two species and increasing covariate complexity
	- 3 With increasing number of species (beware of computing time!)
 - MSOM_SimFun.R: Code for the smulation and fit of data and models ('MSOM_simfit.fun.v2') as well as a functions to do cross-validation
		  and to derivate general parameters from natural parameters (psi.fun). 
		  Previous function to simulate and fit data 'MSOM_SimandFit.fun' in code too. 

# Result processing
 - TwoSpNullModelsResults.R: Code to process the results for all the null models and create associated figures/tables
 - CovariateModelResults.R: Code to process the results for all the covariate complexity models and create associated figures/tables
 - NSpeciesModelResults.R: Code to process the results for all the network complexity models and create associated figures/tables

# Result object description
 - When you simulate a model (e.g. "test" line 70 'Co-occurrence simulation_v2.R') you can access individual data.frames using the `$`:
	- "Det" will be for detecton and "State" for occupancy
	- ".pl" indicates it has been fit with penalised likelihood
	- "time.ellapsed" is how long the function was running for.

# Other Files/"folders"

 - "Results": folder containing the results objects ("rds" and "csv") for all simulations with self-explanatory titles
 - Co-occurrence simulation.R: Old version to simulate code for scenarios (4 scenarios) using the MSOM_SimandFit.fun function. (_Deprecated from 30 Apr 2024 onwards_) 


