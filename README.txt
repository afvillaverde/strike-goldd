This is the README file of the STRIKE-GOLDD 2.1 toolbox
(STRIKE-GOLDD = STRuctural Identifiability taKen as Extended-Generalized Observability with Lie Derivatives and Decomposition)
--------------------------------------------------------------------------------------------------------------------------------
Created:       07/12/2015
Last modified: 08/01/2019
Contact:       Alejandro F. Villaverde (afvillaverde@iim.csic.es; ale.fer.vi@gmail.com)

--------------------------
File list:
--------------------------
- install.m                    add folders to the path
- options.m                    user-defined options and default values
- STRIKE_GOLDD.m               main file

- functions/build_OI_ext.m     Function that builds the generalized observability-identifiability matrix with user-defined number of Lie derivatives
- functions/combin_optim.m     Function that performs combinatorial optimization using the Variable Neighbourhood Search metaheuristic (VNS)
- functions/combos.m           Function that tries to find identifiable combinations of otherwise unidentifiable parameters
- functions/decomp.m           Function that decomposes the model into submodels (either defined by the user, or found via optimization) and analyses them 
- functions/elim_and_recalc.m  Function that determines identifiability of individual parameters one by one,
                               by successive elimination of its column in the identifiability matrix and recalculation of its rank
- functions/graph_model.m      Create a graph showing the relations among model states, outputs, and parameters  
- functions/objective_fun.m    Function that calculates the objective function value in the optimization
                               (as the ratio between number of model outputs and parameters, plus a penalty on the number of states)
