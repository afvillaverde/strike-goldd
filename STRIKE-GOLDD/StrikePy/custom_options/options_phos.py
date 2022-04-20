#==============================================================================
#THE USER CAN DEFINE THE PROBLEM AND SET OPTIONS IN THE FOLLOWING LINES:
#==============================================================================
import sympy as sym
from math import inf
###############################################################################
# (1) NAME OF THE MODEL TO BE STUDIED:
modelname = 'phos'
##############################################################################
# (2) FISPO ANALYSIS OPTIONS:
checkObser = 1    # check state observability, i.e. identifiability of initial conditions (1 = yes; 0 = no).
maxLietime = 100  # max. time allowed for calculating 1 Lie derivative (seconds)
nnzDerU = [inf] # Number of non-zero known input derivatives in each experiment (Rows=inputs;Columns=experiments)
nnzDerW    = [inf] # numbers of nonzero derivatives of the unmeasured inputs (w); may be 'inf'
###############################################################################
# (3) KNOWN/IDENTIFIABLE PARAMETERS (parameters assumed known, or already classified as identifiable):
prev_ident_pars = []
# # An example of use would be:
    # x2 = sym.Symbol('x2')
    # x5 = sym.Symbol('x5')
    # prev_ident_pars = [x2, x5]


