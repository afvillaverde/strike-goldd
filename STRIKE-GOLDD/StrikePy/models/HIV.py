import sympy as sym

#--------------------------------------------------------------------------
# HIV Model with Constant and Time Varying Parameters
# Miao H, Xia X, Perelson AS, Wu H.
# "On identifiability of nonlinear ODE models and applications in viral dynamics."
# SIAM review 53.1 (2011): 3-39.
#--------------------------------------------------------------------------

# 3 states
Tu = sym.Symbol('Tu')
Ti = sym.Symbol('Ti')
V = sym.Symbol('V')
x = [[Tu], [Ti], [V]]

# 2 outputs
h = [[V], [Tu+Ti]]

# 0 known inputs
u = []

# 1 unknown input
eta = sym.Symbol('eta')
w = [eta]

# 5 unknown parameters
lambdA = sym.Symbol('lambdA')
rho = sym.Symbol('rho')
N = sym.Symbol('N')
delta = sym.Symbol('delta')
c = sym.Symbol('c')
p = [[lambdA], [rho], [N], [delta], [c]]

# dynamic equations
f = [[lambdA-rho*Tu-eta*Tu*V], [eta*Tu*V-delta*Ti], [N*delta*Ti-c*V]]


variables_locales = locals().copy()