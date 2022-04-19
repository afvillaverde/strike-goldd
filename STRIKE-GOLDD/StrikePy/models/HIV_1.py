import sympy as sym

#--------------------------------------------------------------------------
# HIV Model with Constant and Time Varying Parameters
# Miao H, Xia X, Perelson AS, Wu H.
# "On identifiability of nonlinear ODE models and applications in viral dynamics."
# SIAM review 53.1 (2011): 3-39.
#-------------------------------------------------------------------------

# 4 states
x1 = sym.Symbol('x1')
x2 = sym.Symbol('x2')
x3 = sym.Symbol('x3')
x4 = sym.Symbol('x4')
x = [[x1], [x2], [x3], [x4]]

# 2 outputs
h = [[x1], [x4]]

# 0 known inputs
u = []

# 0 unknown input
w = []

# 10 unknown parameters
b = sym.Symbol('b')
c = sym.Symbol('c')
d = sym.Symbol('d')
q1 = sym.Symbol('q1')
q2 = sym.Symbol('q2')
k1 = sym.Symbol('k1')
k2 = sym.Symbol('k2')
w1 = sym.Symbol('w1')
w2 = sym.Symbol('w2')
s = sym.Symbol('s')
p = [[b], [c], [d], [q1], [q2], [k1], [k2], [w1], [w2], [s]]

# dynamic equations
f = [[-b*x1*x4-d*x1 + s], [b*q1*x1*x4-k1*x2-w1*x2], [b*q2*x1*x4+k1*x2-w2*x3], [-c*x4+k2*x3]]


variables_locales = locals().copy()