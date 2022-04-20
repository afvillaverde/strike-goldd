import sympy as sym

# 2-compartment linear model analysed in:
# "Input-Dependent Structural Identifiability of Nonlinear Systems"
# AF Villaverde, ND Evans, MJ Chappell, JR Banga
# IEEE Control Systems Letters 3 (2), 272-277
#--------------------------------------------------------------------------
# Modified so that 'b' is assumed to be known

x1, x2, x3, x4, x5, x6 = sym.symbols('x1 x2 x3 x4 x5 x6')
u1, u2 = sym.symbols('u1 u2')

# 2 states
x = [[x1], [x2]]

# 1 output
h = [x1]

# 0 known inputs
u = []

# 1 unknown input
w = [[u1], [u2]]

# 4 unknown parameters
p = [[x3], [x4], [x5]]  # x6 -> 'b' is a known constant

# dynamic equations
f = [[-(x3+x4)*x1+x5*x2+x6*u1], [x4*x1-x5*x2]]


variables_locales = locals().copy()