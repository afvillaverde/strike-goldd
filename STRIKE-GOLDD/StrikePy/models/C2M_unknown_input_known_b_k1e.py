import sympy as sym

# 2-compartment linear model analysed in:
# "Input-Dependent Structural Identifiability of Nonlinear Systems"
# AF Villaverde, ND Evans, MJ Chappell, JR Banga
# IEEE Control Systems Letters 3 (2), 272-277
#--------------------------------------------------------------------------
# Modified so that 'b' and 'kie' are assumed to be known

x1, x2, x4, x5 = sym.symbols('x1 x2 x4 x5')
u1 = sym.symbols('u1')

# 2 states
x = [[x1], [x2]]

# 1 output
h = [x1]

# 0 known inputs
u = []

# 1 unknown input
w = [u1]

# 2 unknown parameters (x4 = k12, x5 = k21)
p = [[x4], [x5]]


# 2 known constants (x3 = k1e, x6 = b)
x3, x6 = sym.symbols('x3 x6')

# dynamic equations
f = [[-(x3+x4)*x1+x5*x2+x6*u1], [x4*x1-x5*x2]]


variables_locales = locals().copy()