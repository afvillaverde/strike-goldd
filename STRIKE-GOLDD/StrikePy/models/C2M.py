import sympy as sym

# 2 states
x1 = sym.Symbol('x1')  # define the symbolic variable x1
x2 = sym.Symbol('x2')  # define the symbolic variable x2
x3 = sym.Symbol('x3')  # define the symbolic variable x3
x4 = sym.Symbol('x4')  # define the symbolic variable x4
x5 = sym.Symbol('x5')  # define the symbolic variable x5
x6 = sym.Symbol('x6')  # define the symbolic variable x6
x = [[x1], [x2]]

# 1 output
h = [x1]

# 1 known input
u1 = sym.Symbol('u1')  # define the symbolic variable u1
u = [u1]

# 0 unknown inputs
w = []

# 4 unknown parameters
p = [[x3], [x4], [x5], [x6]]

# dynamic equations
f = [[-(x3+x4)*x1+x5*x2+x6*u1], [x4*x1-x5*x2]]


variables_locales = locals().copy()
