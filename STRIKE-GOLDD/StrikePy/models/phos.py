import sympy as sym

# 6 states
k1 = sym.Symbol('k1')
k2 = sym.Symbol('k2')
k3 = sym.Symbol('k3')
k4 = sym.Symbol('k4')
k5 = sym.Symbol('k5')
k6 = sym.Symbol('k6')
x1 = sym.Symbol('x1')
x2 = sym.Symbol('x2')
x3 = sym.Symbol('x3')
x4 = sym.Symbol('x4')
x5 = sym.Symbol('x5')
x6 = sym.Symbol('x6')
x = [[x1], [x2], [x3], [x4], [x5], [x6]]

# 2 outputs
h = [[x2], [x3]]

# no known inputs
u = []

# no unknown inputs
w = []

# 5 unknown parameters
p = [[k1], [k2], [k3], [k4], [k5], [k6]]

# dynamic equations
f = [[k4*x6 + k2*x4 - k1*x1*x2], [k3*x4 + k2*x4 + k1*x1*x2], [k5*x6 + k3*x4 - k6*x5*x3], [-k3*x4 - k2*x4 + k1*x1*x2], [k5*x6 + k4*x6 - k6*x5*x3], [-k5*x6 - k4*x6 + k6*x5*x3]]


variables_locales = locals().copy()