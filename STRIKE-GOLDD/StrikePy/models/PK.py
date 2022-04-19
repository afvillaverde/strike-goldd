import sympy as sym

# 4 states
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
k7 = sym.Symbol('k7')
s2 = sym.Symbol('s2')
s3 = sym.Symbol('s3')
x = [[x1], [x2], [x3], [x4]]

# 2 outputs
h = [[s2*x2], [s3*x3]]

# 1 known inputs
u1 = sym.Symbol('u1')
u = [u1]

# no unknown inputs
w = []

# 10 unknown parameters
p = [[k1], [k2], [k3], [k4], [k5], [k6], [k7], [s2], [s3]]

# dynamic equations
f = [[u1-(k1+k2)*x1], [k1*x1-(k3+k6+k7)*x2+k5*x4], [k2*x1+k3*x2-k4*x3], [k6*x2-k5*x4]]


variables_locales = locals().copy()