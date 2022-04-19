import sympy as sym

k3 = sym.Symbol('k3')
k4 = sym.Symbol('k4')
k5 = sym.Symbol('k5')
k6 = sym.Symbol('k6')
k7 = sym.Symbol('k7')
x1 = sym.Symbol('x1')
x2 = sym.Symbol('x2')
x3 = sym.Symbol('x3')
x4 = sym.Symbol('x4')
x5 = sym.Symbol('x5')
a = sym.Symbol('a')
b = sym.Symbol('b')
d = sym.Symbol('d')

# 5 states
x = [[x1], [x2], [x3], [x4], [x5]]

# 1 output
h = [x5]

# 0 known inputs
u = []

# no unknown inputs
w = []

# 5 unknown parameters
p = [[k3], [k4], [k5], [k6], [k7]]

# dynamic equations
f = [[-(k3+k7)*x1+k4*x2], [k3*x1-(k4+a*k5+b*d*k5)*x2+k6*x3+k6*x4+k5*x2*x3+k5*x2*x4], [a*k5*x2-k6*x3-k5*x2*x4], [b*d*k5*x2-k6*x4-k5*x2*x4], [k7*x1]]


variables_locales = locals().copy()