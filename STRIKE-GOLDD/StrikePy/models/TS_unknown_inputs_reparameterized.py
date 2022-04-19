import sympy as sym

# 2 states
x1 = sym.Symbol('x1')
x2 = sym.Symbol('x2')
x6 = sym.Symbol('x6')
x = [[x1], [x2]]

# 2 outputs
h = [[x1], [x2]]

# 0 known inputs
u = []

#1 unknown input
w1= sym.Symbol('w1')
w2= sym.Symbol('w2')
w = [[w1], [w2]]

# 10 unknown parameters
k01 = sym.Symbol('k01')
k1 = sym.Symbol('k1')
ntetr = sym.Symbol('ntetr')
k02 = sym.Symbol('k02')
k2 = sym.Symbol('k2')
nlaci = sym.Symbol('nlaci')
p = [[k01], [k1], [ntetr], [k02], [k2], [nlaci]]

# dynamic equations
f = [
	[k01 + k1/(1+(x2/(1+w1))**ntetr) - x1],
    [k02 + k2/(1+(x1/(1+w2))**nlaci) - x2]
]


variables_locales = locals().copy()