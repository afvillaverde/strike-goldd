import sympy as sym

# Hormonal circuit model with integral feedback
# Originally published in: Karin et al, Mol Syst Biol 2016
# Corresponds to the model in Fig. 1A of: Villaverde & Banga, arXiv:1701.02562

# 2 states
x1 = sym.Symbol('x1')
x2 = sym.Symbol('x2')
x10 = sym.Symbol('x10')
x20 = sym.Symbol('x20')
x = [[x1], [x2]]

# 1 output
h = [x1]

# 1 known input
uu = sym.Symbol('uu')
u0 = sym.Symbol('u0')
u = [uu]

# 2 unknown parameters
p1 = sym.Symbol('p1')
p2 = sym.Symbol('p2')
p = [[p1], [p2]]

# dynamic equation
f = [[u0+uu-p2*x1-p1*x2], [x1-x10]]


variables_locales = locals().copy()
