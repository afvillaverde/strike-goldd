import sympy as sym

# "BIG" (betaIG) model by Topp et al, J Theor Biol 2000
# Originally published in: Karin et al, Mol Syst Biol 2016
# Corresponds to the model in Fig. 1D of: Villaverde & Banga, arXiv:1701.02562

# 3 states
x1 = sym.Symbol('x1')
x2 = sym.Symbol('x2')
x3 = sym.Symbol('x3')
x = [[x1], [x2], [x3]]

# 1 output
h = [x1]

# 1 known input
uu = sym.Symbol('uu')
u0 = sym.Symbol('u0')
u = [uu]

# 5 unknown parameters
p1 = sym.Symbol('p1')
p2 = sym.Symbol('p2')
p3 = sym.Symbol('p3')
p4 = sym.Symbol('p4')
p5 = sym.Symbol('p5')
p = [[p1], [p2], [p3], [p4], [p5]]

# known constants
muplus  = 0.021/(24*60) # turnover of functional mass
muminus = 0.025/(24*60)

# auxiliary functions
rhoG        = x1**2/(p5**2+x1**2)
lambdaplus  = muplus/(1+(8.4/x1)**1.7)
lambdaminus = muminus/(1+(x1/4.8)**8.5)

# dynamic equations
f = [[u0+uu-(p4+p2*x3)*x1], [x2*(lambdaplus-lambdaminus)], [p1*x2*rhoG-p3*x3]]


variables_locales = locals().copy()