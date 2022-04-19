import sympy as sym

# - "BIG" (betaIG) model by: Topp et al, J Theor Biol 2000
# - Originally published in: Karin et al, Mol Syst Biol 2016
# - The two model versions defined here are analysed in:
# Massonis, Banga, Villaverde. "Automatic reformulation method...", 2020.

# 3 states
G = sym.Symbol('G')
beta = sym.Symbol('beta')
I = sym.Symbol('I')
x = [[G], [beta], [I]]

# 1 output
h = [G]

# 1 known input
inputs = sym.Symbol('inputs')
u = [inputs]

# 5 unknown parameters
p1 = sym.Symbol('p1')
si = sym.Symbol('si')
gamma = sym.Symbol('gamma')
c = sym.Symbol('c')
alpha = sym.Symbol('alpha')
p = [[p1], [si], [gamma], [c], [alpha]]

# known constants
muplus  = 0.021/(24*60) # turnover of functional mass
muminus = 0.025/(24*60)

# auxiliary functions
rhoG        = G**2/(alpha**2+G**2)
lambdaplus  = muplus/(1+(8.4/G)**1.7)
lambdaminus = muminus/(1+(G/4.8)**8.5)

# dynamic equations
f = [[inputs-(c+si*I)*G], [beta*(lambdaplus-lambdaminus)], [p1*beta*rhoG-gamma*I]]


variables_locales = locals().copy()