import sympy as sym

# 1 state
x1 = sym.Symbol('x1')
k = sym.Symbol('k')
s1 = sym.Symbol('s1')
s2 = sym.Symbol('s2')

x = [x1]

# 1 output
h = [s1*x1/(1+s2*x1)]

# no known inputs


# 3 unknown parameters
p = [[k], [s1], [s2]]

# dynamic equations
f = [-2*x1*x1*k]


variables_locales = locals().copy()