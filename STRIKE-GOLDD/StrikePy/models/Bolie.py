import sympy as sym

# Model from J.W. Bolie. "Coefficients of normal blood glucose regulation".
# J. Appl. Physiol., 16(5):783-788, 1961.

# 2 states
q1 = sym.Symbol('q1')
q2 = sym.Symbol('q2')
x = [[q1], [q2]]

# 1 output
Vp = sym.Symbol('Vp')
h = [q1/Vp]

# 1 known input
delta = sym.Symbol('delta')
u = [delta]

# 5 unknown parameters
p1 = sym.Symbol('p1')
p2 = sym.Symbol('p2')
p3 = sym.Symbol('p3')
p4 = sym.Symbol('p4')
p = [[p1], [p2], [p3], [p4], [Vp]]

# dynamic equations
f = [[p1*q1-p2*q2+delta], [p3*q2+p4*q1]]


variables_locales = locals().copy()

