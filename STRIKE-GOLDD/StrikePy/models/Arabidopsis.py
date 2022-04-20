import sympy as sym

#--------------------------------------------------------------------------
# File that creates the Arabidopsis thaliana model.
# It stores it in a mat-file named Arabidopsis.mat.
# The model is taken from:
#--------------------------------------------------------------------------
# Locke J, Millar A, Turner M (2005) Modelling genetic networks with noisy
# and varied experimental data: the circadian clock in arabidopsis thaliana.
# Journal of Theoretical Biology 234: 383–393.
#--------------------------------------------------------------------------

# 7 states
x1 = sym.Symbol('x1')
x2 = sym.Symbol('x2')
x3 = sym.Symbol('x3')
x4 = sym.Symbol('x4')
x5 = sym.Symbol('x5')
x6 = sym.Symbol('x6')
x7 = sym.Symbol('x7')
x = [[x1], [x2], [x3], [x4], [x5], [x6], [x7]]

# 2 outputs
h = [[x1], [x4]]

# 1 known input
u1 = sym.Symbol('u1')
u = [u1]

# unknown parameters
a = sym.Symbol('a')
g1 = sym.Symbol('g1')
g2 = sym.Symbol('g2')
k1 = sym.Symbol('k1')
k2 = sym.Symbol('k2')
k3 = sym.Symbol('k3')
k4 = sym.Symbol('k4')
k5 = sym.Symbol('k5')
k6 = sym.Symbol('k6')
k7 = sym.Symbol('k7')
m1 = sym.Symbol('m1')
m2 = sym.Symbol('m2')
m3 = sym.Symbol('m3')
m4 = sym.Symbol('m4')
m5 = sym.Symbol('m5')
m6 = sym.Symbol('m6')
m7 = sym.Symbol('m7')
n1 = sym.Symbol('n1')
n2 = sym.Symbol('n2')
p1 = sym.Symbol('p1')
p2 = sym.Symbol('p2')
p3 = sym.Symbol('p3')
q1 = sym.Symbol('q1')
q2 = sym.Symbol('q2')
r1 = sym.Symbol('r1')
r2 = sym.Symbol('r2')
r3 = sym.Symbol('r3')
r4 = sym.Symbol('r4')

p = [[a], [n1], [r3], [g1], [g2], [k1], [k2], [k3], [k4], [k5], [k6], [k7], [m1], [m2], [m3], [m4], [m5], [m6], [m7], [n2], [p1], [p2], [p3], [q1], [q2], [r1], [r2], [r4]]

# dynamic equation
f = [[n1*x6**a/(g1**a+x6**a)-m1*x1/(k1+x1)+q1*x7*u1], [p1*x1-r1*x2+r2*x3-m2*x2/(k2+x2)], [r1*x2-r2*x3-m3*x3/(k3+x3)], [n2*g2**2/(g2**2+x3**2)-m4*x4/(k4+x4)], [p2*x4-r3*x5+r4*x6-m5*x5/(k5+x5)], [r3*x5-r4*x6-m6*x6/(k6+x6)], [p3-m7*x7/(k7+x7)-(p3+q2*x7)*u1]]


variables_locales = locals().copy()
