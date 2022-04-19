import sympy as sym

k1, k2, k3, k4, k5, k6 = sym.symbols('k1 k2 k3 k4 k5 k6')
ps1, ps2, ps3 = sym.symbols('ps1 ps2 ps3')
s1t, s2t, s3t = sym.symbols('s1t s2t s3t')
KK1, KK2 = sym.symbols('n1 n2')
alpha = sym.Symbol('alpha')
n1, n2 = sym.symbols('n1 n2')

# 3 states
x = [[ps1], [ps2], [ps3]]

# 3 outputs
h = [[ps1], [ps2], [ps3]]

# no known or unknown inputs
u = []
w = []

# 14 unknown parameters
p = [[k1], [k2], [k3], [k4], [k5], [k6], [s1t], [s2t], [s3t], [KK1], [KK2], [n1], [n2], [alpha]]

# dynamic equations
f = [[k1*(s1t-ps1)*(KK1**n1)/(KK1**n1+ps3**n1)-k2*ps1], [k3*(s2t-ps2)*ps1*(1+(alpha*ps3**n2)/(KK2**n2+ps3**n2))-k4*ps2], [k5*(s3t-ps3)*ps2-k6*ps3]]


variables_locales = locals().copy()