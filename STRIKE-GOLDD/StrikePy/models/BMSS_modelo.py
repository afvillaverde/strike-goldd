import sympy as sym

# Hormonal circuit model with integral feedback
# Originally published in: Karin et al, Mol Syst Biol 2016
# Corresponds to the model in Fig. 1A of: Villaverde & Banga, arXiv:1701.02562

# 2 states
x_1 = sym.Symbol('x_1')
s_1 = sym.Symbol('s_1')
h_1 = sym.Symbol('h_1')

x = [[1], [1], [1]]

# 1 output
h = [[1], [1], [1]]

# 1 known input
Y_1 = sym.Symbol('Y_1')
synh_1 = sym.Symbol('synh_1')
Kind_1 = sym.Symbol('Kind_1')
u = []

# 2 unknown parameters
mu_max_1 = sym.Symbol('mu_max_1')
Ks_1 = sym.Symbol('Ks_1')
Ind_1 = 6
p = [[mu_max_1], [Ks_1], [Y_1], [synh_1], [Kind_1]]

# dynamic equation
f = [[mu_max_1*s_1*x_1/(Ks_1 + s_1)], [-mu_max_1*s_1*x_1/(Y_1*(Ks_1 + s_1))], [Ind_1*synh_1/(Ind_1 + Kind_1) - h_1*mu_max_1*s_1/(Ks_1 + s_1)]]


variables_locales = locals().copy()
