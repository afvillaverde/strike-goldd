#--------------------------------------------------------------------------
# File that defines the JAK-STAT model.
# It stores the variables in a mat-file called JAKSTAT.mat.
# The model is taken from:
#--------------------------------------------------------------------------
# Raia V et al. (2011) Dynamic mathematical modeling of IL13-induced
# signaling in Hodgkin and primary mediastinal B-cell lymphoma allows
# prediction of therapeutic targets.
# Cancer Research 71(3):693{704.
#--------------------------------------------------------------------------

import sympy as sym

x1, x2, x3, x4, x5, x6, x8, x10, x11, x13, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, u1 = sym.symbols('x1 x2 x3 x4 x5 x6 x8 x10 x11 x13 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12 t13 t14 t15 t16 t17 t18 t19 t20 t21 t22 t23 u1')

# 10 states:
x = [[x1], [x2], [x3], [x4], [x5], [x6], [x8], [x10], [x11], [x13]]

# 23 parameters (t23 is an initial condition):
p = [[t1], [t2], [t3],[t4],[t5],[t6],[t7],[t8],[t9],[t10],[t11],[t12], [t13],[t14],[t15],[t16],[t17],[t18],[t19],[t20],[t21],[t22]]

# 8 outputs:
h = [[x1 + x3 + x4], [t18*(x3 + x4 + x5 +(0.34-x11))],[t19*(x4 + x5)],[t20*(-x6+2.8)],[t21*x10],[t22*x10*t17/t11],[x13],[-x8+165]]

# 2 known constants:
c1 = 2.265
c2 = 91

# one input:
u = [u1]

# dynamic equations:
f = [[-t1*x1*c1*u1-t5*x1+t6*x2],[t5*x1-t6*x2],[t1*c1*u1*x1-t2*x3*(-x6+2.8)],[t2*x3*(-x6+2.8)-t3*x4],[t3*x4-t4*x5],
    [-t7*x3*x6/(1+t13*x1)-t7*x4*x6/(1+t13*x13)+t8*(-x6+2.8)*c2],[-t9*x8*(-x6+2.8)+t10*(-x8+165)*c2],
    [t11*(-x8+165)],[-t12*c1*u1*x11],[x10*t14/(t15+x10)-t16*x13]]

variables_locales = locals().copy()