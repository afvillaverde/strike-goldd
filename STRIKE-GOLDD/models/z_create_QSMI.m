% Quadratic species-metabolite interaction models analized in:
% "Structural identifiability and observability of microbial community"
% S Díaz-Seoane, E Sellán, AF Villaverde

%%%%%%%%%%%%%%%%%%%%%% Two species and one metabolite %%%%%%%%%%%%%%%%%%%%%
clear;
syms psi11 psi21...
d1 d2...
f1...
das1...
k11 k21...
phi11 phi21...
x1 x2 m1...
% 3 states
x = [x1;x2;m1];
% inputs
u = [];
% parameters
p =[psi11; d1; psi21; d2; f1; das1; k11; k21; phi11; phi21];
% dynamic equations
f = [x1*(psi11*m1-d1);
x2*(psi21*m1-d2);
f1-das1*m1-m1*(k11*x1+k21*x2)+phi11*x1*m1+phi21*x2*m1];
% initial conditions
ics = [];
% known initial conditions
known_ics = [0,0];
% outputs
h = [x1;x2;m1];
save('QSMI_2x1m_hx1x2m1','x','p','h','f','u','ics','known_ics');
h = [x1;x2];
save('QSMI_2x1m_hx1x2','x','p','h','f','u','ics','known_ics');
h = [x1;m1];
save('QSMI_2x1m_hx1m1','x','p','h','f','u','ics','known_ics');
h = [x2;m1];
save('QSMI_2x1m_hx2m1','x','p','h','f','u','ics','known_ics');
h = x1;
save('QSMI_2x1y_hx1','x','p','h','f','u','ics','known_ics');
h = x2;
save('QSMI_2x1y_hx2','x','p','h','f','u','ics','known_ics');
h = m1;
save('QSMI_2x1y_hy1','x','p','h','f','u','ics','known_ics');

%%%%%%%%%%%%%%%%%%%% Three species and one metabolite %%%%%%%%%%%%%%%%%%%%%
clear;
syms psi11 psi21 psi31...
d1 d2  d3...
f1...
das1...
k11 k21 k31...
phi11 phi21 phi31...
x1 x2 x3 m1...
% 4 states
x = [x1;x2;x3;m1];
% inputs
u = [];
% parameters
p =[psi11; d1; psi21; d2; f1; das1; k11; k21; phi11; phi21; psi31; d3; k31; phi31];
% dynamic equations
f = [x1*(psi11*m1-d1);
x2*(psi21*m1-d2);
x3*(psi31*m1-d3);
f1-das1*m1-m1*(k11*x1+k21*x2+k31*x3)+phi11*x1*m1+phi21*x2*m1+phi31*x3*m1];
% initial conditions
ics = [];
% known initial conditions
known_ics = [0,0];
% outputs
h = [x1;x2;x3;m1];
save('QSMI_3x1m_x1x2x3m1','x','p','h','f','u','ics','known_ics');
h = [x1;x2;x3];
save('QSMI_3x1m_hx1x2x3','x','p','h','f','u','ics','known_ics');
h = [x1;m1];
save('QSMI_3x1m_hx1m1','x','p','h','f','u','ics','known_ics');
h = [x2;m1];
save('QSMI_3x1m_hx2m1','x','p','h','f','u','ics','known_ics');
h = [x3;m1];
save('QSMI_3x1m_hx3m1','x','p','h','f','u','ics','known_ics');
h = x1;
save('QSMI_3x1m_hx1','x','p','h','f','u','ics','known_ics');
h = x2;
save('QSMI_3x1m_hx2','x','p','h','f','u','ics','known_ics');
h = x3;
save('QSMI_3x1m_hx3','x','p','h','f','u','ics','known_ics');
h = m1;
save('QSMI_3x1m_hm1','x','p','h','f','u','ics','known_ics');
h = [x1;x2;m1];
save('QSMI_3x1m_hx1x2m1','x','p','h','f','u','ics','known_ics');
h = [x1;x3;m1];
save('QSMI_3x1m_hx1x3m1','x','p','h','f','u','ics','known_ics');
h = [x2;x3;m1];
save('QSMI_3x1m_hx2x3m1','x','p','h','f','u','ics','known_ics');
h = [x1;x2];
save('QSMI_3x1m_hx1x2','x','p','h','f','u','ics','known_ics');
h = [x1;x3];
save('QSMI_3x1m_hx1x3','x','p','h','f','u','ics','known_ics');
h = [x2;x3];
save('QSMI_3x1m_hx2x3','x','p','h','f','u','ics','known_ics');

%%%%%%%%%%%%%%%%%%%%%% Two species and two metabolites %%%%%%%%%%%%%%%%%%%%
clear;
syms psi11 psi21 psi12 psi22...
d1 d2...
f1 f2...
das1 das2...
k11 k21 k12 k22...
phi111 phi211 phi121 phi221 phi112 phi212 phi122 phi222...
x1 x2 m1 m2...
% 4 states
x = [x1;x2;m1;m2];
% inputs
u = [];
% parameters
p =[psi11 psi21 psi12 psi22...
d1 d2...
f1 f2...
das1 das2...
k11 k21 k12 k22...
phi111 phi211 phi121 phi221 phi112 phi212 phi122 phi222];
% dynamic equations
f = [x1*(psi11*m1+psi12*m2-d1);
x2*(psi21*m1+psi22*m2-d2);
f1-das1*m1-m1*(k11*x1+k21*x2)+phi111*x1*m1+phi211*x2*m1+ phi121*x1*m2+phi221*x2*m2;
f2-das2*m2-m2*(k12*x1+k22*x2)+phi112*x1*m1+phi212*x2*m1+phi122*x1*m2+ phi222*x2*m2];
% initial conditions
ics = [];
% known initial conditions
known_ics = [];
% outputs
h = [x1;x2;m1;m2];
save('QSMI_2x2m_hx1x2m1m2','x','p','h','f','u','ics','known_ics');
h = [x1;x2;m1];
save('QSMI_2x2m_hx1x2m1','x','p','h','f','u','ics','known_ics');
h = [x1;x2;m2];
save('QSMI_2x2m_hx1x2m2','x','p','h','f','u','ics','known_ics');
h = [x1;m1;m2];
save('QSMI_2x2m_hx1m1m2','x','p','h','f','u','ics','known_ics');
h = [x2;m1;m2];
save('QSMI_2x2m_hx2m1m2','x','p','h','f','u','ics','known_ics');
h = [x1;x2];
save('QSMI_2x2m_hx1x2','x','p','h','f','u','ics','known_ics');
h = [m1;m2];
save('QSMI_2x2m_hm1m2','x','p','h','f','u','ics','known_ics');
h = [x1;m1];
save('QSMI_2x2m_hx1m1','x','p','h','f','u','ics','known_ics');
h = [x1;m2];
save('QSMI_2x2m_hx1m2','x','p','h','f','u','ics','known_ics');
h = [x2;m1];
save('QSMI_2x2m_hx2m1','x','p','h','f','u','ics','known_ics');
h = [x2;m2];
save('QSMI_2x2m_hx2m2','x','p','h','f','u','ics','known_ics');
h = x1;
save('QSMI_2x2m_hx1','x','p','h','f','u','ics','known_ics');
h = x2;
save('QSMI_2x2m_hx2','x','p','h','f','u','ics','known_ics');
h = m1;
save('QSMI_2x2m_hm1','x','p','h','f','u','ics','known_ics');
h = m2;
save('QSMI_2x2m_hm2','x','p','h','f','u','ics','known_ics');

%%%%%%%%%%%%%%%%%%%%% Three species and two metabolites %%%%%%%%%%%%%%%%%%%
clear;
syms psi11 psi21 psi31 psi12 psi22 psi32...
d1 d2 d3...
f1 f2...
das1 das2...
k11 k21 k31 k12 k22 k32...
phi111 phi211 phi311 phi121 phi221 phi321 phi112 phi212 phi312 phi122 phi222 phi322...
x1 x2 x3 m1 m2

% 5 states
x = [x1;x2;x3;m1;m2];
% inputs
u = [];
% parameters
p =[psi11 psi21 psi31 psi12 psi22 psi32...
d1 d2 d3...
f1 f2...
das1 das2...
k11 k21 k31 k12 k22 k32...
phi111 phi211 phi311 phi121 phi221 phi321 phi112 phi212 phi312 phi122 phi222 phi322];
% dynamic equations
f = [x1*(psi11*m1+psi12*m2-d1);
x2*(psi21*m1+psi22*m2-d2);
x3*(psi31*m1+psi32*m2-d3);
f1-das1*m1-m1*(k11*x1+k21*x2+k31*x3)+ phi111*x1*m1 + phi211*x2*m1 + phi311*x3*m1 + phi121*x1*m2 + phi221*x2*m2 + phi321*x3*m2 ;
f2-das2*m2-m2*(k12*x1+k22*x2+k32*x3)+ phi112*x1*m1 + phi212*x2*m1 + phi312*x3*m1 + phi122*x1*m2 + phi222*x2*m2 + phi322*x3*m2 ];
% initial conditions
ics = [];
% known initial conditions
known_ics = [];
% outputs
h = [x1;x2;x3;m1;m2];
save('QSMI_3x2m_hx1x2x3m1m2','x','p','h','f','u','ics','known_ics');
h = [x1;x2;m1;m2];
save('QSMI_3x2m_hx1x2m1m2','x','p','h','f','u','ics','known_ics');
h = [x1;x3;m1;m2];
save('QSMI_3x2m_hx1x3m1m2','x','p','h','f','u','ics','known_ics');
h = [x2;x3;m1;m2];
save('QSMI_3x2m_hx2x3y1m2','x','p','h','f','u','ics','known_ics');
h = [x1;x2;x3;m1];
save('QSMI_3x2m_hx1x2x3m1','x','p','h','f','u','ics','known_ics');
h = [x1;x2;x3;m2];
save('QSMI_3x2m_hx1x2x3m2','x','p','h','f','u','ics','known_ics');
h = [x1;x2;x3];
save('QSMI_3x2m_hx1x2x3','x','p','h','f','u','ics','known_ics');
h = [x1;x2;m1];
save('QSMI_3x2m_hx1x2m1','x','p','h','f','u','ics','known_ics');
h = [x1;x3;m1];
save('QSMI_3x2m_hx1x3m1','x','p','h','f','u','ics','known_ics');
h = [x2;x3;m1];
save('QSMI_3x2m_hx2x3m1','x','p','h','f','u','ics','known_ics');
h = [x1;x2;m2];
save('QSMI_3x2m_hx1x2m2','x','p','h','f','u','ics','known_ics');
h = [x1;x3;m2];
save('QSMI_3x2m_hx1x3m2','x','p','h','f','u','ics','known_ics');
h = [x2;x3;m2];
save('QSMI_3x2m_hx2x3m2','x','p','h','f','u','ics','known_ics');
h = [x1;m1;m2];
save('QSMI_3x2m_hx1m1m2','x','p','h','f','u','ics','known_ics');
h = [x2;m1;m2];
save('QSMI_3x2m_hx2m1m2','x','p','h','f','u','ics','known_ics');
h = [x3;m1;m2];
save('QSMI_3x2m_hx3m1m2','x','p','h','f','u','ics','known_ics');
h = [m1;m2];
save('QSMI_3x2m_hm1m2','x','p','h','f','u','ics','known_ics');
h = [x1;m1];
save('QSMI_3x2m_hx1m1','x','p','h','f','u','ics','known_ics');
h = [x2;m1];
save('QSMI_3x2m_hx2m1','x','p','h','f','u','ics','known_ics');
h = [x3;m1];
save('QSMI_3x2m_hx3m1','x','p','h','f','u','ics','known_ics');
h = [x1;m2];
save('QSMI_3x2m_hx1m2','x','p','h','f','u','ics','known_ics');
h = [x2;m2];
save('QSMI_3x2m_hx2m2','x','p','h','f','u','ics','known_ics');
h = [x3;m2];
save('QSMI_3x2m_hx3m2','x','p','h','f','u','ics','known_ics');
h = [x1;x2];
save('QSMI_3x2m_hx1x2','x','p','h','f','u','ics','known_ics');
h = [x1;x3];
save('QSMI_3x2m_hx1x3','x','p','h','f','u','ics','known_ics');
h = [x2;x3];
save('QSMI_3x2m_hx2x3','x','p','h','f','u','ics','known_ics');
h = x1;
save('QSMI_3x2m_hx1','x','p','h','f','u','ics','known_ics');
h = x2;
save('QSMI_3x2m_hx2','x','p','h','f','u','ics','known_ics');
h = x3;
save('QSMI_3x2m_hx3','x','p','h','f','u','ics','known_ics');
h = m1;
save('QSMI_3x2m_hm1','x','p','h','f','u','ics','known_ics');
h = m2;
save('QSMI_3x2m_hm2','x','p','h','f','u','ics','known_ics');