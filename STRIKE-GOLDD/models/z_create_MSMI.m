% Species-metabolite interaction models  with Monod growth kinetics 
% analized in:
% "Structural identifiability and observability of microbial community"
% S Díaz-Seoane, E Sellán, AF Villaverdeclear;

%%%%%%%%%%%%%%%%%%%%%% Two species and one metabolite %%%%%%%%%%%%%%%%%%%%%
clear;
syms x1 x2 m1 ...
v11 v21...
k11 k21...
d1 d2...
vas11 vas21...
kas11 kas21...
f1...
das1...
phi111 phi121...
% 3 states
x = [x1;x2;m1];
% inputs
u = [];
% parameters
p = [v11; v21; k11; k21; d1; d2; vas11; vas21; kas11; kas21; f1; das1; phi111; phi121];
% dynamic equations
f = [x1*(v11*m1/(k11+m1*x1)-d1);
x2*(v21*m1/(k21+m1*x2)-d2);
-m1*(vas11*x1/(kas11+m1*x1)+vas21*x2/(kas21+m1*x2)-das1)+f1+phi111*x1*m1/(kas11+x1*m1)+phi121*x2*m1/(kas21+x2*m1)];
% initial conditions
ics = [];
% known initial conditions
known_ics = [0,0];
% outputs
h = [x1;x2;m1];
save('MMSMI_2x1m_hx1x2m1','x','p','h','f','u','ics','known_ics');
h = [x1;x2];
save('MMSMI_2x1m_hx1x2','x','p','h','f','u','ics','known_ics');
h = [x1;m1];
save('MMSMI_2x1m_hx1m1','x','p','h','f','u','ics','known_ics');
h = [x2;m1];
save('MMSMI_2x1m_hx2m1','x','p','h','f','u','ics','known_ics');
h = [x1];
save('MMSMI_2x1m_hx1','x','p','h','f','u','ics','known_ics');
h = [x2];
save('MMSMI_2x1m_hx2','x','p','h','f','u','ics','known_ics');
h = [m1];
save('MMSMI_2x1m_hm1','x','p','h','f','u','ics','known_ics');

%%%%%%%%%%%%%%%%%%%%% Two species and two metabolites %%%%%%%%%%%%%%%%%%%%%
clear;
syms x1 x2 m1 m2 ...
v11 v21 v12 v22...
k11 k21 k12 k22...
d1 d2...
vas11 vas21 vas12 vas22...
kas11 kas21 kas12 kas22...
f1 f2...
das1 das2...
phi111 phi121 phi112 phi122 phi211 phi221 phi212 phi222
% 3 states
x = [x1;x2;m1;m2];
% inputs
u = [];
% parameters
p = [v11; v21; v12; v22; k11; k21; k12; k22; d1; d2; vas11; vas21; vas12; vas22;...
kas11; kas21; kas12; kas22; f1; f2; das1; das2;...
phi111; phi121; phi112; phi122; phi211; phi221; phi212; phi222];
% dynamic equations
f = [x1*(v11*m1/(k11+m1*x1)+v12*m2/(k12+m2*x1)-d1);
x2*(v21*m1/(k21+m1*x2)+v22*m2/(k22+m2*x2)-d2);
-m1*(vas11*x1/(kas11+m1*x1)+vas21*x2/(kas21+m1*x2)+das1)+f1+phi111*x1*m1/(kas11+x1*m1)+phi121*x2*m1/(kas21+x2*m1)+phi112*x1*m2/(kas12+x1*m2)+phi122*x2*m2/(kas22+x2*m2);
-m2*(vas12*x1/(kas12+m2*x1)+vas22*x2/(kas22+m2*x2)+das2)+f2+phi211*x1*m1/(kas11+x1*m1)+phi221*x2*m1/(kas21+x2*m1)+phi212*x1*m2/(kas12+x1*m2)+phi222*x2*m2/(kas22+x2*m2)];
% initial conditions
ics = [];
% known initial conditions
known_ics = [0,0];
% outputs
h = [x1;x2;m1;m2];
save('MMSMI_2x2m_hx1x2m1m2','x','p','h','f','u','ics','known_ics');
h = [x1;x2;m1];
save('MMSMI_2x2m_hx1x2m1','x','p','h','f','u','ics','known_ics');
h = [x1;x2;m2];
save('MMSMI_2x2m_hx1x2m2','x','p','h','f','u','ics','known_ics');
h = [x1;m1;m2];
save('MMSMI_2x2m_hx1m1m2','x','p','h','f','u','ics','known_ics');
h = [x2;m1;m2];
save('MMSMI_2x2m_hx2m1m2','x','p','h','f','u','ics','known_ics');
h = [x1;x2];
save('MMSMI_2x2m_hx1x2','x','p','h','f','u','ics','known_ics');
h = [x1;m1];
save('MMSMI_2x2m_hx1m1','x','p','h','f','u','ics','known_ics');
h = [x1;m2];
save('MMSMI_2x2m_hx1m2','x','p','h','f','u','ics','known_ics');
h = [x2;m1];
save('MMSMI_2x2m_hx2m1','x','p','h','f','u','ics','known_ics');
h = [x2;m2];
save('MMSMI_2x2m_hx2m2','x','p','h','f','u','ics','known_ics');
h = [m1;m2];
save('MMSMI_2x2m_hm1m2','x','p','h','f','u','ics','known_ics');
h = [x1];
save('MMSMI_2x2m_hx1','x','p','h','f','u','ics','known_ics');
h = [x2];
save('MMSMI_2x2m_hx2','x','p','h','f','u','ics','known_ics');
h = [m1];
save('MMSMI_2x2m_hm1','x','p','h','f','u','ics','known_ics');
h = [m2];
save('MMSMI_2x2m_hm2','x','p','h','f','u','ics','known_ics');

%%%%%%%%%%%%%%%%%%%% Three species and two metabolites %%%%%%%%%%%%%%%%%%%%
clear;
syms x1 x2 x3 m1 m2 ...
v11 v21 v12 v22 v31 v32...
k11 k21 k12 k22 k31 k32...
d1 d2 d3...
vas11 vas21 vas12 vas22 vas31 vas32...
kas11 kas21 kas12 kas22 kas31 kas32...
f1 f2...
das1 das2...
phi111 phi121 phi112 phi122 phi211 phi221 phi212 phi222...
phi131 phi132  phi231 phi232
% 3 states
x = [x1;x2;x3;m1;m2];
% inputs
u = [];
% parameters
p = [v11; v21; v12; v22; v31; v32; k11; k21; k12; k22; k31; k32;...
    d1; d2; d3; vas11; vas21; vas12; vas22; vas31; vas32;...
kas11; kas21; kas12; kas22; kas31; kas32; f1; f2; das1; das2;...
phi111; phi121; phi112; phi122; phi211; phi221; phi212; phi222;...
phi131; phi132;  phi231; phi232];
% dynamic equations
f = [x1*(v11*m1/(k11+m1*x1)+v12*m2/(k12+m2*x1)-d1);
x2*(v21*m1/(k21+m1*x2)+v22*m2/(k22+m2*x2)-d2);
x3*(v31*m1/(k31+m1*x3)+v32*m2/(k32+m2*x3)-d3);
-m1*(vas11*x1/(kas11+m1*x1)+vas21*x2/(kas21+m1*x2)+vas31*x3/(kas31+m1*x3)+das1)+f1+phi111*x1*m1/(kas11+x1*m1)+phi121*x2*m1/(kas21+x2*m1)+phi112*x1*m2/(kas12+x1*m2)+phi122*x2*m2/(kas22+x2*m2)+phi131*x3*m1/(kas31+x3*m1)+phi132*x3*m2/(kas32+x3*m2);
-m2*(vas12*x1/(kas12+m2*x1)+vas22*x2/(kas22+m2*x2)+vas32*x3/(kas32+m2*x3)+das2)+f2+phi211*x1*m1/(kas11+x1*m1)+phi221*x2*m1/(kas21+x2*m1)+phi212*x1*m2/(kas12+x1*m2)+phi222*x2*m2/(kas22+x2*m2)+phi231*x3*m1/(kas31+x3*m1)+phi232*x3*m2/(kas32+x3*m2)];
% initial conditions
ics = [];
% known initial conditions
known_ics = [0,0];
% outputs
h = [x1;x2;x3;m1;m2];
save('MMSMI_3x2m_hx1x2x3m1m2','x','p','h','f','u','ics','known_ics');
h = [x1;x2;x3;m1];
save('MMSMI_3x2m_hx1x2x3m1','x','p','h','f','u','ics','known_ics');
h = [x1;x2;x3;m2];
save('MMSMI_3x2m_hx1x2x3m2','x','p','h','f','u','ics','known_ics');
h = [x1;x2;m1;m2];
save('MMSMI_3x2m_hx1x2m1m2','x','p','h','f','u','ics','known_ics');
h = [x1;x3;m1;m2];
save('MMSMI_3x2m_hx1x3m1m2','x','p','h','f','u','ics','known_ics');
h = [x2;x3;m1;m2];
save('MMSMI_3x2m_hx2x3m1m2','x','p','h','f','u','ics','known_ics');
h = [x1;x2;x3];
save('MMSMI_3x2m_hx1x2x3','x','p','h','f','u','ics','known_ics');
h = [x1;x2;m1];
save('MMSMI_3x2m_hx1x2m1','x','p','h','f','u','ics','known_ics');
h = [x1;x2;m2];
save('MMSMI_3x2m_hx1x2m2','x','p','h','f','u','ics','known_ics');
h = [x1;x3;m1];
save('MMSMI_3x2m_hx1x3m1','x','p','h','f','u','ics','known_ics');
h = [x1;x3;m2];
save('MMSMI_3x2m_hx1x3m2','x','p','h','f','u','ics','known_ics');
h = [x1;m1;m2];
save('MMSMI_3x2m_hx1m1m2','x','p','h','f','u','ics','known_ics');
h = [x2;x3;m1];
save('MMSMI_3x2m_hx2x3m1','x','p','h','f','u','ics','known_ics');
h = [x2;x3;m2];
save('MMSMI_3x2m_hx2x3m2','x','p','h','f','u','ics','known_ics');
h = [x2;m1;m2];
save('MMSMI_3x2m_hx2m1m2','x','p','h','f','u','ics','known_ics');
h = [x3;m1;m2];
save('MMSMI_3x2m_hx3m1m2','x','p','h','f','u','ics','known_ics');
h = [x1;x2];
save('MMSMI_3x2m_hx1x2','x','p','h','f','u','ics','known_ics');
h = [x1;x3];
save('MMSMI_3x2m_hx1x3','x','p','h','f','u','ics','known_ics');
h = [x1;m1];
save('MMSMI_3x2m_hx1m1','x','p','h','f','u','ics','known_ics');
h = [x1;m2];
save('MMSMI_3x2m_hx1m2','x','p','h','f','u','ics','known_ics');
h = [x2;x3];
save('MMSMI_3x2m_hx2x3','x','p','h','f','u','ics','known_ics');
h = [x2;m1];
save('MMSMI_3x2m_hx2m1','x','p','h','f','u','ics','known_ics');
h = [x2;m2];
save('MMSMI_3x2m_hx2m2','x','p','h','f','u','ics','known_ics');
h = [x3;m1];
save('MMSMI_3x2m_hx3m1','x','p','h','f','u','ics','known_ics');
h = [x3;m2];
save('MMSMI_3x2m_hx3m2','x','p','h','f','u','ics','known_ics');
h = [m1;m2];
save('MMSMI_3x2m_hm1m2','x','p','h','f','u','ics','known_ics');
h = [x1];
save('MMSMI_3x2m_hx1','x','p','h','f','u','ics','known_ics');
h = [x2];
save('MMSMI_3x2m_hx2','x','p','h','f','u','ics','known_ics');
h = [x3];
save('MMSMI_3x2m_hx3','x','p','h','f','u','ics','known_ics');
h = [m1];
save('MMSMI_3x2m_hm1','x','p','h','f','u','ics','known_ics');
h = [m2];
save('MMSMI_3x2m_hm2','x','p','h','f','u','ics','known_ics');