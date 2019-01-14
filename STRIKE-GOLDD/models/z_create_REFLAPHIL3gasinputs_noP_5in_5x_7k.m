%--------------------------------------------------------------------------
% The dissolved gas concentrations are inputs (u3/Co2 and u4/O2).
% The file is adapted from a model provided by Flavia Neddermeyer (TU Berlin), 
% Modellbildungsarbeit 2017 mit Philipp Kunz 
%--------------------------------------------------------------------------
clear all;

% 5 states:
syms x1 x2 x3 x4 x5
x = [x1; x2; x3; x4; x5];

% 8 parameters:
syms p1 p2 p3 p4 p5 p6 p7 p8
p = [p1 p2 p3 p4 p5 p6 p7 p8].';

% 5 outputs:
h = x;

% known constants:
syms k1 k2 k3 k4 k5 k6 k7

% u1 (cN*rN); u2 (cP * rP); u3 (pCO2), u4 (pO2); u5 (Korrekturfluide+rN+rP):
syms u1 u2 u3 u4 u5;
u = [u1 u2 u3 u4 u5].';

% dynamic equations:
f = [ 
k1*(p4-1.*p4*u3)*(1.+2./p5)*u3/(1.+k2*u3+k3*u3^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4/(1.+k4*u4+(k4*u4/p6)^p7)*x2/x5/(x2/x5+p8)*x1;
-k1*p3*p4*(1.+2./p5)*u3/(1.+k2*u3+k3*u3^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4/(1.+k4*u4+(k4*u4/p6)^p7)*x2/x5/(x2/x5+p8)*x1+u1;
 -k1*p2*p4*(1.+2./p5)*u3/(1.+k2*u3+k3*u3^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4/(1.+k4*u4+(k4*u4/p6)^p7)*x2/x5/(x2/x5+p8)*x1+u2;
k5*p3*p4*(1.+2./p5)*u3/(1.+k2*u3+k3*u3^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4/(1.+k4*u4+(k4*u4/p6)^p7)*x2/x5/(x2/x5+p8)*x1+u3*x5*p1;
 u5+k7*p4*(1.+2./p5)*u3/(1.+k2*u3+k3*u3^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4/(1.+k4*u4+(k4*u4/p6)^p7)*x2/x5/(x2/x5+p8)*x1-k6;
];

% initial conditions:
ics  = [1.05, 8.59, 22.93, 0, 10];   

% which initial conditions are known:
known_ics = [1,1,1,1,1]; 

save('REFLAPHIL3gasinputs_noP_5in_5x_7k','x','p','u','h','f','ics','known_ics');
