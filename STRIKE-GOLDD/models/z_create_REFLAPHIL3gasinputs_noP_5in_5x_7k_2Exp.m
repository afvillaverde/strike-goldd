%--------------------------------------------------------------------------
% The dissolved gas concentrations are inputs (u3/Co2 and u4/O2).
% The file is adapted from a model provided by Flavia Neddermeyer (TU Berlin), 
% Modellbildungsarbeit 2017 mit Philipp Kunz 
%--------------------------------------------------------------------------
clear all;

% 5 --> 10 states:
syms x1Exp1 x2Exp1 x3Exp1 x4Exp1 x5Exp1 x1Exp2 x2Exp2 x3Exp2 x4Exp2 x5Exp2
x = [x1Exp1 x2Exp1 x3Exp1 x4Exp1 x5Exp1 x1Exp2 x2Exp2 x3Exp2 x4Exp2 x5Exp2].';

% 8 parameters:
syms p1 p2 p3 p4 p5 p6 p7 p8
p = [p1 p2 p3 p4 p5 p6 p7 p8].';

% 5 outputs:
h = x;

% known constants:
syms k1 k2 k3 k4 k5 k6 k7

% u1 (cN*rN); u2 (cP * rP); u3 (pCO2), u4 (pO2); u5 (Korrekturfluide+rN+rP):
syms u1Exp1 u2Exp1 u3Exp1 u4Exp1 u5Exp1 u1Exp2 u2Exp2 u3Exp2 u4Exp2 u5Exp2
u = [u1Exp1 u2Exp1 u3Exp1 u4Exp1 u5Exp1 u1Exp2 u2Exp2 u3Exp2 u4Exp2 u5Exp2].';

% dynamic equations:
f = [ 
k1*(p4-1.*p4*u3Exp1)*(1.+2./p5)*u3Exp1/(1.+k2*u3Exp1+k3*u3Exp1^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4Exp1/(1.+k4*u4Exp1+(k4*u4Exp1/p6)^p7)*x2Exp1/x5Exp1/(x2Exp1/x5Exp1+p8)*x1Exp1;
-k1*p3*p4*(1.+2./p5)*u3Exp1/(1.+k2*u3Exp1+k3*u3Exp1^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4Exp1/(1.+k4*u4Exp1+(k4*u4Exp1/p6)^p7)*x2Exp1/x5Exp1/(x2Exp1/x5Exp1+p8)*x1Exp1+u1Exp1;
 -k1*p2*p4*(1.+2./p5)*u3Exp1/(1.+k2*u3Exp1+k3*u3Exp1^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4Exp1/(1.+k4*u4Exp1+(k4*u4Exp1/p6)^p7)*x2Exp1/x5Exp1/(x2Exp1/x5Exp1+p8)*x1Exp1+u2Exp1;
k5*p3*p4*(1.+2./p5)*u3Exp1/(1.+k2*u3Exp1+k3*u3Exp1^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4Exp1/(1.+k4*u4Exp1+(k4*u4Exp1/p6)^p7)*x2Exp1/x5Exp1/(x2Exp1/x5Exp1+p8)*x1Exp1+u3Exp1*x5Exp1*p1;
 u5Exp1+k7*p4*(1.+2./p5)*u3Exp1/(1.+k2*u3Exp1+k3*u3Exp1^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4Exp1/(1.+k4*u4Exp1+(k4*u4Exp1/p6)^p7)*x2Exp1/x5Exp1/(x2Exp1/x5Exp1+p8)*x1Exp1-k6;
k1*(p4-1.*p4*u3Exp2)*(1.+2./p5)*u3Exp2/(1.+k2*u3Exp2+k3*u3Exp2^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4Exp2/(1.+k4*u4Exp2+(k4*u4Exp2/p6)^p7)*x2Exp2/x5Exp2/(x2Exp2/x5Exp2+p8)*x1Exp2;
-k1*p3*p4*(1.+2./p5)*u3Exp2/(1.+k2*u3Exp2+k3*u3Exp2^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4Exp2/(1.+k4*u4Exp2+(k4*u4Exp2/p6)^p7)*x2Exp2/x5Exp2/(x2Exp2/x5Exp2+p8)*x1Exp2+u1Exp2;
 -k1*p2*p4*(1.+2./p5)*u3Exp2/(1.+k2*u3Exp2+k3*u3Exp2^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4Exp2/(1.+k4*u4Exp2+(k4*u4Exp2/p6)^p7)*x2Exp2/x5Exp2/(x2Exp2/x5Exp2+p8)*x1Exp2+u2Exp2;
k5*p3*p4*(1.+2./p5)*u3Exp2/(1.+k2*u3Exp2+k3*u3Exp2^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4Exp2/(1.+k4*u4Exp2+(k4*u4Exp2/p6)^p7)*x2Exp2/x5Exp2/(x2Exp2/x5Exp2+p8)*x1Exp2+u3Exp2*x5Exp2*p1;
 u5Exp2+k7*p4*(1.+2./p5)*u3Exp2/(1.+k2*u3Exp2+k3*u3Exp2^2/p5^2)*(1.+p7/p6*(p7-1.)^((1.-1.*p7)/p7))*u4Exp2/(1.+k4*u4Exp2+(k4*u4Exp2/p6)^p7)*x2Exp2/x5Exp2/(x2Exp2/x5Exp2+p8)*x1Exp2-k6;
 ];

% initial conditions:
ics  = [1.05, 8.59, 22.93, 0, 10,1.05, 8.59, 22.93, 0, 10];   

% which initial conditions are known:
known_ics = [1,1,1,1,1,1,1,1,1,1]; 

save('REFLAPHIL3gasinputs_noP_5in_5x_7k_2Exp','x','p','u','h','f','ics','known_ics');
