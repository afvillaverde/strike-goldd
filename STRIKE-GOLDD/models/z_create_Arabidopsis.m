%--------------------------------------------------------------------------
% File that creates the Arabidopsis thaliana model. 
% It stores it in a mat-file named Arabidopsis.mat.
% The model is taken from:
%--------------------------------------------------------------------------
% Locke J, Millar A, Turner M (2005) Modelling genetic networks with noisy
% and varied experimental data: the circadian clock in arabidopsis thaliana. 
% Journal of Theoretical Biology 234: 383–393.
%--------------------------------------------------------------------------
clear all;

syms x1 x2 x3 x4 x5 x6 x7...
     a...
     g1 g2 k1 k2 k3 k4 k5 k6 k7 m1 m2 m3 m4 m5 m6 m7 n1 n2...
     p1 p2 p3 q1 q2 r1 r2 r3 r4...
     u1
 
% states:
x    = [x1 x2 x3 x4 x5 x6 x7].';  

% outputs:
h    = [x1 x4].'; 

% parameters:
p    = [a n1 r3 g1 g2 k1 k2 k3 k4 k5 k6 k7 m1 m2 m3 m4 m5 m6 m7 n2...
     p1 p2 p3 q1 q2 r1 r2 r4].'; 
 
% inputs:
u    = u1;
 
% dynamic equations: 
f    = [ n1*x6^a/(g1^a+x6^a)-m1*x1/(k1+x1)+q1*x7*u1;        
         p1*x1-r1*x2+r2*x3-m2*x2/(k2+x2);    
         r1*x2-r2*x3-m3*x3/(k3+x3);
         n2*g2^2/(g2^2+x3^2)-m4*x4/(k4+x4);
         p2*x4-r3*x5+r4*x6-m5*x5/(k5+x5);
         r3*x5-r4*x6-m6*x6/(k6+x6);
         p3-m7*x7/(k7+x7)-(p3+q2*x7)*u1];
     
% initial conditions:
ics  = [0,0,0,0,0,0,0];    

% which initial conditions are known:
known_ics = [1 1 1 1 1 1 1];

save('Arabidopsis','x','h','p','f','u','ics','known_ics'); 