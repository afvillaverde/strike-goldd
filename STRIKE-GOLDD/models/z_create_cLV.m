% Composite Lotka-Volterra models introducer by:
% "Compositional Lotka-Volterra describes microbial dynamics in the simplex"
% TA Joseph, L Shenhav, JB Xavier, E Halperin, and I Pe’er
% PLOS Comput. Biol., vol. 16, no. 5, p. e1007917, 2020
% Analized in:
% "Structural identifiability and observability of microbial community"
% S Díaz-Seoane, E Sellán, AF Villaverde

%%%%%%%%%%%%%%%%%%% Two species no external disturbance %%%%%%%%%%%%%%%%%%%
clear;
syms g1...
A11 A12...
pi1...
% 1 state (pi2 = 1 - pi1)
x = pi1;
% 0 input
u = [];
% parameters
p =[g1; A11; A12];
% dynamic equations
f = pi1*(g1+A11*pi1+A12*(1-pi1)-pi1*(g1+A11*pi1+A12*(1-pi1)));
% initial conditions
ics = [];
% known initial conditions
known_ics = [];
% 1 output
h = pi1;
save('cLV_2pi_0u_hpi1','x','p','h','f','u','ics','known_ics');

%%%%%%%%%%%%%%%%%%%%% Two species external disturbance %%%%%%%%%%%%%%%%%%%%
clear;
syms g1...
A11 A12...
B11...
u1...
pi1...
% 1 state (pi2 = 1 - pi1)
x = pi1;
% 1 input
u = u1;
% parameters
p =[g1; A11; A12; B11];
% dynamic equations
f = pi1*(g1+A11*pi1+A12*(1-pi1)+B11*u1-pi1*(g1+A11*pi1+A12*(1-pi1)+B11*u1));
% initial conditions
ics = [];
% known initial conditions
known_ics = [];
% 1 output
h = pi1;
save('cLV_2pi_1u_hpi1','x','p','h','f','u','ics','known_ics');

%%%%%%%%%%%%%%%%%% Three species no external disturbance %%%%%%%%%%%%%%%%%%
clear;
syms g1 g2...
A11 A12 A13 A21 A22 A23...
pi1 pi2...
% 2 states (pi3 = 1 - pi1 - pi2)
x = [pi1;pi2];
% inputs
u = [];
% parameters
p =[g1; g2; A11; A12; A13; A21; A22; A23];
% dynamic equations
f = [pi1*(g1+A11*pi1+A12*pi2+A13*(1-pi1-pi2)-pi1*(g1+A11*pi1+A12*pi2+A13*(1-pi1-pi2))-pi2*(g2+A21*pi1+A22*pi2+A23*(1-pi1-pi2)));
pi2*(g2+A21*pi1+A22*pi2+A23*(1-pi1-pi2)-pi1*(g1+A11*pi1+A12*pi2+A13*(1-pi1-pi2))-pi2*(g2+A21*pi1+A22*pi2+A23*(1-pi1-pi2)))];
% initial conditions
ics = [];
% known initial conditions
known_ics = [];
% outputs
h = [pi1;pi2];
save('cLV_3pi_0u_hpi1pi2','x','p','h','f','u','ics','known_ics');
h = pi1;
save('cLV_3pi_0u_hpi1','x','p','h','f','u','ics','known_ics');
h = pi2;
save('cLV_3pi_0u_hpi2','x','p','h','f','u','ics','known_ics');

%%%%%%%%%%%%%%%%%%%% Three species external disturbance %%%%%%%%%%%%%%%%%%%
clear;
syms g1 g2...
A11 A12 A13 A21 A22 A23...
B11 B21...
u1...
pi1 pi2...
% 2 states (pi3 = 1 - pi1 - pi2)
x = [pi1;pi2];
% inputs
u = [u1];
% parameters
p =[g1; g2; A11; A12; A13; A21; A22; A23; B11; B21];
% dynamic equations
f = [pi1*(g1+A11*pi1+A12*pi2+A13*(1-pi1-pi2)+B11*u1-pi1*(g1+A11*pi1+A12*pi2+A13*(1-pi1-pi2)+B11*u1)-pi2*(g2+A21*pi1+A22*pi2+A23*(1-pi1-pi2)+B21*u1));
pi2*(g2+A21*pi1+A22*pi2+A23*(1-pi1-pi2)+B21*u1-pi1*(g1+A11*pi1+A12*pi2+A13*(1-pi1-pi2)+B11*u1)-pi2*(g2+A21*pi1+A22*pi2+A23*(1-pi1-pi2)+B21*u1))];
% initial conditions
ics = [];
% known initial conditions
known_ics = [];
% outputs
h = [pi1;pi2];
save('cLV_3pi_1u_hpi1pi2','x','p','h','f','u','ics','known_ics');
h = pi1;
save('cLV_3pi_1u_hpi1','x','p','h','f','u','ics','known_ics');
h = pi2;
save('cLV_3pi_1u_hpi2','x','p','h','f','u','ics','known_ics');
