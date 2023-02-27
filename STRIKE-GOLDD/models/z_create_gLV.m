% Generalized Lotka-Volterra models analized in:
% "Structural identifiability and observability of microbial community"
% S Díaz-Seoane, E Sellán, AF Villaverde

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Two species %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
syms r1 r2 ...
x1 x2 ...
beta11 beta12 beta21 beta22 ...
% 2 states
x = [x1; x2];
% inputs
u = [];
% parameters
p =[r1; r2; beta11; beta12; beta21; beta22];
% dynamic equations
f = [r1*x1+beta11*x1^2+beta12*x1*x2;
r2*x2+beta21*x1*x2+beta22*x2^2];
% initial conditions
ics = [];
% known initial conditions
known_ics = [0,0];
% outputs
h = [x1;x2];
save('gLV_2x_hx1x2','x','p','h','f','u','ics','known_ics');
h = x1;
save('gLV_2x_hx1','x','p','h','f','u','ics','known_ics');
h = x2;
save('gLV_2x_hx2','x','p','h','f','u','ics','known_ics');
% relative
h = [x1/(x1+x2);x2/(x1+x2)];
save('gLV_rel_2x_hx1x2','x','p','h','f','u','ics','known_ics');
h = x1/(x1+x2);
save('gLV_rel_2x_hx1','x','p','h','f','u','ics','known_ics');
h = x2/(x1+x2);
save('gLV_rel_2x_hx2','x','p','h','f','u','ics','known_ics');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Three species %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
syms r1 r2 r3...
x1 x2 x3...
beta11 beta12 beta13 beta21 beta22 beta23 beta31 beta32 beta33...
% 2 states
x = [x1; x2; x3];
% inputs
u = [];
% parameters
p =[r1; r2; r3; beta11; beta12; beta13; beta21; beta22; beta23; beta31; beta32; beta33];
% dynamic equations
f = [r1*x1+beta11*x1^2+beta12*x1*x2+beta13*x1*x3;
r2*x2+beta21*x1*x2+beta22*x2^2+beta23*x2*x3;
r3*x3+beta31*x1*x3+beta32*x2*x3+beta33*x3^2];
% initial conditions
ics = [];
% known initial conditions
known_ics = [0,0];
% outputs
h = [x1;x2;x3];
save('gLV_3x_hx1x2x3','x','p','h','f','u','ics','known_ics');
h = [x1;x2];
save('gLV_3x_hx1x2','x','p','h','f','u','ics','known_ics');
h = [x1;x3];
save('gLV_3x_hx1x3','x','p','h','f','u','ics','known_ics');
h = [x2;x3];
save('gLV_3x_hx2x3','x','p','h','f','u','ics','known_ics');
h = x1;
save('gLV_3x_hx1','x','p','h','f','u','ics','known_ics');
h = x2;
save('gLV_3x_hx2','x','p','h','f','u','ics','known_ics');
h = x3;
save('gLV_3x_hx3','x','p','h','f','u','ics','known_ics');
% relative
h = [x1/(x1+x2+x3);x2/(x1+x2+x3);x3/(x1+x2+x3)];
save('gLV_rel_3x_hx1x2x3','x','p','h','f','u','ics','known_ics');
h = [x1/(x1+x2+x3);x2/(x1+x2+x3)];
save('gLV_rel_3x_hx1x2','x','p','h','f','u','ics','known_ics');
h = [x1/(x1+x2+x3);x3/(x1+x2+x3)];
save('gLV_rel_3x_hx1x3','x','p','h','f','u','ics','known_ics');
h = [x2/(x1+x2+x3);x3/(x1+x2+x3)];
save('gLV_rel_3x_hx2x3','x','p','h','f','u','ics','known_ics');
h = x1/(x1+x2+x3);
save('gLV_rel_3x_hx1','x','p','h','f','u','ics','known_ics');
h = x2/(x1+x2+x3);
save('gLV_rel_3x_hx2','x','p','h','f','u','ics','known_ics');
h = x3/(x1+x2+x3);
save('gLV_rel_3x_hx3','x','p','h','f','u','ics','known_ics');
