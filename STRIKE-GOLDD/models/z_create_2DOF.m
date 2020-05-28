% Model from Maes, K., Chatzis, M., Lombaert, G., 2019. 
% "Observability of nonlinear systems with unmeasured inputs."
% Mech. Syst. Signal Process. 130, 378–394..

clear;

% 4 states
syms x1 x2 dx1 dx2
x = [x1;x2;dx1;dx2];

% 1 known input and 1 unknown input
syms F1 F2;
u = F1;
w = F2;

% 3 unknown parameters 
syms k1 k2 dk1 m1 m2 c1 c2
p =[k1; dk1; m2];

% 2 outputs
h = [x1; 1/m2*( k2*(x1-x2)+c2*(dx1-dx2)+F2 )];

% initial conditions
ics  = []; 
known_ics = [0,0,0,0];

% dynamic equations
f = [dx1;
     dx2;
     1/m1*( -(k1+dk1*x1)*x1+k2*(x2-x1)-c1*dx1+c2*(dx2-dx1)+F1 );
     1/m2*( k2*(x1-x2)+c2*(dx1-dx2)+F2 )];

save('2DOF','x','p','u','w','h','f','ics','known_ics');
