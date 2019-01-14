clear;

% 3 states
syms x1 x2 x3
x = [x1; x2; x3];

% 1 outputs
h = x3;

% 1 input
syms u u0;

% 5 unknown parameters 
syms p1 p2 p3 p4 p5
p =[p1; p2; p3; p4; p5];

% known constants
muplus  = 0.021/(24*60); % turnover of functional mass
muminus = 0.025/(24*60);

% auxiliary functions
syms rhoG  lambdaplus lambdaminus 
rhoG        = x1^2/(p5^2+x1^2); 
lambdaplus  = muplus/(1+(8.4/x1)^1.7);   
lambdaminus = muminus/(1+(x1/4.8)^8.5);

% dynamic equations
f = [u0+u-(p4+p2*x3)*x1;
     x2*(lambdaplus-lambdaminus);     
     p1*x2*rhoG-p3*x3];

% initial conditions
ics  = []; 
known_ics = [0,0,0];

save('1D_out_I_par_agcsp','x','p','h','f','ics','known_ics');