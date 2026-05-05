clear;

% 4 states
syms S E I R
x = [S; E; I; R];

% 6 unknown parameters
syms beta epsilon rho mu d
p = [beta; epsilon; rho; mu; d];

% no inputs
u = [];

% 1 output
h = mu*I;

% dynamic equations
f = [
    -beta*I*S;

     beta*I*S - epsilon*E;

     epsilon*E - (rho + mu)*I;

     rho*I - d*R
];

% initial conditions
ics = [];
known_ics = [];

save('SEIR16','x','p','h','f','u','ics','known_ics');