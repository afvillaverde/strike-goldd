% - "BIG" (betaIG) model by: Topp et al, J Theor Biol 2000 
% - Originally published in: Karin et al, Mol Syst Biol 2016
% - The two model versions defined here are analysed in: 
% Massonis, Banga, Villaverde. "Automatic reformulation method...", 2020.

clear;

% 3 states
syms G beta I
x = [G; beta; I];

% 1 output
h = G;

% 1 input
syms input;

% 5 unknown parameters 
syms p1 si gamma c alpha
p =[p1; si; gamma; c; alpha];

% known constants
muplus  = 0.021/(24*60); % turnover of functional mass
muminus = 0.025/(24*60);

% auxiliary functions
syms rhoG  lambdaplus lambdaminus 
rhoG        = G^2/(alpha^2+G^2); 
lambdaplus  = muplus/(1+(8.4/G)^1.7);   
lambdaminus = muminus/(1+(G/4.8)^8.5);

% dynamic equations
f = [input-(c+si*I)*G;
     beta*(lambdaplus-lambdaminus);     
     p1*beta*rhoG-gamma*I];

% initial conditions
ics  = []; 
known_ics = [0,0,0];

u = input;
w = [];
save('BIG_known_input','x','p','h','u','w','f','ics','known_ics');
u = [];
w = input;
save('BIG_unknown_input','x','p','h','u','w','f','ics','known_ics');
