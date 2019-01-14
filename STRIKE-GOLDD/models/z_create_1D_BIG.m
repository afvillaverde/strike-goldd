% "BIG" (betaIG) model by Topp et al, J Theor Biol 2000 
% Originally published in: Karin et al, Mol Syst Biol 2016
% Corresponds to the model in Fig. 1D of: Villaverde & Banga, arXiv:1701.02562

clear;

% 3 states
syms x1 x2 x3
x = [x1; x2; x3];

% 1 output
h = x1;

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

save('1D_BIG','x','p','h','u','f','ics','known_ics');
% h = x3;
% save('1D_BIG_output_I','x','p','h','u','f','ics','known_ics');
% h = x2;
% save('1D_BIG_output_Beta','x','p','h','u','f','ics','known_ics');
% h = [x2;x3];
% save('1D_BIG_output_BetaI','x','p','h','u','f','ics','known_ics');
% h = [x1;x2];
% save('1D_BIG_output_BetaG','x','p','h','u','f','ics','known_ics');
% h = [x1;x3];
% save('1D_BIG_output_IG','x','p','h','u','f','ics','known_ics');

