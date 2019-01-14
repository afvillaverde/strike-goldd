%--------------------------------------------------------------------------
% File that creates the Pitavastatin model with pseudo-steady state
% assumption, and stores it in a mat-file named Pitavastatin_ss.mat.
% The model is described in eqs. (109)-(110) in the paper:
%--------------------------------------------------------------------------
% Grandjean TR, Chappell MJ, Yates JW, Evans ND. (2014) Structural 
% identifiability analyses of candidate models for in vitro Pitavastatin
% hepatic uptake.
% Comput Methods Programs Biomed. 2014;114(3):e60-e69
%--------------------------------------------------------------------------
clear

syms x1 x3...                   % states
     k D1...                    % known constants
     r3 k3 k4 Vm Km T0          % unknown parameters
 
% states:
x    = [x1 x3].';              

% output:
h    = k*( T0*x1/(Km+x1) + x3); 

% no input:
u    = [];

% unknown parameters:
p    = [r3 k3 k4 Vm Km T0].';

% dynamic equations:
f    = [ k3*x3 - r3*x1 - Vm*x1/(Km+x1);        
         r3*x1 - (k3+k4)*x3 + Vm*x1/(Km+x1)];  

% initial conditions: 
ics = [D1, 0];

% which initial conditions are known:
known_ics = [1, 1]; 

save('Pitavastatin_ss','x','h','u','p','f','ics','known_ics');