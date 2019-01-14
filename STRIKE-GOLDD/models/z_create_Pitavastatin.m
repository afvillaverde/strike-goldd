%--------------------------------------------------------------------------
% File that creates the Pitavastatin model. 
% It stores it in a mat-file named Pitavastatin.mat.
% The model is described in eqs. (103)-(105) in the paper:
%--------------------------------------------------------------------------
% Grandjean TR, Chappell MJ, Yates JW, Evans ND. (2014) Structural 
% identifiability analyses of candidate models for in vitro Pitavastatin
% hepatic uptake.
% Comput Methods Programs Biomed. 2014;114(3):e60-e69
%--------------------------------------------------------------------------
clear

syms x1 x2 x3...                % states
     k D ...                    % known constants
     k1 k2 k3 k4 r1 r3 T0       % unknown parameters

% states:
x    = [x1 x2 x3].';   

% output:
h    = k*(x2 + x3);   

% no input:
u    = [];

% parameters:
p    = [k1 k2 k3 k4 r1 r3 T0].';

% dynamic equations:
f    = [ k3*x3 - r3*x1 - k1*x1*(T0-x2) + r1*x2;        
         k1*x1*(T0-x2) - (r1+k2)*x2;
         r3*x1 - (k3+k4)*x3 + k2*x2];  

% initial conditions: 
ics = [D, 0, 0];

% which initial conditions are known:
known_ics = [1, 1, 1]; 

save('Pitavastatin','x','h','u','p','f','ics','known_ics');