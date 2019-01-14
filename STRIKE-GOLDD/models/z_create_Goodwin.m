%--------------------------------------------------------------------------
% File that creates the Goodwin oscillator model. 
% It stores it in a mat-file named Goodwin.mat.
% The model is taken from:
%--------------------------------------------------------------------------
% Goodwin BC (1965) Oscillatory behavior in enzymatic control processes. 
% Advances in enzyme regulation, 3, 425-437.
%--------------------------------------------------------------------------
clear all;

syms x1 x2 x3 ...
     a b AA sigma alpha beta gamma delta
 
% states: 
x    = [x1; x2; x3]; 

% outputs:
h    = [x1];    

% no input:
u    = [];

% parameters:
p    = [a; b; AA; sigma; alpha; beta; gamma; delta]; 

% dynamic equations:
f    = [-b*x1 + a/(AA+x3^sigma);        
        alpha*x1-beta*x2;    
        gamma*x2-delta*x3]; 
    
% initial conditions:    
ics  = [0.3617, 0.9137, 1.3934];  

% which initial conditions are known:
known_ics = [1,1,1];

save('Goodwin','x','h','u','p','f','ics','known_ics');