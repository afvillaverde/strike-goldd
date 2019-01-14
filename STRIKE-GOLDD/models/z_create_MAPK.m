%--------------------------------------------------------------------------
% File that creates the model of the MAPK cascade with mixed feedback. 
% It stores it in a mat-file named MAPK.mat.
% The model is taken from:
%--------------------------------------------------------------------------
% Nguyen LK, Degasperi A, Cotter P, Kholodenko BN (2015) DYVIPAC: an 
% integrated analysis and visualisation framework to probe 
% multi-dimensional biological networks. 
% Scientific Reports 5.
%--------------------------------------------------------------------------
clear all;

syms k1 k2 k3 k4 k5 k6 ...
    ps1 ps2 ps3 ...
    s1t s2t s3t ...
    KK1 KK2 n1 n2 alpha  ...
    
% states:
x    = [ps1; ps2; ps3]; 

% outputs:
h    = x;  

% no input:
u    = [];

% parameters:
p    = [k1; k2; k3; k4; k5; k6;...
        s1t; s2t; s3t; KK1; KK2; n1; n2; alpha];  
    
% dynamic equations:    
f    = [k1*(s1t-ps1)*(KK1^n1)/(KK1^n1+ps3^n1)-k2*ps1;        
        k3*(s2t-ps2)*ps1*(1+(alpha*ps3^n2)/(KK2^n2+ps3^n2))-k4*ps2 ;    
        k5*(s3t-ps3)*ps2-k6*ps3 ]; 
    
% initial conditions:    
ics  = [];

% which initial conditions are known:
known_ics = [0,0,0];

save('MAPK','x','h','u','p','f','ics','known_ics');