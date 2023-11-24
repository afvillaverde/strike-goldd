%--------------------------------------------------------------------------
% File that creates the model of the NFKB signalling pathway. 
% It stores it in a mat-file named NFKB.mat.
% The model is taken from:
%--------------------------------------------------------------------------
% Lipniacki T, Paszek P, Brasier AR, Luxon B, Kimmel M. (2004)
% Mathematical model of NF-B regulatory module. 
% J Theor Biol. 2004; 228(2):195{215
%--------------------------------------------------------------------------
clear all;

syms t1 t2 c3a c4a c5 k1 k2 k3 kprod kdeg i1 e2a i1a...
    a1 a2 a3 c1a c2a c5a c6a c1 c2 c3 c4 kv e1a c1c c2c c3c...
    y1 y2 y3 y4 y5 y6...
    x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15...
    u1 

% states:
x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15].';

% outputs:
y1 = x7;
y2 = x10+x13;
y3 = x9;
y4 = x1+x2+x3;
y5 = x2;
y6 = x12;
h = [y1 y2 y3 y4 y5 y6].';

% parameters:
p = [t1,t2,c3a,c4a,c5,k1,k2,k3,kprod,kdeg,i1,e2a,i1a,...
    a1,a2,a3,c1a,c2a,c5a,c6a,c1,c2,c3,c4,kv,e1a,c1c,c2c,c3c].';

% input:
u = u1;

% dynamic equations:
f = [kprod-kdeg*x1-k1*x1*u1; 
    -k3*x2-kdeg*x2-a2*x2*x10+t1*x4-a3*x2*x13+t2*x5+(k1*x1-k2*x2*x8)*u1;
    k3*x2-kdeg*x3+k2*x2*x8*u1;
    a2*x2*x10-t1*x4;
    a3*x2*x13-t2*x5;
    c6a*x13-a1*x6*x10+t2*x5-i1*x6;
    i1*kv*x6-a1*x11*x7;
    c4*x9-c5*x8;
    c2+c1*x7-c3*x9;
    -a2*x2*x10-a1*x10*x6+c4a*x12-c5a*x10-i1a*x10+e1a*x11;
    -a1*x11*x7+i1a*kv*x10-e1a*kv*x11;
    c2a+c1a*x7-c3a*x12;
    a1*x10*x6-c6a*x13-a3*x2*x13+e2a*x14;
    a1*x11*x7-e2a*kv*x14;
    c2c+c1c*x7-c3c*x15];

% initial conditions:
ics  = [];  

% which initial conditions are known:
known_ics = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];   

save('NFKB','x','p','h','f','u','ics','known_ics');