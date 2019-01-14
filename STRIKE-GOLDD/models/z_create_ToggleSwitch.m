%--------------------------------------------------------------------------
% Genetic toggle switch by
% Lugagne JB, Carrillo SS, Kirch M, Köhler A, Batt G, Hersen P. 
% "Balancing a genetic toggle switch by real-time feedback control and periodic forcing." 
% Nature Communications. 2017 Nov 17;8(1):1671.
%--------------------------------------------------------------------------
clear all;

% 2 states:
syms x1 x2
x = [x1 x2].';

% 10 parameters:
syms k01 k1 tatc natc ntetr ...
     k02 k2 tiptg niptg nlaci
p = [k01 k1 tatc natc ntetr ...
     k02 k2 tiptg niptg nlaci].';
% Nominal values: p = [ 0.032*0.9726/(0.1386*31.94*0.0165) ...
%                       8.30*0.9726/(0.1386*31.94*0.0165) ...
%                       11.65 ...
%                       2 ...
%                       2 ...
%                       0.119*1.17/(0.1386*30*0.0165) ...
%                       2.06*1.17/(0.1386*30*0.0165) ...
%                       0.0906 ...
%                       2 ...
%                       2 ];

% 2 outputs:
h = x;

% 2 inputs:
syms u1 u2;
u = [u1 u2];

% Unknown inputs:
w = [];

% dynamic equations:
f = [ 
	k01 + k1/(1+(x2/(1+(u1/tatc)^natc))^ntetr) - x1;
    k02 + k2/(1+(x1/(1+(u2/tiptg)^niptg))^nlaci) - x2;
];

% initial conditions:
ics  = [];   

% which initial conditions are known:
known_ics = [0,0]; 

% Known input case (with parameterization as defined above):
save('ToggleSwitch_known_inputs','x','p','u','w','h','f','ics','known_ics');

% Unknown input case (reparameterization):
u = [];
syms w1 w2
w = [w1 w2];
p = [k01 k1 ntetr ...
     k02 k2 nlaci].';
f = [ 
	k01 + k1/(1+(x2/(1+w1))^ntetr) - x1;
    k02 + k2/(1+(x1/(1+w2))^nlaci) - x2;
];
save('ToggleSwitch_unknown_inputs_reparameterized','x','p','u','w','h','f','ics','known_ics');