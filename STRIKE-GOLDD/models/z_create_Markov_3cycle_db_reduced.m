% Markov model analysed in:
% "Input-Dependent Structural Identifiability of Nonlinear Systems"
% AF Villaverde, ND Evans, MJ Chappell, JR Banga
% IEEE Control Systems Letters 3 (2), 272-277
%--------------------------------------------------------------------------
clear
syms x1 x2 x3 a12 a21 b12 b21 a23 a32 b23 b32 a13 a31 b13 b31 u1 

% 1 output:
h = x1; 

% 1 input:
u = u1;

% 10 unknown parameters:
p = [a12 a21 b12 b21 a23 a32 b23 b32 a13 b13].';

% aux expressions:
k12 = exp(a12 + b12*u);
k21 = exp(a21 + b21*u);
k23 = exp(a23 + b23*u);
k32 = exp(a32 + b32*u);
k13 = exp(a13 + b13*u);

% Use detailed balance to remove one transition, k31.
% It should have k12*k23*k31 = k21*k13*k32,
% so a12+a23+a31 = a21+a13+a32 and similar for the b's
a31 = a21 + a13 + a32 - a12 - a23;
b31 = b21 + b13 + b32 - b12 - b23;
k31 = exp(a31 + b31*u);

% Matrix:
M = [0 k12 k13; k21 0 k23; k31 k32 0];
for i = 1:length(M)
   M(i,i) = -sum(M(i,:));
end

% Reduce model: remove x3:
x3 = 1 - x1 - x2;
x  = [x1 x2 x3].';
f  = M.'*x;

% 2 states:
x  = [x1 x2].';

% dynamic equations:
f  = f(1:2);

% initial conditions:
ics       = [];
known_ics = [0 0];

save('Markov_3cycle_db_reduced','x','h','p','f','u','ics','known_ics');