% Markov models with 3 states, 8, 12 or 10 parameters, 1 output, 1 input.
% The last one ('Markov_3cycle_db_reduced') is analysed in:
% "Input-Dependent Structural Identifiability of Nonlinear Systems"
% AF Villaverde, ND Evans, MJ Chappell, JR Banga
% IEEE Control Systems Letters 3 (2), 272-277

clear

syms x1 x2 x3 a12 a21 b12 b21 a23 a32 b23 b32 a13 a31 b13 b31 u1
x    = [x1 x2 x3].';   
h    = x1;               
u    = u1;
ics  = [];
known_ics = [0 0 0];

p   = [a12 a21 b12 b21 a23 a32 b23 b32].'; 
k12 = exp(a12 + b12*u);
k21 = exp(a21 + b21*u);
k23 = exp(a23 + b23*u);
k32 = exp(a32 + b32*u);

M = [0 k12 0; k21 0 k23; 0 k32 0];
for i = 1:length(M)
   M(i,i) = -sum(M(i,:));
end
f = M.'*x;
save('Markov_cco','x','h','p','f','u','ics','known_ics');

h = x2;
save('Markov_coc','x','h','p','f','u','ics','known_ics');

% Make the 3-cycle, i.e. include link from state 1 to state 3:
h   = x1; % puts it back to x1, not important
p   = [a12 a21 b12 b21 a23 a32 b23 b32 a13 a31 b13 b31].';
k13 = exp(a13 + b13*u);
k31 = exp(a31 + b31*u);
M = [0 k12 k13; k21 0 k23; k31 k32 0];
for i = 1:length(M)
   M(i,i) = -sum(M(i,:));
end
f = M.'*x;
save('Markov_3cycle','x','h','p','f','u','ics','known_ics');

% use detailed balance to remove one transition, k31.
% should have k12*k23*k31 = k21*k13*k32
% so a12+a23+a31 = a21+a13+a32 and similar for the b's
a31 = a21 + a13 + a32 - a12 - a23;
b31 = b21 + b13 + b32 - b12 - b23;
k31 = exp(a31 + b31*u);
p    = [a12 a21 b12 b21 a23 a32 b23 b32 a13 b13].';
% Matrix as before:
M = [0 k12 k13; k21 0 k23; k31 k32 0];
for i = 1:length(M)
   M(i,i) = -sum(M(i,:));
end
f = M.'*x;
save('Markov_3cycle_db','x','h','p','f','u','ics','known_ics');

%remove x3 from the dynamical system:
x3 = 1 - x1 - x2;
x  = [x1 x2 x3].';
f = M.'*x;
x = [x1 x2].';
f = f(1:2);
save('Markov_3cycle_db_reduced','x','h','p','f','u','ics','known_ics');
