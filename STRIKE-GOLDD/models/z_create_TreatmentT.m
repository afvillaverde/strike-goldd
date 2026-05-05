clear;

% 4 states
syms S In Tr
x = [S; In; Tr];

% no inputs
u = [];

% 5 unknown parameters
syms b d a g nu N
p = [b; d; a; g; nu];

% 2 outputs
h = [
    Tr
];

% dynamic equations
f = [
    -b*S*In/N - d*b*S*Tr/N;

     b*S*In/N + d*b*S*Tr/N - (a + g)*In;

     g*In - nu*Tr;

];

% initial conditions
ics = [];
known_ics = [];

save('TreatmentT','x','p','h','f','u','ics','known_ics');