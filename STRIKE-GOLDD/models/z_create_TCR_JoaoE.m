clear;

% 3 states
syms S T A
x = [S; T; A];

% no inputs
u = [];

% unknown parameters
syms lambda phi s k ki L
p = [lambda; phi; s; k; ki; L];

% 1 output
h = S/(lambda + 1) + ((T + A)*lambda/(lambda + 1));

% dynamic equations
f = [
    -lambda * phi * (S - T) + s*(1 - S);
     phi * (S - T) + s*(1 - T) - k*(T^5)*(L^5);
     k*(T^5)*(L^5) - ki*A
];

% initial conditions
ics = [];
known_ics = [];

save('TCR_JoaoE','x','p','h','f','u','ics','known_ics');