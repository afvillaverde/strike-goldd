% FitzHughNagumo model, as analysed in the high-gain observer paper.

clear;

% 2 states (the 3rd 'state' is the unknown input)
syms x1 x2 
x = [x1; x2];

% 1 output
h = x1;

% 1 known input
syms u

% no unknown input
w = [];

% unknown parameters
syms alpha
p = alpha;

% 2 known constants
syms epsilon eta

% auxiliary functions
syms f1 f2 g1 g2
f1 = -x2;
g1 = x1-x1^3/3;
f2 = epsilon*alpha*x1;
g2 = -epsilon*(x2+eta);

% dynamic equations
f = [f1+g1;
     f2+g2+epsilon*alpha*x1];

% initial conditions
ics  = []; 
known_ics = [0,0];

save('FHN_observer','x','p','h','u','w','f','ics','known_ics')
