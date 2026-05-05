clear;

% 7 states
syms beta cry zea beta10 OHbeta10 betaio OHbetaio
x = [beta; cry; zea; beta10; OHbeta10; betaio; OHbetaio];

% 2 outputs
h = [
    beta;
    beta10
];

% no inputs
u = [];

% 7 unknown parameters
syms kbeta kcryOH kcrybeta kzea kbeta10 kOHbeta10
p = [kbeta; kcryOH; kcrybeta; kzea; kbeta10; kOHbeta10];

% dynamic equations
f = [
    -kbeta*beta;

    -kcryOH*cry - kcrybeta*cry;

    -kzea*zea;

     kbeta*beta + kcryOH*cry - kbeta10*beta10;

     kcrybeta*cry + kzea*zea - kOHbeta10*OHbeta10;

     kbeta*beta + kcrybeta*cry + kbeta10*beta10;

     kcryOH*cry + kzea*zea + kOHbeta10*OHbeta10
];

% initial conditions
ics = [];
known_ics = [];

save('Bruno2016','x','p','h','f','u','ics','known_ics');