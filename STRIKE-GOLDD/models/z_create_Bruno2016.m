% File that creates the model from Bruno et al, J Exp Bot 2016

clear all;

%--------------------------------------------------------------------------
% Model features common to all experiments:
%--------------------------------------------------------------------------

% 7 states
syms beta cry zea beta10 OHbeta10 betaio OHbetaio 
x = [beta cry zea beta10 OHbeta10 betaio OHbetaio];

% 7 unknown parameters (6 kinetic constants + 1 scaling) (+ init conditions):
syms kbeta kcryOH kcrybeta kzea kbeta10 kOHbeta10 szea
p = [kbeta kcryOH kcrybeta kzea kbeta10 kOHbeta10];

% no inputs:
u = [];
w = [];

% dynamic equations: 
f = [-kbeta*beta;
-kcryOH*cry-kcrybeta*cry;
-kzea*zea;
kbeta*beta+kcryOH*cry-kbeta10*beta10;
kcrybeta*cry+kzea*zea-kOHbeta10*OHbeta10;
kbeta*beta+kcrybeta*cry+kbeta10*beta10;
kcryOH*cry+kzea*zea+kOHbeta10*OHbeta10];

%--------------------------------------------------------------------------
% Model features varying across experiments:
%--------------------------------------------------------------------------

%--Data set 1 ('beta1'):---------------------------------------------------

% 2 outputs:
h = [beta; beta10];

% initial conditions:
ics       = [1,0,0,0,0,0,0];    
known_ics = [0,1,1,1,1,1,1];

save('bruno_1_beta1','x','h','p','f','u','w','ics','known_ics'); 


%--Data set 2 ('beta2'):---------------------------------------------------

% 2 outputs:
h = [beta; beta10];

% initial conditions:
ics       = [1,0,0,0,0,0,0];    
known_ics = [0,1,1,1,1,1,1];

% parameters: all rate parameters multiplied by szea
p = [kbeta kcryOH kcrybeta kzea kbeta10 kOHbeta10 szea];

% dynamic equations: all rate parameters multiplied by szea
f = [-szea*(kbeta*beta);
-szea*(kcryOH*cry-kcrybeta*cry);
-szea*kzea*zea;
szea*(kbeta*beta+kcryOH*cry-kbeta10*beta10);
szea*(kcrybeta*cry+kzea*zea-kOHbeta10*OHbeta10);
szea*(kbeta*beta+kcrybeta*cry+kbeta10*beta10);
szea*(kcryOH*cry+kzea*zea+kOHbeta10*OHbeta10)];

save('bruno_2_beta2','x','h','p','f','u','w','ics','known_ics'); 


%--Data set 3 ('cry'):-----------------------------------------------------

% 3 outputs:
h = [cry; beta10; OHbeta10];

p = [kbeta kcryOH kcrybeta kzea kbeta10 kOHbeta10];

% dynamic equations: 
f = [-kbeta*beta;
-kcryOH*cry-kcrybeta*cry;
-kzea*zea;
kbeta*beta+kcryOH*cry-kbeta10*beta10;
kcrybeta*cry+kzea*zea-kOHbeta10*OHbeta10;
kbeta*beta+kcrybeta*cry+kbeta10*beta10;
kcryOH*cry+kzea*zea+kOHbeta10*OHbeta10];

% initial conditions:
ics       = [0,1,0,0,0,0,0];    
known_ics = [1,0,1,1,1,1,1];

save('bruno_3_cry','x','h','p','f','u','w','ics','known_ics'); 


%--Data set 4 ('zea'):-----------------------------------------------------

% 2 outputs:
h = [zea; OHbeta10];

% initial conditions:
ics       = [0,0,1,0,0,0,0];    
known_ics = [1,1,0,1,1,1,1];

% parameters: all rate parameters multiplied by szea
p = [kbeta kcryOH kcrybeta kzea kbeta10 kOHbeta10 szea];

% dynamic equations: all rate parameters multiplied by szea
f = [-szea*(kbeta*beta);
-szea*(kcryOH*cry-kcrybeta*cry);
-szea*kzea*zea;
szea*(kbeta*beta+kcryOH*cry-kbeta10*beta10);
szea*(kcrybeta*cry+kzea*zea-kOHbeta10*OHbeta10);
szea*(kbeta*beta+kcrybeta*cry+kbeta10*beta10);
szea*(kcryOH*cry+kzea*zea+kOHbeta10*OHbeta10)];

save('bruno_4_zea','x','h','p','f','u','w','ics','known_ics'); 


%--Data set 5 ('beta10'):-----------------------------------------------------

% 1 output:
h = beta10;

% parameters:
p = [kbeta kcryOH kcrybeta kzea kbeta10 kOHbeta10];

% dynamic equations: 
f = [-kbeta*beta;
-kcryOH*cry-kcrybeta*cry;
-kzea*zea;
kbeta*beta+kcryOH*cry-kbeta10*beta10;
kcrybeta*cry+kzea*zea-kOHbeta10*OHbeta10;
kbeta*beta+kcrybeta*cry+kbeta10*beta10;
kcryOH*cry+kzea*zea+kOHbeta10*OHbeta10];

% initial conditions:
ics       = [0,0,0,1,0,0,0];    
known_ics = [1,1,1,0,1,1,1];

save('bruno_5_beta10','x','h','p','f','u','w','ics','known_ics'); 


%--Data set 6 ('OHbeta10'):-----------------------------------------------------

% 1 output:
h = OHbeta10;

% parameters:
p = [kbeta kcryOH kcrybeta kzea kbeta10 kOHbeta10];

% dynamic equations: 
f = [-kbeta*beta;
-kcryOH*cry-kcrybeta*cry;
-kzea*zea;
kbeta*beta+kcryOH*cry-kbeta10*beta10;
kcrybeta*cry+kzea*zea-kOHbeta10*OHbeta10;
kbeta*beta+kcrybeta*cry+kbeta10*beta10;
kcryOH*cry+kzea*zea+kOHbeta10*OHbeta10];

% initial conditions:
ics       = [0,0,0,0,1,0,0];    
known_ics = [1,1,1,1,0,1,1];

save('bruno_6_OHbeta10','x','h','p','f','u','w','ics','known_ics')


%-- 5 outputs, all parameters known ('five_out'):---------------------------------------------------

% 5 outputs:
h = [beta; cry; zea; beta10; OHbeta10];

% initial conditions:
ics       = [1,1,1,1,1,0,0];    
known_ics = [0,0,0,0,0,1,1];

% parameters: all known
p = [];

% dynamic equations: all rate parameters multiplied by szea
f = [-szea*(kbeta*beta);
-szea*(kcryOH*cry-kcrybeta*cry);
-szea*kzea*zea;
szea*(kbeta*beta+kcryOH*cry-kbeta10*beta10);
szea*(kcrybeta*cry+kzea*zea-kOHbeta10*OHbeta10);
szea*(kbeta*beta+kcrybeta*cry+kbeta10*beta10);
szea*(kcryOH*cry+kzea*zea+kOHbeta10*OHbeta10)];

save('bruno_five_outputs','x','h','p','f','u','w','ics','known_ics'); 
