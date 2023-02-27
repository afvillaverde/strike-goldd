% A phage cocktail model analized in:
% "Structural identifiability and observability of microbial community"
% S Díaz-Seoane, E Sellán, AF Villaverdeclear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% General version %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
syms r1 r2 m Kc E KD beta alpha fi v KI KN...
S R I PS PR...
Pc roSt roRt...
% 5 states
x = [S;R;PS;I;PR];
% inputs
u = [roSt,roRt];
% parameters
p =[r1; r2; m; Kc; E; KD; beta; alpha; fi; v; KI; KN; Pc];
% dynamic equations
f = [r1*S*(1-((S+R)/Kc))*(1-m)-S*fi*PS/(1+PS/Pc)-((E*I*S)/(1+(S+R)/KD));
r2*R*(1-((S+R)/Kc))+m*r1*S*(1-((S+R)/Kc))-R*fi*PR/(1+PR/Pc)-((E*I*R)/(1+(S+R)/KD));
beta*S*fi*PS/(1+PS/Pc)-fi*S*PS-v*PS+roSt;
alpha*I*(1-(I/KI))*((S+R)/(S+R+KN));
beta*R*fi*PR/(1+PR/Pc)-fi*R*PR-v*PR+roRt];
% initial conditions
ics = [];
% known initial conditions
known_ics = [0,0];
% outputs
h = [S;R;PS;I;PR];
save('PC_hSRPsIPr','x','p','h','f','u','ics','known_ics');
h = [S;R];
save('PC_hSR','x','p','h','f','u','ics','known_ics');
h = [S;R;I];
save('PC_hSRI','x','p','h','f','u','ics','known_ics');
h = [PS;PR];
save('PC_hPsPr','x','p','h','f','u','ics','known_ics');
h = [PS;I;PR];
save('PC_hPsIPr','x','p','h','f','u','ics','known_ics');
h = [S;PS;I];
save('PC_hSPsI','x','p','h','f','u','ics','known_ics');
h = [R;I;PR];
save('PC_hRIPr','x','p','h','f','u','ics','known_ics');
h = [S;R;PS;PR];
save('PC_hSRPsPr','x','p','h','f','u','ics','known_ics');
h = [S;I];
save('PC_hSI','x','p','h','f','u','ics','known_ics');
h = [R;I];
save('PC_hRI','x','p','h','f','u','ics','known_ics');
h = I;
save('PC_hI','x','p','h','f','u','ics','known_ics');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scaled version %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
syms r1 r2 m kCD kPD E1 beta alpha Y v q kND...
x1 x2 x3 x4 x5...
u1t u2t...
% 5 states
x = [x1;x2;x3;x4;x5];
% inputs
u = [u1t,u2t];
% parameters
p =[r1; r2; m; kCD; kPD; kND; E1; beta; alpha; Y; v; q];
% dynamic equations
f = [r1*x1*(1-((x1+x2)/kCD))*(1-m)-kPD*x1*Y*x3/(1+x3)-((E1*x4*x1)/(1+x1+x2));
r2*x2*(1-((x1+x2)/kCD))+r1*x1*(1-((x1+x2)/kCD))*m-kPD*x2*Y*x5/(1+x5)-((E1*x4*x2)/(1+x1+x2));
beta*x1*Y*x3/(1+x3)-Y*x1*x3-v*x3+q*u1t;
alpha*x4*(1-x4)*((x1+x2)/(x1+x2+kND));
beta*x2*Y*x5/(1+x5)-Y*x2*x5-v*x5+q*u2t];
% initial conditions
ics = [];
% known initial conditions
known_ics = [0,0];
% outputs
h = [x1;x2;x3;x4;x5];
save('PC_rel_hx1x2x3x4x5','x','p','h','f','u','ics','known_ics');
h = [x1;x2];
save('PC_rel_hx1x2','x','p','h','f','u','ics','known_ics');
h = [x1;x2;x4];
save('PC_rel_hx1x2x4','x','p','h','f','u','ics','known_ics');
h = [x3;x5];
save('PC_rel_hx3x5','x','p','h','f','u','ics','known_ics');
h = [x3;x4;x5];
save('PC_rel_hx3x4x5','x','p','h','f','u','ics','known_ics');
h = [x1;x3;x4];
save('PC_rel_hx1x3x4','x','p','h','f','u','ics','known_ics');
h = [x2;x4;x5];
save('PC_rel_hx2x4x5','x','p','h','f','u','ics','known_ics');
h = [x1;x2;x3;x5];
save('PC_rel_hx1x2x3x5','x','p','h','f','u','ics','known_ics');
h = [x1;x4];
save('PC_rel_hx1x4','x','p','h','f','u','ics','known_ics');
h = [x2;x4];
save('PC_rel_hx2x4','x','p','h','f','u','ics','known_ics');
h = [x4];
save('PC_rel_hx4','x','p','h','f','u','ics','known_ics');