clear all;
syms x1 x2 x3 x4 x5 Gin D kover kdeg yh ya yg kg Kg ka Ka kpts kacs ....
     Oa Og n m l rgupP raupP raoverP rgupC raupC raoverC
    
% states:
x    = [x1; x2; x3; x4; x5]; %G=x1 A=x2 BP=x3 BC=x4 H=x5

% outputs:
h    =  [x1; x2; x3+x5];  

% inputs:
u    = [Gin];

% parameters:
%p = [yg; ya; kover; kg; Kg; ka; Ka; n; m; Oa; Og; l; kpts; kacs];
p = [yg; ya; kover; kg; ka; n; m; Oa; Og];

% productor
rgupP = kg*(x1/(x1+Kg))*((Oa^n)/((x2^n)+(Oa^n)));
raupP = ka*(x2/(x2+Ka))*((Og^m)/((rgupP^m)+(Og^m)));
raoverP = kover*(rgupP-l);

%cleaner
rgupC = kpts*(x1/(x1+Kg))*((Oa^n)/((x2^n)+(Oa^n)));
raoverC = kover*(rgupC-l);
raupC = (ka*(x2/(x2+Ka))*((Og^m)/((rgupC^m)+(Og^m))))+(kacs*(x2/(x2+kacs)));


% dynamic equations:    
f    = [-(rgupP)*x3-rgupC*x4+D*(Gin-x1);        
       (raoverP-raupP)*x3+(raoverC-raupC)*x4-D*x2;    
       (1-yh)*(yg*rgupP+ya*(raupP-raoverP))*x3-kdeg*x3-D*x3; 
       (yg*rgupC+ya*(raupC-raoverC))*x4-kdeg*x4-D*x4;
       yh*(yg*rgupP+ya*(raupP-raoverP))*x3-kdeg*x5-D*x5];

    
% initial conditions:    
ics  = [19.5;0;0.1;0.1;0];

% which initial conditions are known:
known_ics = [1;1;1;1;1];

save('microbial_community','x','h','u','p','f','ics','known_ics');
