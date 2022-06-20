% -------------------------------------------------------------------------
% computation using a newton operator
%--------------------------------------------------------------------------

function [sol1,sol2]=ExpMatrixSeries(Mat,LocalOrder,Matsize,MyPrime)

if LocalOrder==1
    sol1=zeros(1,Matsize*Matsize);
    sol1((1:Matsize:Matsize*Matsize)+(0:Matsize-1))=1;
    sol2=zeros(1,Matsize*Matsize);
    sol2((1:Matsize:Matsize*Matsize)+(0:Matsize-1))=1;
else
    NewOrder=2^ceil(log2(LocalOrder)-1);
    aux=trunpoly(Mat,LocalOrder);
    [aux1,aux2]=ExpMatrixSeries(aux,NewOrder,Matsize,MyPrime);

    aux1=trunpoly(aux1,LocalOrder);
    aux2=trunpoly(aux2,LocalOrder);
    
    sol2=multmatpolytrun(-aux1,aux2,NewOrder,Matsize);
    sol2=mod(sol2,MyPrime);
    sol2(end,(1:Matsize:Matsize*Matsize)+(0:Matsize-1))= ...
        sol2(end,(1:Matsize:Matsize*Matsize)+(0:Matsize-1))+2;
    sol2=multmatpolytrun(aux2,sol2,NewOrder,Matsize);
    sol2=mod(sol2,MyPrime);

    
    Matder=dermatpoly(Mat,LocalOrder);
    aux1der=dermatpoly(aux1,LocalOrder);
    sol1=mod(trunpoly( ...
        multmatpolytrun(Matder,aux1,LocalOrder,Matsize)-aux1der ...
        ,LocalOrder),MyPrime);
    sol1=mod(multmatpolytrun(sol2,sol1,LocalOrder,Matsize),MyPrime);
    sol1=intmatpoly(sol1,LocalOrder);
    sol1=ratmod(sol1,MyPrime); 
    sol1(end,(1:Matsize:Matsize*Matsize)+(0:Matsize-1))= ...
        sol1(end,(1:Matsize:Matsize*Matsize)+(0:Matsize-1))+1;
    sol1=mod(multmatpolytrun(aux1,sol1,LocalOrder,Matsize),MyPrime);


end
