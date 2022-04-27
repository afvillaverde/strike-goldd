%computation of the variation of constants
function VarOfCte=calcVarOfCte(HomSol,InvHomSolInvA,sndmem_esp,nr1,nc1,nr2,nc2,Order,Myprime)

aux=mod(multmatpolytrun( HomSol, ...
        ratmod(-intmatpoly(mod( ...
        multmatpolytrun(InvHomSolInvA,sndmem_esp,Order,nr1,nc1,nr2,nc2) ...
        ,Myprime),Order),Myprime) ...
        ,Order,nr1,nc1,nr2,nc2),Myprime);
t=sym('t');

aux=t.^(Order-1:-1:0)*aux;

VarOfCte=sym('VarOfCte',[nr1,nc2]);
ind=1;
for i=1:nr1
    for j=1:nc2
        VarOfCte(i,j)=aux(ind);
        ind=ind+1;
    end
end