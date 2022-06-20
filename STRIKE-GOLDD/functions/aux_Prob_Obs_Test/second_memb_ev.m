%--------------------------------------------------------------------------
% decomposition in submatrix of the SLP
%--------------------------------------------------------------------------

function [invsep_esp,logder_esp,sndmem_esp,system_esp]= ...
    second_memb_ev(SLP_ev,t,Myprime,n,q,Order)

SLP_ev=taylor(SLP_ev,t,'order',Order);  

%allocate memory
invsep_esp=zeros(Order,n*n);
logder_esp=zeros(Order,n*n);
sndmem_esp=zeros(Order,n);
system_esp=zeros(Order,n*(n+q));

% inverse of separent
ind=1;
for i=1:n
    for j=2:n+1
        if SLP_ev(i,j)==0
            aux=0;
        else
            aux=coeffs(SLP_ev(i,j),'All');
            aux=ratmod(aux,Myprime);
        end
        invsep_esp(Order-length(aux)+1:1:Order,ind)=aux;
        ind=ind+1;
    end
end

% logarithmic derivative
ind=1;
for i=1:n
    for j=n+2:2*n+1
        if SLP_ev(i,j)==0
            aux=0;
        else
            aux=coeffs(SLP_ev(i,j),'All');
        end
        aux=ratmod(aux,Myprime);
        logder_esp(Order-length(aux)+1:1:Order,ind)=aux;
        ind=ind+1;
    end
end

% second member
ind=1;
for i=1:n
    if SLP_ev(i,1)==0
        aux=0;
    else
        aux=coeffs(SLP_ev(i,1),'All');
    end
    aux=ratmod(aux,Myprime);
    sndmem_esp(Order-length(aux)+1:1:Order,ind)=aux;
    ind=ind+1;
end

% system
ind=1;
for i=1:n
    for j=2+2*n:1+2*n+n+q
        if SLP_ev(i,j)==0
            aux=0;
        else
            aux=coeffs(SLP_ev(i,j),'All');
        end
        aux=ratmod(aux,Myprime);
        system_esp(Order-length(aux)+1:1:Order,ind)=aux;
        ind=ind+1;
    end
end

