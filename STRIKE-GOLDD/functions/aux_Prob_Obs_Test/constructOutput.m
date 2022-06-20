%--------------------------------------------------------------------------
% construction of the output
%--------------------------------------------------------------------------

function Output=constructOutput(h,x,p,m,n,q,VSMOut)

aux=jacobian(h,[x;p]);
dSysOndV=aux(1:m,1:n);
if ~isempty(p)
    dSysOndP=aux(1:m,n+1:n+q);
    aux=dSysOndV*VSMOut(1:n,1:q)+dSysOndP;
    Output=[aux,dSysOndV*VSMOut(1:n,q+1:n+q)];
else
    Output=dSysOndV*VSMOut(1:n,1:n+q);
end