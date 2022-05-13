% construction of the output
function Output=constructOutput(h,x,p,w,m,n,q,nw,VSMOut)

Mhd=jacobian(h,[x;p;w]);
dSysOndV=Mhd(1:m,1:n);
aux=[];
if ~isempty(p)
    dSysOndP=Mhd(1:m,n+1:n+q);
    aux=dSysOndV*VSMOut(1:n,1:q)+dSysOndP;
end
if ~isempty(w)
    dSysOndW=Mhd(1:m,n+q+1:n+q+nw);
    aux=[aux,dSysOndV*VSMOut(1:n,q+1:q+nw)+dSysOndW];
end
Output=[aux,dSysOndV*VSMOut(1:n,q+nw+1:n+q+nw)];