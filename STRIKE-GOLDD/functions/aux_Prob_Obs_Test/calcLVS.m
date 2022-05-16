%==========================================================================
% Construction of linear variational sistem for SLP:
% Gives a matrix with the system, the inverse of A, the inverse of A&*B
% (logarithmic derivative) and partial derivatives of System w.r.t. initial
% conditions end parameters.
function LVS=calcLVS(x,xd,p,f,VSM,VSMd,n,q)

tic
%--------------------------------------------------------------------------
% Initial caculations

% system into vector field
P = xd - f;

% System into polynomial formulation:
% Gives a polinomial formulation of the function P

for i=1:length(P)
    P(i)=rational2polinomial(P(i));
end

%--------------------------------------------------------------------------
% Computation of the partial derivative w.r.t. the states, the state
% derivatives end the parameters

Mpd=jacobian(P,[xd;x;p]);

%--------------------------------------------------------------------------
% Computation of the linear variational sistem

Mpdxd=Mpd(1:n,1:n);
Mpdx=Mpd(1:n,n+1:2*n);

invMpdxd=Mpdxd;
invMpdxd=invMpdxd-diag(diag(invMpdxd))+diag(1./diag(Mpdxd));

logder=invMpdxd*Mpdx;

LVS=[P,invMpdxd,logder];
if ~isempty(p)
    Mpdp=Mpd(1:n,2*n+1:2*n+q);
    LVS=[LVS,Mpdxd*VSMd(:,1:q)+Mpdx*VSM(:,1:q)+Mpdp];
    LVS=[LVS,Mpdxd*VSMd(:,q+1:q+n)+Mpdx*VSM(:,q+1:q+n)];
else
    LVS=[LVS,Mpdxd*VSMd+Mpdx*VSM];
end