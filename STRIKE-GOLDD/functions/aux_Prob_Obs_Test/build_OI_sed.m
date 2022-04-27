%--------------------------------------------------------------------------
% Function that builds the generalized observability-identifiability
% matrix of the model specified in the 'options' script. the matrix is
% built using Sedoglavic's probabilistic approach.
%--------------------------------------------------------------------------

function [onx]=build_OI_sed(Myprime,x,p,u,h,f,n,q,nu,m)

% variable t
t=sym('t');

% random values for parameters, varibles and known imputs
Sol = randi([0,Myprime],n,1); 
p_esp = randi([0,Myprime],q,1);

if nu~=0
    u_esp_coef = randi([0,Myprime],nu,2+n+q);
    u_esp=sum(t.^[1+n+q:-1:0].*u_esp_coef,2); %#ok<NBRAK> 
else
    u_esp=[];
end

% vector with derivates of the states and states
xd=sym('xd',[n 1]);
for i=1:n
    xd(i) = sym([char(x(i)) 'd']);
end

% variational sensitibity matrix and derivative
VSMd = sym('VSMd_%d_%d',[n,n+q]); 
VSM = sym('VSM_%d_%d',[n,n+q]); 


% computation of linear variational system
LVS=calcLVS(x,xd,p,f,VSM,VSMd,n,q);

%computation of SLP
SLP=subs(LVS,p,p_esp);
SLP=subs(SLP,u,u_esp);

% sensitivity matrix
SM=[zeros(n,q),eye(n)];
Order=2;

%evaluation of the power series Sol(...) on SLP
[invsep_esp,logder_esp,sndmem_esp]=second_memb_ev(Sol,SM,xd,x,VSM,VSMd,SLP,t,Myprime,n,q,Order);

%==========================================================================
% Loop

% initialization
OneMoreLoop=1;
Bound=n+q;

ind=0;

fprintf('\n >>> Current order of calculations: \n');

while OneMoreLoop==1

    ind=ind+1;

    if Order==Bound
        OneMoreLoop=0;
    end

    HomSol=calcHomSol(logder_esp,Order,n,Myprime);

    InvHomSolInvA=calcInvHomSolInvA(HomSol,invsep_esp,Order,n,Myprime,t);

    VarOfCte=calcVarOfCte(HomSol,InvHomSolInvA,sndmem_esp,n,n,n,1,Order,Myprime);

    Sol=Sol+VarOfCte;

    OldOrder=Order;
    if 2*Order<Bound
        Order=2*Order;
    else
        Order=Bound;
    end
    NewOrder=Order;
    
    [invsep_esp,logder_esp,sndmem_esp,system_esp]=second_memb_ev(Sol,SM,xd,x,VSM,VSMd,SLP,t,Myprime,n,q,Order);

    Order=OldOrder;

    VarOfCte=calcVarOfCte(HomSol,InvHomSolInvA,system_esp,n,n,n,n+q,Order,Myprime);

    SM=SM+VarOfCte;

    Order=NewOrder;

    disp(Order)

end

%==========================================================================
% Construct output

VSMOut=sym('VSMOut',[n,n+q]);
Output=contructOutput(h,x,p,m,n,q,VSMOut);

SLPOut=OutputSLP(Output,p,p_esp,u,u_esp,x,Sol,VSMOut,SM,t,Order,Myprime);


%==========================================================================
% Jacobian Matrix
onx=constructonx(SLPOut,OldOrder);



