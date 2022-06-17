%--------------------------------------------------------------------------
% Function that builds the generalized observability-identifiability
% matrix of the model specified in the 'options' script. the matrix is
% built using Sedoglavic's probabilistic approach.
%--------------------------------------------------------------------------

function [onx]=build_OI_sed(Myprime,x,p,u,h,f,n,q,nu,m,opts)

% variable t
t=sym('t');

% random values for parameters, state varibles and known imputs
Sol = randi([0,Myprime],n,1); 
p_esp = randi([0,Myprime],q,1);

if nu~=0
    if length(opts.nnzDerU) == 1
        if opts.nnzDerU == inf
            u_esp_coef = randi([0,Myprime],nu,1+n+q);
            u_esp=sum(t.^[n+q:-1:0].*u_esp_coef,2); %#ok<NBRAK>
        else
            u_esp_coef=randi([0,Myprime],nu,opts.nnzDerU+1);
            u_esp=sum(t.^[opts.nnzDerU:-1:0].*u_esp_coef,2); %#ok<NBRAK> 
        end
    else
        u_esp=sym('u_esp',[nu,1]);
        for ind_u=1:nu
            if opts.nnzDerU(ind_u) == inf
                u_esp_coef = randi([0,Myprime],nu,1+n+q);
                u_esp=sum(t.^[n+q:-1:0].*u_esp_coef,2); %#ok<NBRAK>
            else
                u_esp_coef=randi([0,Myprime],1,opts.nnzDerU(ind_u)+1);
                u_esp(ind_u)=sum(t.^[opts.nnzDerU(ind_u):-1:0].* ...
                    u_esp_coef,2); %#ok<NBRAK> 
            end
        end 
    end
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

% computation of SLP
SLP=subs(LVS,p,p_esp);
SLP=subs(SLP,u,u_esp);

% sensitivity matrix
SM=[zeros(n,q),eye(n)];
Order=2;

% derivation of variables and sensibility matrix
Sold=diff(Sol,t);
SMd=diff(SM,t);

% substitutions with initial conditions
SLP_ev=subs(SLP,[xd;x],[Sold;Sol]);
SLP_ev=subs(SLP_ev,VSM,SM);
SLP_ev=subs(SLP_ev,VSMd,SMd);

% check for mesured values
k=conj(symvar(SLP_ev))';
[ch, original_index] = ismember(t,k); 
if ch~=0
    k(original_index)=[]; 
end
k_esp=[];
if ~isempty(k)
    k_esp = randi([0,Myprime],length(k),1);
    SLP=subs(SLP,k,k_esp);
    SLP_ev=subs(SLP_ev,k,k_esp);
end

% Obtein power series Sol(...) on SLP with polynomials as vectors on
% columns of a matrix
[invsep_esp,logder_esp,sndmem_esp]= ...
    second_memb_ev(SLP_ev,t,Myprime,n,q,Order);

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
    
    %----------------------------------------------------------------------
    % Homogeneous resolution

    % integration of polynomials in a matrix
    HomSol=-intmatpoly(logder_esp,Order);

    % computation using a newton operator
    HomSol=ExpMatrixSeries(HomSol,Order,n,Myprime);

    %----------------------------------------------------------------------
    % Computation of InvHomSolInvA

    % computation using a newton operator
    InvHomSolInvA=InverseMatrixSeries(HomSol,Order,n,Myprime);
    
    % multiplication of InvHomSolInvA and inssep_esp
    InvHomSolInvA= ...
        mod(multmatpolytrun(InvHomSolInvA,invsep_esp,Order,n),Myprime);

    %----------------------------------------------------------------------
    % Computation of variation of constants for states
    VarOfCte= ...
        calcVarOfCte(HomSol,InvHomSolInvA,sndmem_esp,n,n,n,1,Order,Myprime);

    Sol=Sol+VarOfCte;

    %----------------------------------------------------------------------
    % New order
    OldOrder=Order;
    if 2*Order<Bound
        Order=2*Order;
    else
        Order=Bound;
    end
    NewOrder=Order;

    %----------------------------------------------------------------------
    % derivation of variables and sensibility matrix
    Sold=diff(Sol,t);
    SMd=diff(SM,t);

    %----------------------------------------------------------------------
    % substitutions with actualiced values
    SLP_ev=subs(SLP,[xd;x],[Sold;Sol]);
    SLP_ev=subs(SLP_ev,VSM,SM);
    SLP_ev=subs(SLP_ev,VSMd,SMd);
    
    %----------------------------------------------------------------------
    % Obtein power series Sol(...) on SLP with polynomials as vectors on
    % columns of a matrix
    [invsep_esp,logder_esp,sndmem_esp,system_esp]= ...
        second_memb_ev(SLP_ev,t,Myprime,n,q,Order);

    Order=OldOrder;

    %----------------------------------------------------------------------
    % Computation of variation of constants for variational systems
    VarOfCte=calcVarOfCte(HomSol,InvHomSolInvA,system_esp,n,n,n,n+q, ...
        Order,Myprime);

    SM=SM+VarOfCte;

    Order=NewOrder;

    disp(Order)

end

%==========================================================================
% Construct output
VSMOut=sym('VSMOut',[n,n+q]);
Output=constructOutput(h,x,p,m,n,q,VSMOut);

SLPOut=OutputSLP(Output,p,p_esp,u,u_esp,k,k_esp,x,Sol,VSMOut,SM,t, ...
    Order,Myprime);


%==========================================================================
% Jacobian Matrix
onx=constructonx(SLPOut,OldOrder);



