%--------------------------------------------------------------------------
% evaluation on the equations with the specific values to obtain and 
% evaluate the SLP
%--------------------------------------------------------------------------

function SLPOut=OutputSLP(Output,p,p_esp,u,u_esp,k,k_esp,x,Sol, ...
    VSMOut,SM,t,Order,Myprime)

SLPOut=subs(Output,p,p_esp); 
SLPOut=subs(SLPOut,u,u_esp); 
SLPOut=subs(SLPOut,k,k_esp);
SLPOut=subs(SLPOut,x,Sol);
SLPOut=subs(SLPOut,VSMOut,SM);

SLPOut=taylor(SLPOut,t,'order',Order); 
SLPOut=polmod(SLPOut,Myprime); 