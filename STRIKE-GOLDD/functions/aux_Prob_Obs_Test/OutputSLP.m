%--------------------------------------------------------------------------
% evaluation on the equations with the specific values to obtain and 
% evaluate the SLP
%--------------------------------------------------------------------------

function SLPOut=OutputSLP(Output,p,p_esp,u,u_esp,k,k_esp,x,Sol, ...
    VSMOut,SM,t,Order,Myprime,logic)

SLPOut=subs(Output,p,p_esp); 
SLPOut=subs(SLPOut,u,u_esp); 
SLPOut=subs(SLPOut,k,k_esp);
SLPOut=subs(SLPOut,x,Sol);
SLPOut=subs(SLPOut,VSMOut,SM);

kout=conj(symvar(SLPOut))';
[ch, original_index] = ismember(t,kout); 
if ch~=0
    kout(original_index)=[]; 
end
if ~isempty(kout)
    kout_esp = randi([0,Myprime],length(kout),1);
    SLPOut=subs(SLPOut,kout,kout_esp);
end

SLPOut=taylor(SLPOut,t,'order',Order); 
SLPOut=polmod(SLPOut,Myprime,logic); 
