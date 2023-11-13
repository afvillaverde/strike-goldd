%=========================================================================%
%= Function that checks the Accesibility Rank Condition ===============%
function ctrl_ARC(modelname,opts,n,nu,f1,f2,jac_f1,jac_f2,x,x0)

tARC = tic;
fprintf('\n Computing the ACCESIBILITY RANK CONDITION (ARC)...\n Computing Lie bracket number');
last_distr=sym('distr',[n,sum(1:nu)]);
last_distr(:,1:nu)=Lie_bracket(f1,f2,x,n,jac_f1);
last_index=nu;
for j=1:nu-1
    new_index=last_index+nu-j;
    last_distr(:,last_index+1:new_index)=Lie_bracket(f2(:,j),f2(:,j+1:end),x,n,jac_f2(:,(j-1)*n+1:j*n));
end
last_distr  =  sum(last_distr,2) ;
Cx(:,1:nu+2)  = [f1,f2,last_distr];
%****** Prov. test: evaluate at x = 0 *************************************
increaseLie = 1;
iLie        = 1;
lasttime = 0;
while increaseLie == 1 && lasttime < opts.maxtime
    fprintf(' %d... ',iLie)
    new_brackets = [Lie_bracket(f1,last_distr,x,n,jac_f1), ... % Construction 
                    Lie_bracket(f2,last_distr,x,n,jac_f2)];    % of brackets
    last_distr = sum(new_brackets,2);
    Cx= [Cx,last_distr]; %#ok<AGROW> 
    if last_distr == zeros(n,1)
        fprintf('   => The ARC is not met');    
        fprintf(' => it does not inform about the accessibility of the model %s.\n', modelname);
        return
    end
    if iLie>=ceil((n-1)/nu) % If nu>1, it may not be necessary to calculate n-1 Lie brackets. The minimum we need to calculate is (n-1)/nu
        ctrl_rank   = rank(subs(Cx,x,x0));
        if ctrl_rank == n 
            increaseLie = 0;
        end
    end
    iLie = iLie+1    ;   
    lasttime = toc(tARC);
end      

if iLie-1>=ceil((n-1)/nu)
    if lasttime >= opts.maxtime
        fprintf('\n    => More Lie brackets would be needed to see if the model is accessible.');
        fprintf('\n    However, the maximum computation time allowed is reached.');
        fprintf('\n    You can increase it by changing <<opts.maxtime>> (currently opts.maxtime = %d)\n',opts.maxtime);
        return
    end
end
fprintf('\n rank(C(x)) = %d',ctrl_rank);
if ctrl_rank == n
    fprintf('   => Model %s is accessible.\n', modelname);    
else
    if lasttime >= opts.maxtime
        fprintf('\n    => More Lie brackets would be needed to see if the model is accessible.');
        fprintf('\n    However, the maximum computation time allowed is reached.');
        fprintf('\n    You can increase it by changing <<opts.maxtime>> (currently opts.maxtime = %d)\n',opts.maxtime);
    else
        fprintf('   => The ARC is not met');    
        fprintf(' => it does not inform about the accessibility of the model %s.\n', modelname);   
    end
    fprintf(' (ARC computation completed in %d seconds.)\n ',toc(tARC));
end
end