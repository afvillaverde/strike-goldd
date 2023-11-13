%=========================================================================%
%= Function that checks the Linearization Condition for controllability ==%
function ctrl_LC(modelname,opts,n,f1,f2,x,x0)

tLC = tic;
fprintf('\n Computing the LINEARIZATION CONDITION (LC)... ');
A = jacobian(f1,x);
B = f2;
A = subs(A,x,x0);
B = subs(B,x,x0);
%**************************************************************************
C_array = B;
new_term = B;

if opts.numericLC == 1
%--- Replace symbolic arrays with random numbers:
    Myprime = nextprime(10^6);
    vars = symvar([A,B]);
    vars_num = randi([0,Myprime],size(vars));
    A_num = subs(A,vars,vars_num);
    new_term_num = subs(new_term,vars,vars_num);
    C_array_num = subs(C_array,vars,vars_num);
    %--- Calculate the Controllability array numerically:
    for i=1:n-1
%         fprintf('\n i = %d... ',i) % Uncomment for debugging
        new_term_num = A_num*new_term_num;
        C_array_num = [C_array_num,new_term_num]; %#ok<AGROW> 
    end
    %--- Calculate rank:    
    lin_rank = rank(C_array_num); 
else
    %--- Calculate the Controllability array symbolically:
    for i=1:n-1
        new_term = A*new_term;
        C_array = [C_array,new_term]; %#ok<AGROW> 
    end
    %--- Calculate rank:    
    lin_rank = rank(C_array);
end
%--- Print result:
fprintf('\n rank(C(x)) = %d',lin_rank);
if lin_rank == n
    fprintf('   => Model %s is short time locally controllable (STLC).\n', modelname);    
else
    fprintf('   => The LC is not met');    
    fprintf(' => it does not inform about the controllability of the model %s.\n', modelname);   
end
fprintf(' (LC computation completed in %d seconds.)\n ',toc(tLC));

end