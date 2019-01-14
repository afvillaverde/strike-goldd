        
function [xaug,faug] = add_input_der(xaug,faug,u,w,ind,opts)

for num_in=1:numel(w) 
    xaug = [xaug;[sym(strcat(char(w(num_in)),sprintf('_d%d',ind-1)))] ];
    if opts.nnzDerW(num_in)>=ind    
        faug = [faug;[sym(strcat(char(w(num_in)),sprintf('_d%d',ind)))] ];
    else
        faug = [faug;0];
    end 
end

for num_in=1:numel(u) 
    xaug = [xaug;[sym(strcat(char(u(num_in)),sprintf('_d%d',ind-1)))] ];
    if opts.nnzDerU(num_in)>=ind    
        faug = [faug;[sym(strcat(char(u(num_in)),sprintf('_d%d',ind)))] ];
    else
        faug = [faug;0];
    end 
end
