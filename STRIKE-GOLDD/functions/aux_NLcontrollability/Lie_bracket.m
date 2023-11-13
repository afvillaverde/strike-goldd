% Auxiliary function that computes the Lie bracket of (v1,v2) w.r.t. x.
% If v1 and/or v2 contain several vectors, it computes several brackets:
function brackets = Lie_bracket(v1,v2,x,n,jac_v1) 
    num_brcs = size(v1,2)*size(v2,2);
    brackets = zeros(n,num_brcs,'sym'); 
    ibracket = 1;
    for ivec2=1:size(v2,2)
        for ivec1=1:size(v1,2)
            brackets(:,ibracket) = jacobian(v2(:,ivec2),x)*v1(:,ivec1)...
            -jac_v1(:,(ivec1-1)*n+1:ivec1*n)*v2(:,ivec2); 
            ibracket = ibracket +1;
        end
    end
end