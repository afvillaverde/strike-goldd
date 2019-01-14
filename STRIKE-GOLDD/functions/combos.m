%--------------------------------------------------------------------------
% Function that tries to find identifiable combinations of otherwise
% unidentifiable parameters
%--------------------------------------------------------------------------
function [parpde,stringpde] = combos(p_un,onx,n)

global p

[~, indices] = ismember(p_un, p);
indices = n + indices;
unid_onx = onx(:,indices);
parpde = null(unid_onx);
[numpars,numpde] = size(parpde);

for j=1:numpde
    stringpde{j} = '';
    firstflag = 1;
    for i=1:numpars
        if (parpde(i,j) ~= 0) && firstflag == 1
            stringpde{j} = strcat(stringpde{j},'(',char(parpde(i,j)),')*dF/d',char(p_un(i)));
            firstflag = 0;
        else if (parpde(i,j) ~= 0)
            stringpde{j} = strcat(stringpde{j},' + (',char(parpde(i,j)),')*dF/d',char(p_un(i)));
            end
        end
    end
    stringpde{j} = strcat(stringpde{j},' = 0');
end
