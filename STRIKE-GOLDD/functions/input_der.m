        
function [wlvector,wlvector_dot,ulvector,ulvector_dot] = input_der(u,w,ind,opts)

uw = [u;w];  
nu = numel(u);
nw = numel(w);

% Create array of derivatives of measured & unmeasured inputs:
for num_in=1:numel(uw) 
    uwlarray(num_in,:) = [uw(num_in);[sym(strcat(char(uw(num_in)),sprintf('_d')),[1;ind+1])].'];
end
% Create vector of derivatives of measured inputs:
if nu>0
    ulvector = reshape(uwlarray(1:nu,1:end-1),[],1);            
end
% Create vector of derivatives of unmeasured inputs:
if nw>0
    wlvector = reshape(uwlarray(nu+1:end,1:end-1),[],1);          
end
% Set to zero the input derivatives of order > maximum allowed:
input_der_prov = [uwlarray,zeros(numel(uw),1)];
for ind_u=1:nu
    input_der_prov(ind_u,(opts.nnzDerU(ind_u)+2):end)=0;
end
for ind_w=1:nw
    input_der_prov(nu+ind_w,(opts.nnzDerW(ind_u)+2):end)=0;
end
uwlarray = input_der_prov(:,1:end-1);
% Replace the corresponding terms in the derivative vectors:
if nu>0
    ulvector_dot = reshape(uwlarray(1:nu,2:end),[],1);
end
if nw>0
    wlvector_dot = reshape(uwlarray(nu+1:end,2:end),[],1);
end