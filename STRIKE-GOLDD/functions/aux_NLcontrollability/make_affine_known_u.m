% - Convert model into affine in inputs form
% - Called by:
%   * ORC_DF
%   * ctrl-analysis
% - Outputs:
%   * is_affine   (= 1 if the system is affine, = 0 otherwise);
%   * varargout{1} = fxw;
%   * varargout{2} = fu;
%   * varargout{3} = hxw;
%   * varargout{4} = hu;
function [is_affine,varargout] = make_affine_known_u(modelname)

fprintf('\n >>> Converting model into affine-in-inputs form...');

%==========================================================================
% Initialize output variables in case to account for early termination:
is_affine    = 0;
varargout{1} = [];
varargout{2} = [];
if nargout >3
    varargout{3} = [];
    varargout{4} = [];
end

%==========================================================================
% Load model and read its dimensions:
load(modelname) %#ok<*LOAD>  %,'x','p','w','u','f','h'

% Number of states:
ns = numel(x);   
% Number of outputs:
if ~exist('h','var')
    h=x;
    warning('Output assumed to be x')
end
    m  = numel(h);
% Number of known inputs:
if exist('u','var') && numel(u)>0 
    nu = numel(u);
else
    u  = [];
    nu = 0;
end
% Number of unknown inputs:
if exist('w','var') && numel(w)>0 
    nw = numel(w);
else
    w  = [];
    nw = 0;
end
% Number of unknown parameters:
if exist('p','var') && numel(p)>0 
    np = numel(p);
else
    p  = [];
    np = 0;
end


%==========================================================================
% Check that the system is affine in known inputs and store coefficients:

% Initialize known inputs coefficients of system dynamics: 
fu = sym(zeros(ns,nu));
% Initialize known inputs coefficients of output dynamics:
hu = sym(zeros(m,nu));

syms coeff_u coefh_u
for j=1:nu
    % States:
    for i=1:ns
        try
            coeff_u = coeffs(f(i),u(j),'all');
            if numel(coeff_u)>2
                warning('Unable to convert system to affine in inputs form. System dynamics is not affine in known inputs.')
                return
            end
        catch
            warning('Unable to convert system to affine in inputs form. System dynamics is not affine in known inputs.')
            return
        end
        variables = symvar(coeff_u);
        nz_variables = simplify(variables-u(j)*ones(1,numel(variables)));
        if numel(find(nz_variables==0))~=0
            warning('Unable to convert system to affine in inputs form. System dynamics is not affine in known inputs.')
            return
        end
        % Store uk coefficients of system dynamics: 
        if numel(coeff_u)==2, fu(i,j)=coeff_u(1);end
    end
    % Outputs:
    for i=1:m
        try
            coefh_u = coeffs(h(i),u(j),'all');
            if numel(coefh_u)>2
                warning('Unable to convert system to affine in inputs form. Output dynamics is not affine in known inputs.')
                return
            end
        catch
            warning('Unable to convert system to affine in inputs form. Output dynamics is not affine in known inputs.')
            return
        end
        variables = symvar(coefh_u);
        nz_variables = simplify(variables-u(j)*ones(1,numel(variables)));
        if numel(find(nz_variables==0))~=0
            warning('Unable to convert system to affine in inputs form. Output dynamics is not affine in known inputs.')
            return
        end
        % Store uk coefficients of output dynamics:
        if numel(coefh_u)==2, hu(i,j)=coefh_u(1);end
    end
end

% State-unknown inputs 0-augmented system dynamics:
if size(u,2)>1, u = transpose(u); end
fxw = simplify(f-fu*u);
% State-unknown inputs output dynamics:
hxw = simplify(h-hu*u);

%==========================================================================
% Save affine model and return output(s):
affine_modelname = sprintf('affine_%s',modelname);
fullaffinename   = strcat(pwd,filesep,'models',filesep,affine_modelname);
if exist('ics','var')==0,ics=[];end
if exist('known_ics','var')==0,known_ics=[];end
save(fullaffinename,'x','u','w','p','f','h','fu','fxw','hu','hxw','ics','known_ics');

fprintf('\n >>> Conversion successful.\n');
is_affine    = 1;
varargout{1} = fxw;
varargout{2} = fu;
if nargout >3
    varargout{3} = hxw;
    varargout{4} = hu;
end

end