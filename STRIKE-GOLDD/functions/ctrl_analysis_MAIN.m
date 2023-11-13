% Controllability and Reachability analysis of Nonlinear Models
function ctrl_analysis_MAIN(varargin)
tStart = tic;

%==========================================================================
%=== Read options, add folders to path, and load model: ===================
switch nargin
    case 0
        [modelname,~,opts,~,~] = options;
    case 1
        copyfile(varargin{1},"current_options.m");
        [modelname,~,opts,~,~] = current_options;
        delete("current_options.m");
    case 4
        tests     = varargin{1};
        modelname = varargin{2};        
        opts.LC   = tests(1); 
        opts.ARC  = tests(2); 
        opts.LARC = tests(3); 
        opts.GSC  = tests(4);
        opts.numericLC = varargin{3};       % 0 or 1
        opts.maxtime = varargin{4};
    case 5
        tests     = varargin{1};
        modelname = varargin{2};        
        opts.LC   = tests(1); 
        opts.ARC  = tests(2); 
        opts.LARC = tests(3); 
        opts.GSC  = tests(4);
        opts.numericLC = varargin{3};       % 0 or 1
        opts.maxtime = varargin{4};
        x0 = varargin{5};
end
load(modelname) %#ok<*LOAD>
fprintf('\n Analysing model %s ', modelname);
fprintf('\n Preprocessing model for controllability/accessibility analysis...');

%==========================================================================
%=== Check if the model is affine in the control inputs: ==================
if exist('u','var') && numel(u)>0
    affine_modelname   = sprintf('affine_%s',modelname);
    affine_model_file  = strcat(pwd,filesep,'models',filesep,strcat(affine_modelname,'.mat'));
    exist_affine_model = exist(affine_model_file,'file');
    %--- If the model is not written using the affine notation, rewrite it:
    if exist_affine_model ~= 2
        [is_affine,fxw,fu] = make_affine_known_u(modelname);        
        if is_affine == 0
            fprintf('\n %s is not affine in the inputs. ', modelname);
            fprintf('Its controllability cannot be analysed. \n\n');
            return
        end            
    end 
    %--- Load the model written in the affine-in-control-inputs notation:
    modelname = affine_modelname;
    load(affine_modelname) %#ok<*LOAD>
    n  = numel(x);
    nu = numel(u);
    f1 = fxw; % = 'f(x)' in eq. (2.1) in the paper
    f2 = fu;  % = 'g(x)' in eq. (2.1) in the paper
    fprintf('\n Model %s has %d states and %d inputs.\n ', modelname,n,nu);
else
    fprintf('\n %s is an uncontrolled system. ', modelname);
    fprintf('Its controllability cannot be analysed. \n\n');
    return
end
%=========================================================================%
%=== Check for simbolic parameters and replace them ======================%
k=conj(symvar([f1,f2]))';
for i=1:length(x) % Check for parameters
    [ch, original_index] = ismember(x(i),k); 
    if ch~=0
        k(original_index)=[]; 
    end
end
if ~isempty(k) % Substitution of parameters
    k_esp = randi([1,11],length(k),1);
    f1=subs(f1,k,k_esp);
    f2=subs(f2,k,k_esp);
    fprintf('Parameters replaced with random values \n');
    for i=1:length(k)
        fprintf('%s = %d\n',k(i),k_esp(i));
    end
end
fprintf('(Preprocessing completed in %d seconds.)\n ',toc(tStart));
%=========================================================================%
%===================== Compute jacobians of f1 and f2 ====================%
t_jac=tic;
jac_f1 = jacobian(f1,x);
jac_f2=sym('d',[n,n*nu]);
for i=1:nu
    jac_f2(:,(i-1)*n+1:i*n)=jacobian(f2(:,i),x);
end
toc(t_jac)
%=========================================================================%
%=== Compute the equilibrium point if one is not given ===================%
logic_eq=0;
if nargin < 5
    x0=struct2cell(solve(f1==0,x,'Real',true));
    x0=subs(x,x,x0);
    logic_eq=1;
end
%=========================================================================%
%=== Check the conditions indicated in 'opts': ===========================%
n_x0=length(x0)/n;
for i=1:n_x0
    x0_i=x0(i:n_x0:(n-1)*n_x0+i);
    if logic_eq==1
        if n_x0>1
            fprintf('More than one equilibrium point');
        end
        fprintf('Condition check for the equilibrium point: \n');
        fprintf('    %s \n',x0_i);
    else
        fprintf('Condition check for the point: \n');
        fprintf('    %d \n',x0_i);
        fprintf('If the point is not an equilibrium point only the LARC results are valid. \n');
    end
    if opts.LC == 1
        ctrl_LC(modelname,opts,n,f1,f2,x,x0_i)
    end
    if opts.ARC == 1
        ctrl_ARC(modelname,opts,n,nu,f1,f2,jac_f1,jac_f2,x,x0_i)
    end
    if opts.LARC == 1 || opts.GSC == 1
        ctrl_LARC_GSC(modelname,opts,n,nu,f1,f2,jac_f1,jac_f2,x,x0_i)
    end
end

%==========================================================================
%=== Save results: ========================================================
resultsname = sprintf('ctrl_results_%s_%s',modelname,date);
fullresultsname = strcat(pwd,filesep,'results',filesep,resultsname);
save(fullresultsname);

totaltime = toc(tStart);
fprintf('\n Total execution time: %d ',totaltime);
fprintf('\n--------------------------------------------------------\n\n');

end % END function
