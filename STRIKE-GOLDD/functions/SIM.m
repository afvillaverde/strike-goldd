%--------------------------------------------------------------------------
% Function that analyses the identifiability and observability of
% nonlinear models with an implementation of the method presented in: 
% Castro, M., and De Boer, R.J. "Testing structural identifiability by a 
% simple scaling method." PLoS Comput Biol 16(11): e1008248..
%--------------------------------------------------------------------------

function SIM(modelname,opts,prev_ident_pars,nmf)

%==========================================================================
% Load model:

% Convert model to multi-experiment form (Optional):
if opts.multiexp == 1
    ME_analysis(modelname,opts);
    modelname=strcat(modelname,'_',num2str(opts.multiexp_numexp),'Exp');
end

load(modelname)  %#ok<LOAD> 
fprintf('\n Analyzing the %s model with the SIM algorithm...\n',modelname);

timeSIM=tic;
%==========================================================================
% Initialize variables:
identifiables  = [];       % identifiable parameters.
nonidentif     = [];       % unidentifiable parameters.
obs_states     = [];       % observable states.
unobs_states   = [];       % unobservable states.
obs_inputs     = [];       % observable inputs.
unobs_inputs   = [];       % unobservable inputs.
% lastrank       = NaN;
% unidflag       = 0;
% isFISPO        = 0;
    
% ==========================================================================
% Remove parameters that have already been classified as identifiable:
if exist('prev_ident_pars','var') == 1
    for np=1:numel(prev_ident_pars)
        [~, original_index] = ismember(prev_ident_pars(np),p); 
        p(original_index)=[];                
    end
end
    
%==========================================================================
% Dimensions of the problem: 
m    = numel(h); % number of outputs 
n    = numel(x); % number of states
q    = numel(p); % number of unknown parameters

if size(h,2)>size(h,1),h=h.';end
if size(x,2)>size(x,1),x=x.';end
if size(p,2)>size(p,1),p=p.';end
if exist('w','var')
    nw = numel(w); % number of unknown inputs
else 
    nw = 0;
    w = [];
end
if exist('u','var')
    nu = numel(u); % number of known inputs
else 
    nu = 0;
    u = [];
end

fprintf('\n >>> The model contains:\n %d states:\n %s',n,char(x));
fprintf('\n %d outputs:\n %s',m,char(h));
fprintf('\n %d known inputs:\n %s',nu,char(u));
fprintf('\n %d unknown inputs:\n %s',nw,char(w));
fprintf('\n %d parameters:\n %s',q,char(p));


%==========================================================================
% Check which states are directly measured, if any: 
saidas     = cell(m,1);
ismeasured = zeros(1,n); 
for i=1:m
    saidas{i} = char(h(i));
end
for i=1:n
    ismeasured(i) = sum(strcmp(char(x(i)),saidas));
end

meas_x_indices = find(ismeasured);       % indices of the measured states
unmeas_x_indices = find((1-ismeasured)); % indices of the unmeasured states
meas_x = x(meas_x_indices);              % names of the measured states

%==================================================================
% Scaling factors
scaling_params = arrayfun(@(i) sym(['u_p' num2str(i)]), 1:q, 'UniformOutput', false);
scaling_states = arrayfun(@(i) sym(['u_x' num2str(i)]), unmeas_x_indices, 'UniformOutput', false);

% Parameter scaling: u_pi * parameter
p_scaled = cell(q, 1);
for i = 1:q
    p_scaled{i} = scaling_params{i} * p(i);
end

% Scaling of unmeasured states: u_xi * state
x_scaled = x;
for i = 1:length(unmeas_x_indices)
    idx = unmeas_x_indices(i);
    x_scaled = subs(x_scaled, x(idx), scaling_states{i} * x(idx));
end

%==========================================================================
% Expand system equations and separate into summands
fexpand = expand(f);           % Expand products and remove parentheses
f_terms = cell(size(fexpand)); % Cell array to store individual terms

% Separate the summands in each component as functionally independent terms
for i = 1:length(fexpand)
    f_terms{i} = split_terms_abs(fexpand(i));
end

f_terms2 = f_terms; % Duplicate array for substitutions of parameters and unmeasured states
                    % and division by scaling factors in unmeasured states equations
   
for i=1:length(f_terms)
    for j=1:length(f_terms{i})
        f_terms{i}{j}=str2sym(f_terms{i}{j});
        idx_in_unmeas = find(unmeas_x_indices == i);
        if ~isempty(idx_in_unmeas)
            f_terms2{i}{j} = (f_terms{i}{j}) / scaling_states(idx_in_unmeas);
        end
    end  
end

% Substitute parameters and unmeasured states with their scaled versions
p_scaled_vec = vertcat(p_scaled{:});  % Convert cell array into a symbolic column vector
for i = 1:length(f_terms2)
    for j=1:length(f_terms2{i})
        if ischar(f_terms2{i}{j})
            f_terms2{i}{j} = str2sym(f_terms2{i}{j});
        end
        term = sym(f_terms2{i}{j});
        vars_in_term = symvar(term);
        % Substitute parameters 
        if any(ismember(vars_in_term, p))
            term = subs((term), p, p_scaled_vec);
        end
        % Substitute unmeasured states
        if any(ismember(vars_in_term, x(unmeas_x_indices)))
            term = subs((term), x(unmeas_x_indices), x_scaled(unmeas_x_indices));
        end
        f_terms2{i}{j} = term; 
    end
end

%==========================================================================
% Generate scaling equations from original and scaled terms
n=0;
scaling_eqs_storage = {};  % Cell to store all generated equations
for i = 1:numel(f_terms)
    for j = 1:numel(f_terms2{i})
        n=n+1;
        orig_term = f_terms{i}{j};
        scaled_term = f_terms2{i}{j};
        scaling_eqs = (orig_term) == scaled_term;
        scaling_eqs_storage{n} = scaling_eqs;
    end
end

%==========================================================================
% Simplify system equations by removing common factors
scaling_simplified_storage={};
for i=1:numel(scaling_eqs_storage)
    eq=scaling_eqs_storage{i};
    eq_str = string(eq);   
    parts = strsplit(eq_str, '==');
    lhs = str2sym(parts(1)); 
    rhs = str2sym(parts(2)); 
    if isequal(lhs,rhs)
        eq_simpl=(lhs/rhs == rhs/lhs);
    else
        [num_lhs, den_lhs] = numden(lhs);
        [num_rhs, den_rhs] = numden(rhs);
        new_lhs = num_lhs * den_rhs;
        new_rhs = num_rhs * den_lhs;
        eq2= (new_lhs == new_rhs);
        lhs_fact = factor(new_lhs);
        rhs_fact = factor(new_rhs);
        while true         
            % Find common factors
            factores_comunes = intersect(lhs_fact, rhs_fact);
            if isempty(factores_comunes)
                break;
            end
            divisor = prod(factores_comunes);
            % Divide both sides of the equation by the common factor
            new_lhs= new_lhs/divisor;
            new_rhs= new_rhs/divisor;
            eq_simpl = (new_lhs == new_rhs);
            lhs_fact = factor(new_lhs);
            rhs_fact = factor(new_rhs);
        end
    end 
    scaling_simplified_storage{i} = eq_simpl;
end

%==========================================================================
% Substitute measured states with 1 to avoid redundancy
vars = [scaling_params, scaling_states];
vars_sym = [vars{:}];
ecuaciones_importantes = subs([scaling_simplified_storage{:}], x(meas_x_indices), ones(numel(meas_x_indices), 1));

% Separate equations into two groups: 
% ec_impor   → equations of standard (algebraic) type
% ec_imporc  → equations involving more complex terms (exponentials, logs, trigonometric functions, etc.)
ec_impor = {};
ec_imporc = {};

for i=1:numel(ecuaciones_importantes)
    eq=ecuaciones_importantes(i);
    eq_str = string(eq);   
    parts = strsplit(eq_str, '==');
    lhs = str2sym(parts(1)); 
    rhs = str2sym(parts(2)); 
    es= simplify(lhs-rhs);
   
    if ~isequal(char(es), '0')
        % Equations with exponential, logarithmic, or trigonometric functions
         if has(rewrite(es,'exp'), 'exp') ||... 
             has(rewrite(es,'exp'), 'log') ||... 
             has(rewrite(es,'exp'), 'sin') ||... 
             has(rewrite(es,'exp'), 'cos') ||...
             has(rewrite(es,'exp'), 'tan') ||...
             has(rewrite(es,'exp'), 'asin') ||...
             has(rewrite(es,'exp'), 'acos') ||...
             has(rewrite(es,'exp'), 'atan')
            ec_imporc{end+1} = ecuaciones_importantes(i);
         else
            ec_impor{end+1} = ecuaciones_importantes(i);
         end 
    end
end

% Solve complex equations separately
if ~isempty(ec_imporc)
    sol2 = solve([ec_imporc{:}], vars_sym, 'ReturnConditions', true); 
else
    sol2=0;
end 

% Solve standard equations
sol = solve([ec_impor{:}], vars_sym, 'ReturnConditions', true);  

%==========================================================================
% Classify parameters based on scaling solutions
% Create symbolic variables and associate them with original parameters
scaling_params2 = cell(1,length(scaling_params));
for i = 1:length(scaling_params)
    scaling_params2{i}.symbol = sym(['u_p' num2str(i)]); % Scaling factor symbol
    scaling_params2{i}.original = p(i);                  % Original parameter
    var = scaling_params2{i}.symbol;
    var_name = char(var);  
    ovar = scaling_params2{i}.original;
    ovar_name = char(ovar);
    % If the scaling solution equals 1 → parameter is structurally identifiable
    if isfield(sol, var_name) && (all(sol.(var_name) == sym(1))) ||... 
        isfield(sol2, var_name) && (all(sol2.(var_name) == sym(1))) ;
        identifiables{end+1} = ovar; 
    end 
    % If the scaling solution differs from 1 → parameter is unidentifiable
    if isfield(sol, var_name) && (all(sol.(var_name) ~= sym(1))) ||...
        isfield(sol2, var_name) && (all(sol.(var_name) ~= sym(1)));
        if any(has(f,ovar_name) )
            nonidentif{end+1}=ovar; end
    end
end

%==========================================================================
% Classify states based on scaling solutions
% Create symbolic variables and associate them with unmeasured states
n=1;
scaling_states2 = cell(1,length(scaling_states));
for i = unmeas_x_indices
    scaling_states2{n}.symbol = sym(scaling_states(n)); % Scaling factor symbol
    scaling_states2{n}.original_name = x(i);            % Original state name
    var = scaling_states{n};
    var_name = char(var);
    ovar = scaling_states2{n}.original_name;
    n=n+1;
    if isfield(sol, var_name) && (all(sol.(var_name) == sym(1))) || isfield(sol2, var_name) && (all(sol.(var_name) == sym(1))); %isequal((sol.(var_name)), sym(1))
        obs_states{end+1} = ovar; 
    end
    if isfield(sol, var_name) && (all(sol.(var_name) ~= sym(1))) || isfield(sol2, var_name) && (all(sol.(var_name) ~= sym(1)));    
        unobs_states{end+1} = ovar;
    end
end

%==========================================================================
% ------------------------ 
%        RESULTS
% ------------------------
fprintf('\n\n ------------------------ \n');
fprintf(' >>> RESULTS SUMMARY:\n');
fprintf(' ------------------------ \n');

fprintf('\n >>> NOTE THAT THE SIM ALGORITHM PROVIDES A NECESSARY BUT NOT SUFFICIENT CONDITION FOR IDENTIFIABILITY AND OBSERVABILITY.');
fprintf('\n     THEREFORE, CLAIMS ABOUT IDENTIFIABILITY OR OBSERVABILITY SHOULD BE TAKEN ONLY AS POSSIBILITIES, NOT GUARANTEES.\n');

	if numel(identifiables) == numel(p)
		fprintf('\n >>> The model is structurally identifiable:');
		if numel(identifiables)>0, fprintf('\n     All its parameters are structurally identifiable:\n      %s \n',strjoin(cellfun(@char, identifiables, 'UniformOutput', false), ', ')); end
        if numel(obs_states) == numel(x) && numel(obs_states)>0,fprintf('\n     All its unknown inputs are observable.'); end
	else    
			fprintf('\n >>> The model is structurally unidentifiable.');
            if numel(identifiables)>0, fprintf('\n >>> These parameters are identifiable:\n      %s \n',strjoin(cellfun(@char, identifiables, 'UniformOutput', false), ', ')); end
            if numel(nonidentif)>0,fprintf('\n >>> These parameters are unidentifiable:\n      %s \n',strjoin(cellfun(@char, nonidentif, 'UniformOutput', false), ', ')); end 
    end
		
	if numel(obs_states)>0 && numel(obs_states) ~= numel(x), fprintf('\n >>> These states are observable (and their initial conditions, if considered unknown, are identifiable):\n      %s ',strjoin(cellfun(@char, obs_states, 'UniformOutput', false), ', ')); end 
	if numel(unobs_states)>0, fprintf('\n >>> These states are unobservable (and their initial conditions, if considered unknown, are unidentifiable):\n      %s ',strjoin(cellfun(@char, unobs_states, 'UniformOutput', false), ', ')); end
  
    if numel(meas_x)>0,        fprintf('\n >>> These states are directly measured:\n      %s ',char(meas_x)); end
	if numel(obs_inputs)>0,    fprintf('\n >>> These unmeasured inputs are observable:\n      %s ',char(obs_inputs)); end
	if numel(unobs_inputs)>0,  fprintf('\n >>> These unmeasured inputs are unobservable:\n      %s ',char(unobs_inputs)); end
	if numel(u)>0,             fprintf('\n >>> These inputs are known:\n      %s ',char(u)); end
    
%==========================================================================
% Execution time
totaltime = toc(timeSIM);
fprintf('\n Total execution time: %d \n\n',totaltime);

%==========================================================================
% Save results:
resultsname = sprintf('id_results_%s_%s',modelname,date);   
fullresultsname = strcat(nmf,filesep,'results',filesep,resultsname);
save(fullresultsname);


end

