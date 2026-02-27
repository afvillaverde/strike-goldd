
function results = SIM_PDE(modelFile)
%% SIM method for analysing Identifiability of PDE models
timeSIM = tic;
Identifiable_parameters   = {};
unidentifiable_parameters = {};
observable_states         = {};
unobservable_states       = {};
Identifiable_combinations = {};
observed_state = {};

%% -------------------------------------------------
% Load required variables from .mat file
% (States, parameters, equations, observations, optional equations)
load(modelFile, 'x', 'p', 'f', ...
                'observed_vars', 'observation_eq', 'opt_eq');
states         = x;
params         = p;
eq             = f;

%% -------------------------------------------------
% Extract model name 
[~, modelName, ~] = fileparts(modelFile);

%% -------------------------------------------------
% Scaling states
stateNames = cellfun(@(v) regexprep(char(v), '\(.*\)', ''), ...
                     states, 'UniformOutput', false);

scaling_states = sym("u_" + string(stateNames));

% Scaling parameters
ParNames = cellfun(@(v) regexprep(char(v), '\(.*\)', ''), ...
                     params, 'UniformOutput', false);

param_scaling = sym("u_" + string(ParNames));

%% -------------------------------------------------
% Separate equation terms
for i = 1:length(eq)
    exp_eq = expand(eq{i});
    sep_eq = children(exp_eq);

    if length(sep_eq) == 2 && sep_eq{1}*sep_eq{2} == eq{i}
        separeted_eq{i} = {sep_eq{1}*sep_eq{2}};
    else
        separeted_eq{i} = sep_eq;
    end
end

%% -------------------------------------------------
% Apply scaling
for i = 1:length(eq)
    scaled_eq{i} = separeted_eq{i} / scaling_states(i);
end

scaled_params = params .* param_scaling;
scaled_states = states .* scaling_states;

for i = 1:length(scaled_eq)
    scaled_eq{i} = subs(scaled_eq{i}, params, scaled_params);
    scaled_eq{i} = subs(scaled_eq{i}, states, scaled_states);
end

%% -------------------------------------------------
% Construct invariance equations
invariance_eqs = sym([]);

for k = 1:length(separeted_eq)

    sep_k   = separeted_eq{k};
    scale_k = scaled_eq{k};

    for j = 1:length(sep_k)
        invariance_eqs(end+1,1) = sep_k{j} - scale_k(j) == 0;
    end

end

%% -------------------------------------------------
% Observed states scaling = 1
E = zeros(length(scaling_states),1);

for i = 1:length(states)
    for j = 1:length(observed_vars)
        if isequal(states{i}, observed_vars{j})
            invariance_eqs = [invariance_eqs; scaling_states(i) == 1];
            E(i) = 1;
        end
    end
end

%% -------------------------------------------------
% Observation equations
for i = 1:length(observation_eq)
    scaled_obs_eq = subs(observation_eq{i}, params, scaled_params);
    scaled_obs_eq = subs(scaled_obs_eq, states, scaled_states);

    invariance_eqs = [invariance_eqs;
                      observation_eq{i} - scaled_obs_eq == 0];
end

%% -------------------------------------------------
% Optional equations
for i = 1:length(opt_eq)
    invariance_eqs = [invariance_eqs; opt_eq{i}];
end

%% -------------------------------------------------
% Solve scaling system
sol = solve(invariance_eqs, ...
            [param_scaling, scaling_states], ...
            "ReturnConditions", true, ...
            "MaxDegree", 4, ...
            "IgnoreAnalyticConstraints", true, ...
            "IgnoreProperties", true);

%% -------------------------------------------------
% Parameter identifiability
for i = 1:length(param_scaling)

    u_val = sol.(char(param_scaling(i)));

    if isequal(u_val, sym(1))
        Identifiable_parameters = [Identifiable_parameters; params{i}];
    else
        unidentifiable_parameters = [unidentifiable_parameters; params{i}];
    end
end

%% -------------------------------------------------
% State observability
for i = 1:length(scaling_states)

    u_val = sol.(char(scaling_states(i)));

    if isequal(u_val, sym(1)) && E(i) == 0
        observable_states = [observable_states; states{i}];

    elseif ~isequal(u_val, sym(1))
        unobservable_states = [unobservable_states; states{i}];
    elseif E(i) == 1
        observed_state = [observed_state;states{i}];
       
    end
end
%% -------------------------------------------------
% Parameter combination identifiability
for i = 1:length(param_scaling)

    u_val_i = sol.(char(param_scaling(i)));

    for j = i+1:length(param_scaling)

        u_val_j = sol.(char(param_scaling(j)));

        m1 = simplify(u_val_i * u_val_j);
        m2 = simplify(u_val_i / u_val_j);

        if isequal(m1, sym(1)) && ...
           ~isequal(u_val_i, sym(1)) && ...
           ~isequal(u_val_j, sym(1))

            Identifiable_combinations{end+1} = ...
                [char(params{i}) '*' char(params{j})];

        elseif isequal(m2, sym(1)) && ...
               ~isequal(u_val_i, sym(1)) && ...
               ~isequal(u_val_j, sym(1))

            Identifiable_combinations{end+1} = ...
                [char(params{i}) '/' char(params{j})];

        end
    end
end

%% -------------------------------------------------
% Execution time
totaltime = toc(timeSIM);

%% -------------------------------------------------
% Return results as struct
results.ModelName                 = modelName;
results.Identifiable_parameters   = Identifiable_parameters;
results.Unidentifiable_parameters = unidentifiable_parameters;
results.Observable_states         = observable_states;
results.Unobservable_states       = unobservable_states;
results.ScalingSolution           = sol;
results.ExecutionTime             = totaltime;
results.Identifiable_combinations = Identifiable_combinations;
results.Observed_states           = observed_state;


%% -------------------------------------------------
% Display results

fprintf('\n=====================================\n');
fprintf('SIM PDE Identifiability Analysis for Model: %s\n', modelName);
fprintf('=====================================\n');

fprintf('\n >>> NOTE THAT THE SIM ALGORITHM PROVIDES A NECESSARY BUT NOT SUFFICIENT CONDITION FOR IDENTIFIABILITY AND OBSERVABILITY.');
fprintf('\n     THEREFORE, CLAIMS ABOUT IDENTIFIABILITY OR OBSERVABILITY SHOULD BE TAKEN ONLY AS POSSIBILITIES, NOT GUARANTEES.\n\n');

% Identifiable parameters
disp("Identifiable parameters:");
disp(results.Identifiable_parameters);

% Unidentifiable parameters
disp("Unidentifiable parameters:");
disp(results.Unidentifiable_parameters);

% Directly measured states
disp("Observed states (directly measured):");
disp(results.Observed_states);

% Observable states 
disp("Observable states:");
disp(results.Observable_states);

% Unobservable states
disp("Unobservable states:");
disp(results.Unobservable_states);

% Identifiable parameter combinations
disp("Identifiable combinations:");
disp(results.Identifiable_combinations);


%% -------------------------------------------------
% Display total execution time
fprintf('\nExecution time: %.4f seconds\n', results.ExecutionTime);


end
