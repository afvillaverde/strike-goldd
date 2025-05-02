% sg2strucID.m
% This script converts .mat model files, defined in the STRIKE-GOLDD format, 
% to structured .txt files that can be analysed by StrucID
% Authored by Mahmoud Shams Falavarjani, initial version April 2025.

function sg2strucID
    % Create the GUI figure
    figure('Name', 'MAT to TXT Converter', 'Position', [500, 400, 350, 150]);
    uicontrol('Style', 'pushbutton', 'String', 'Select .mat File', ...
              'Position', [80 60 200 40], 'Callback', @convertMatToTxt);
end

function convertMatToTxt(~, ~)
    [file,path] = uigetfile('*.mat','select mat file');
    if isequal(file,0)
        return;
    end

    % Load symbolic variables from the file
    s = load(fullfile(path, file));

    % Extract model components
    x = s.x;
    p = s.p;
    h = s.h;
    f = s.f;
    u = s.u;

    % Determine input variable names
    if isa(u, 'sym')
        uVars = symvar(u);
    else
        uVars = sym([]);
    end

    % find input u0 and remove it from equations(for a model that has u0)
    uVarsFiltered = uVars(~strcmp(string(uVars), 'u0'));
    if any(strcmp(string(symvar(f)), 'u0'))
        u0 = sym('u0');
        f = subs(f, u0, 0);
    end

    % Identify algebraic variables not in x, p, or u 
    allKnownVars = [x(:); p(:); uVarsFiltered(:)];
    allFVars = unique(symvar(f));
    algebraicVars = setdiff(allFVars, allKnownVars);

    % Start writing strings
    txt = {};

    % Write 'Algebraic Rules!' 
    txt{end+1} = 'Algebraic Rules!';
    for i = 1:length(algebraicVars)
        val = rand()*5;  % Assign a random value <= 5
        txt{end+1} = sprintf('%s = %.3f', char(algebraicVars(i)), val);
    end
    txt{end+1} = '';

    % Write 'ODEs (define the individual ODE equations - 1 per line)!
    txt{end+1} = 'ODEs (define the individual ODE equations - 1 per line)!';
    for i = 1:length(f)
        eq = f(i);
        eqStr = char(eq);
        txt{end+1} = sprintf('d%s/dt = %s;', char(x(i)), eqStr);
    end
    txt{end+1} = '';

    % Write 'Input variables! 
    txt{end+1} = 'Input variables!';
    for i = 1:length(uVarsFiltered)
        txt{end+1} = sprintf('%s = ', char(uVarsFiltered(i)));
    end
    txt{end+1} = '';

    % Write 'Measured Outputs (define the measured sensors - 1 per line)!
    txt{end+1} = 'Measured Outputs (define the measured sensors - 1 per line)!';
    for i = 1:length(h)
            txt{end+1} = sprintf('y%d = %s', i, char(h(i)));
    end
    txt{end+1} = '';

    % Write 'Parameter names and values (define all the system parameters - 1 per line, OPTIONAL - define known parameter values)!'
    txt{end+1} = 'Parameter names and values (define all the system parameters - 1 per line, OPTIONAL - define known paramter values)!';
    for i = 1:length(p)
        txt{end+1} = sprintf('%s = ', char(p(i)));
    end
    txt{end+1} = '';

    % Write 'State names and initial values (define all the model state names - 1 per line, OPTIONAL - define known initial values)!'
    txt{end+1} = 'State names and initial values (define all the model state names - 1 per line, OPTIONAL - define known initial values)!';
    for i = 1:length(x)
         txt{end+1} = sprintf('%s = ', char(x(i)));
    end
    txt{end+1} = '';

    % Write 'Analyse (list the unknown parameter and initial conditions which should be included into the structural identifiability analysis)!'
    txt{end+1} = 'Analyse (list the unknown parameter and initial conditions which should be included into the structural identifiability analysis)!';

    % Save to txt file
    [~, name] = fileparts(file);
    txtFile = fullfile(path, [name, '.txt']);
    fid = fopen(txtFile, 'w');
    for i = 1:length(txt)
        fprintf(fid, '%s\n', txt{i});
    end
    fclose(fid);

    msgbox(['TXT file saved to ', txtFile], 'Success');
end
