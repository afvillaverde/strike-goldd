% Main file for analysing NonLinear controllability & accessibility

clear; clc;
addpath(genpath('models'))
addpath(genpath('functions'))
% ======================== Explanation of variables =======================
% 
% ctrl_analysis_MAIN([opts.LC,opts.ARC,opts.LARC,opts.GSC],'model_name',...
%                     opts.numericLC, opts.maxtime, x0(optional));
%     opts.LC        --> 1 if you want to check the Linearization Condition
%                        0 if not
%     opts.ARC       --> 1 if you want to check the Accesibility Rank 
%                          Condition
%                        0 if not
%     opts.LARC      --> 1 if you want to check the Lie Algebraic Rank 
%                          Condition
%                        0 if not
%     opts.GSC       --> 1 if you want to check Sussmann's General 
%                          Sufficient Condition
%                        0 if not
%     'model_name'   --> Name of a .mat file placed in the 'models' folder
%     opts.numericLC --> 0 if you want to check the Linearization
%                          Condition symbolically
%                        1 for numeric computation
%     opts.maxtime   --> max time for the each test in seconds
%     x0             --> Specific initial point. If no  point is given it 
%                        would try to compute equilibrium points.

% ================================ Example ================================
%
% You can rewrite the following line to run the analisys on your model,
% more examples are given in the script 'ctrl_reach_unit_tests_script.m'
% found in: functions/aux_NLcontrollability
ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_53', 0, 10); 

