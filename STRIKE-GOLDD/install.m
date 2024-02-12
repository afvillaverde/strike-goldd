%==========================================================================
% File that installs STRIKE-GOLDD: 
% it adds the necessary paths and creates the results folder if needed                                                       
%==========================================================================
strikegoldd_path = fileparts(mfilename('fullpath'));
addpath(strikegoldd_path)
addpath(fullfile(strikegoldd_path,'functions'))
addpath(fullfile(strikegoldd_path,'functions',filesep,'aux_AutoRepar'))
addpath(fullfile(strikegoldd_path,'functions',filesep,'aux_Lie_symmetry',filesep,'lie_sym_functions'))
addpath(fullfile(strikegoldd_path,'functions',filesep,'aux_Lie_symmetry',filesep,'lie_sym_models'))
addpath(fullfile(strikegoldd_path,'functions',filesep,'aux_Prob_Obs_Test'))
addpath(fullfile(strikegoldd_path,'models'))
addpath(fullfile(strikegoldd_path,'models',filesep,'1D_BIG_model_variants'))
if exist(fullfile(strikegoldd_path,'results'),'dir') 
    addpath(fullfile(strikegoldd_path,'results'))
else
    mkdir(fullfile(strikegoldd_path,'results'))
    addpath(fullfile(strikegoldd_path,'results'))
    fprintf('\n A results folder has been created. ');
end
fprintf('\n STRIKE-GOLDD folders have been added to the path.\n You can now use the toolbox. \n\n');
