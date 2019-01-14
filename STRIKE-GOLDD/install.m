strikegoldd_path = fileparts(mfilename('fullpath'));
addpath(strikegoldd_path)
addpath(fullfile(strikegoldd_path,'functions'))
addpath(fullfile(strikegoldd_path,'models'))
addpath(fullfile(strikegoldd_path,'models',filesep,'1D_BIG_model_variants'))
addpath(fullfile(strikegoldd_path,'results'))
fprintf('\n StrIkE-GOLDD folders added to the path \n\n');