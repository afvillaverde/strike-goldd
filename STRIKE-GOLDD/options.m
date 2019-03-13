%==========================================================================
% THE USER CAN DEFINE THE PROBLEM AND SET OPTIONS IN THE FOLLOWING LINES:                                                       
%==========================================================================

function [modelname,paths,opts,submodels,prev_ident_pars] = options()

%%% (1) MODEL: 
modelname = 'C2M'; 

%%% (2) PATHS:
paths.meigo     = '/home/alexandre/MEIGO_AMICI_Nov2016/MEIGO_Nov2016/MEIGO';      
paths.models    = strcat(pwd,filesep,'models');
paths.results   = strcat(pwd,filesep,'results');
paths.functions = strcat(pwd,filesep,'functions');
                            
%%% (3) IDENTIFIABILITY OPTIONS:
opts.numeric    = 0;       % calculate rank numerically (= 1) or symbolically (= 0)
opts.replaceICs = 0;       % replace states with known initial conditions (= 1) or use generic values (= 0) when calculating rank
opts.checkObser = 1;       % check identifiability of initial conditions (1 = yes; 0 = no).
opts.checkObsIn = 1;       % check input observability (1 = yes; 0 = no).
opts.findcombos = 0;       % try to find identifiable combinations? (1 = yes; 0 = no).
opts.unidentif  = 0;       % use method to try to establish unidentifiability instead of identifiability, when using decomposition. 
opts.forcedecomp= 0;       % always decompose model (1 = yes; 0 = no).
opts.decomp     = 0;       % decompose model if the whole model is too large (1 = yes; 0 = no: instead, calculate rank with few Lie derivatives).
opts.decomp_user= 0;       % when decomposing model, use submodels specified by the user (= 1) or found by optimization (= 0). 
opts.maxLietime = 1000;    % max. time allowed for calculating 1 Lie derivative.
opts.maxOpttime = 300;     % max. time allowed for every optimization (if optimization-based decomposition is used).
opts.maxstates  = 6;       % max. number of states in the submodels (if optimization-based decomposition is used).
opts.nnzDerU    = 1;       % numbers of nonzero derivatives of the measured inputs (u); may be 'inf'
opts.nnzDerW    = 1;       % numbers of nonzero derivatives of the unmeasured inputs (w); may be 'inf'
opts.nnzDerIn   = opts.nnzDerU; % deprecated option

%%% (4) User-defined submodels for decomposition (may be left = []):  
submodels = []; 
%%% Submodels are specified as a vector of states, as e.g.:
%         submodels{1}  = [1 2];
%         submodels{2}  = [1 3];
        
%%% (5) Parameters already classified as identifiable may be entered below.
prev_ident_pars = [];
%%% They must first be declared as symbolic variables. For example:
%         syms lambda rho c
%         prev_ident_pars = [lambda rho c];
end