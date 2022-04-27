%==========================================================================
% THE USER MUST DEFINE THE PROBLEM AND SET OPTIONS IN THE FOLLOWING LINES:                                                       
%==========================================================================

function [modelname,paths,opts,submodels,prev_ident_pars] = options()

%=================== BASIC OPTIONS ("CHOOSE WHAT TO DO") ==================
%%% (1) CHOOSE MODEL TO ANALYSE: 
modelname ='C2M';    % Name of a .mat file placed in the 'models' folder 
  
%%% (2) CHOOSE TYPE OF ANALYSIS:
opts.algorithm = 2;  % Choose one of the following:
   %--- Structural Identifiability and Observability (SIO) algorithms:
      % 1: FISPO (default) -- applicable to nonlinear models in general
      % 2: Prob_Obs_Test   -- applicable to rational models
      % 3: ORC-DF          -- applicable to systems affine in the inputs
   %--- Symmetries and Reparameterization algorithms:
      % 4: Lie_Symmetry    -- search for symmetries
      % 5: AutoRepar       -- automatic reparameterization

%================== ADDITIONAL OPTIONS FOR SIO ANALYSIS ===================     
%%% (3) MAIN STRUCTURAL IDENTIFIABILITY & OBSERVABILITY (SIO) OPTIONS:
opts.maxLietime = 100;     % In FISPO, max. time allowed for calculating 1 Lie derivative.
opts.nnzDerU    = inf;     % In FISPO, numbers of nonzero derivatives of the measured inputs (u); may be 'inf'
opts.nnzDerW    = 1;       % In FISPO and ORC-DF, numbers of nonzero derivatives of the unmeasured inputs (w); may be 'inf'
opts.numeric    = 0;       % In FISPO and ORC-DF, calculate rank numerically (= 1) or symbolically (= 0)
opts.replaceICs = 0;       % In FISPO and ORC-DF, replace states with specific initial conditions (= 1) or use generic values (= 0) when calculating rank
prev_ident_pars = [];      % parameters assumed known, or already classified as identifiable. Example:
% syms p2 p5
% prev_ident_pars = [p2 p5];

%%% (4) OPTIONS FOR AFFINE SYSTEMS (apply only to the ORC-DF algorithm):
opts.affine_tStage       = 1000;  % max. computation time for the last iteration.
opts.affine_kmax         = 4;     % max. number of iterations.
opts.affine_parallel     = 0;     % use parallel toolbox (=1) or not (=0) to calculate partial ranks.
opts.affine_workers      = 4;     % number of workers for parallel pool.
opts.affine_graphics     = 1;     % display graphics (=1) or not (=0)
opts.affine_delete_model = 1;     % delete affine model when finished (=1) or not (=0).

%%% (5) DECOMPOSITION OPTIONS (apply only to the FISPO algorithm):
opts.forcedecomp= 0;       % always decompose model (1 = yes; 0 = no).
opts.decomp     = 0;       % decompose model if the whole model is too large (1 = yes; 0 = no: instead, calculate rank with few Lie derivatives).
opts.decomp_user= 0;       % when decomposing model, use submodels specified by the user (= 1) or found by optimization (= 0). 
opts.maxOpttime = 30;      % max. time allowed for every optimization (if optimization-based decomposition is used).
opts.maxstates  = 6;       % max. number of states in the submodels (if optimization-based decomposition is used).
submodels       = [];      % User-defined submodels for decomposition (may be left = [], or specified as a vector of states, as e.g.:)
% submodels{1}  = [1 2];
% submodels{2}  = [1 3];

%%% (6) MULTI-EXPERIMENT OPTIONS (apply to the 3 SIO algorithms):
opts.multiexp              = 0;     % Execute multi-experiment analysis (=1) or not (=0).
opts.multiexp_numexp       = 2;     % Number of experiments.
opts.multiexp_user_nnzDerU = 0;     % Set manually the number of non-zero known input derivatives in each experiment (=1) or not (=0).
opts.multiexp_nnzDerU      = [0 1]; % Number of non-zero known input derivatives in each experiment (Rows=inputs;Columns=experiments).
opts.multiexp_user_nnzDerW = 0;     % Set manually the number of non-zero unknown input derivatives in each experiment (=1) or not (=0).
opts.multiexp_nnzDerW      = [1 0]; % Number of non-zero unknown input derivatives in each experiment (Rows=inputs;Columns=experiments).
opts.multiexp_user_ics     = 0;     % Set manually the initial conditions for each experiment (=1) or not (=0).
%- Multi-experiment initial conditions (Rows=variables;Columns:experiments):
opts.multiexp_ics       = [ [1,0,0,1,0,1,0,0,0].', [1,0,0,1,0,1,0,0,0].' ];
%- Which initial conditions are replaced (Rows=variables;Columns=experiments):
opts.multiexp_known_ics = [ [0,1,1,0,1,0,1,1,1].', [0,1,1,0,1,0,1,1,1].' ];

%=========== LIE SYMMETRIES & REPARAMETERIZATION OPTIONS: =================  
opts.ansatz     = 2; % Type of Ansatz: %  uni -> Univariate (1)
                                       %  par -> Partially variate (2)
                                       %  multi -> Multivariate (3)
opts.degree     = 2; % Degree of Ansatz Polynomial
opts.tmax       = 4; % Maximum degree of Lie series
opts.ode_n      = 1; % (Only in Matlab R2020a and later:) use ode solver (=1) or not (=0)
opts.use_existing_results = 0; % if the model has already been analysed with STRIKE-GOLDD (=1) or not (=0)
opts.results_file = 'id_results_1D_BIG_p_16-Oct-2020.mat'; % .mat file to use if use_existing_results = 1                              

%================================= PATHS ==================================
paths.meigo     = '/.../MEIGO';      
paths.models    = strcat(pwd,filesep,'models');
paths.results   = strcat(pwd,filesep,'results');
paths.functions = strcat(pwd,filesep,'functions');

%======= DEPRECATED OPTIONS (FOR COMPATIBILITY WITH OLDER VERSIONS) =======
opts.unidentif  = 0;            % use method to try to establish unidentifiability instead of identifiability when using decomposition. 
opts.nnzDerIn   = opts.nnzDerU; % deprecated option
opts.checkObser = 1;            % check state observability, i.e. identifiability of initial conditions (1 = yes; 0 = no).
opts.checkObsIn = 1;            % check input observability (1 = yes; 0 = no).
opts.affine     = 0;            % use the ORC-DF algorithm for affine control systems (=1) or not(=0).

end
