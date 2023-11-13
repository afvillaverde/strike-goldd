% Test suite for the controllability/accessibility analysis scripts.
% It analyses a battery of models. Uncomment the lines you want to run.
% Usage: 
%     tests           = varargin{1}; % vector (0/1) that selects the tests:
%       --> opts.LC   = tests(1); 
%       --> opts.ARC  = tests(2); 
%       --> opts.LARC = tests(3); 
%       --> opts.GSC  = tests(4);
%     modelname       = varargin{2}; % .mat file in the 'models' folder 
%     opts.numericLC  = varargin{3}; % 0 or 1
%     opts.maxtime    = varargin{4}; % max time for the tests in seconds
%     x0              = varargin{5}; % Specific initial point. If no
                                     % point is given it would try to 
                                     % compute equilibrium points
%==========================================================================
clear; clc;
addpath('models')

% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_11',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_31',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_31_odd',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_31_even',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_32',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_33',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_41',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_51',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_52',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_53',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_61',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_61b',0,10);


% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_30_VIDYASAGAR',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_35_VIDYASAGAR',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_45_VIDYASAGAR',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_Sussmann87',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_Sussmann83',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_1_Willigenburg',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_2_Willigenburg',0,10);
% ctrl_analysis_MAIN([1,1,1,1],'ctrl_ex_Lenhart_1',0,10);


% ctrl_analysis_MAIN([1,1,1,1],'1D_BIG',1,100,[1;1;1]); 
% ctrl_analysis_MAIN([1,1,1,1],'Arabidopsis',0,600);
% ctrl_analysis_MAIN([1,1,1,1],'BachmannJAKSTAT',1,1000,ones(25,1));
% ctrl_analysis_MAIN([1,1,1,1],'Bolie',1,600); 
% ctrl_analysis_MAIN([1,1,1,1],'C2M',1,600); 
% ctrl_analysis_MAIN([1,1,1,1],'CSTR_observer',1,100);
% ctrl_analysis_MAIN([1,1,1,1],'EGF_Brown_2inputs',1,600,ones(26,1));
% ctrl_analysis_MAIN([1,1,1,1],'Fujita',1,600,ones(9,1));
% ctrl_analysis_MAIN([1,1,1,1],'HIV_known_u',1,1000); 
% ctrl_analysis_MAIN([1,1,1,1],'microbial_community',1,600,ones(5,1));
% ctrl_analysis_MAIN([1,1,1,1],'NFKB',1,1000,ones(15,1));
% ctrl_analysis_MAIN([1,1,1,1],'NFKB',1,1000);
% ctrl_analysis_MAIN([1,1,1,1],'NFKB_Merkt_affine_u',1,1000,ones(10,1));
% ctrl_analysis_MAIN([1,1,1,1],'PK_known_input',1,100);
% ctrl_analysis_MAIN([1,1,1,1],'RaiaJAKSTAT',1,600,ones(10,1)); 
ctrl_analysis_MAIN([1,1,1,1],'Vajda1989_u',1,600);
