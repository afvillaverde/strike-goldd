%==========================================================================
% StrIkE-GOLDD (Structural Identifiability taken as Extended-Generalized
% Observability with Lie Derivatives):
% a Matlab toolbox for structural identifiability and observability (SIO)
% analysis of nonlinear models.
%--------------------------------------------------------------------------
% % Version 4.2.0
% Contact: Alejandro F. Villaverde (afvillaverde@uvigo.gal)
%==========================================================================

function STRIKE_GOLDD(varargin)
clc
fprintf('\n\n -------------------------------- \n');
fprintf(' >>> STRIKE-GOLDD toolbox 4.2.0 \n');
fprintf(' -------------------------------- \n');

%==========================================================================
% Read options and add folders to path:
tStart = tic;
clearvars -global
global x f p u unidflag w wlvector 
switch nargin
    case 0
        [modelname,paths,opts,prev_ident_pars] = options;
        nmf = pwd;
    case 1
        copyfile(varargin{1},"current_options.m");
        [modelname,paths,opts,prev_ident_pars] = current_options; 
        nmf = pwd;
    case 4  % STRIKE_GOLDD.m is called by AutoRepar -> add extra paths:
        modelname       = varargin{1};
        paths           = varargin{2};
        opts            = varargin{3};
        prev_ident_pars = varargin{4};
        mf   = pwd;
        idcs = strfind(mf,filesep);
        nmf  = mf(1:idcs(end)-1);    
        paths.models    = strcat(nmf,filesep,'models');
        paths.results   = strcat(nmf,filesep,'results');
        paths.functions = strcat(nmf,filesep,'functions');
end
addpath(genpath(paths.models));
addpath(genpath(paths.results));
addpath(genpath(paths.functions));

%==========================================================================
% Read the algorithm chosen in the options file:
switch opts.algorithm
    case 1
        % Run the FISPO algorithm (remainder of this file)
    case 2
        prob_obs_test(modelname,opts,prev_ident_pars,nmf);
        return
    case 3
        ORC_DF(modelname,opts,prev_ident_pars,nmf);
        return
    case 4
        Lie_Symmetry
        return
    case 5
        if nargin < 2 
            AutoRepar
            return
        else % if nargin == 4, STRIKE-GOLDD is being called by AutoRepar => avoid recursive loop
            prob_obs_test(modelname,opts,prev_ident_pars,nmf);
            return
	    end
end

%==========================================================================
% If this line is executed, it means that the FISPO algorithm was chosen.
% Load model:

% Convert model to multi-experiment form (Optional):
if opts.multiexp == 1
    ME_analysis(modelname,opts);
    modelname=strcat(modelname,'_',num2str(opts.multiexp_numexp),'Exp');
end

load(modelname);

fprintf('\n Analyzing the %s model with the FISPO algorithm...\n',modelname);

tic
%==========================================================================
% Initialize variables:
identifiables  = [];       % identifiable parameters.
nonidentif     = [];       % unidentifiable parameters.
obs_states     = [];       % observable states.
unobs_states   = [];       % unobservable states.
obs_inputs     = [];       % observable inputs.
unobs_inputs   = [];       % unobservable inputs.
lastrank       = NaN;
unidflag       = 0;
isFISPO        = 0;
    
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
    if opts.multiexp == 1
        if opts.multiexp_user_nnzDerW
            opts.nnzDerW=opts.multiexp_nnzDerW;
        else
            opts.nnzDerW=repmat(opts.nnzDerW,1,opts.multiexp_numexp);
        end
    end
else 
    nw = 0;
    w = [];
end
if exist('u','var')
    nu = numel(u); % number of known inputs
    if opts.multiexp == 1
        if opts.multiexp_user_nnzDerU
            opts.nnzDerU=opts.multiexp_nnzDerU;
        else
            opts.nnzDerU=repmat(opts.nnzDerU,1,opts.multiexp_numexp);
        end
    end
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
meas_x = x(meas_x_indices);  %#ok<FNDSB> % names of the measured states

%==========================================================================
r  = n+q+nw;        % number of unknown variables to observe / identify
nd = ceil((r-m)/m); % lower bound of the minimum number of Lie derivatives for Oi to have full rank
fprintf('\n\n >>> Building the observability-identifiability matrix requires at least %d Lie derivatives',nd);
fprintf('\n     Calculating derivatives: ');
tic
   
%======================================================================
% Input derivatives:

%- Create array of known inputs and set certain derivatives to zero:
if numel(u)>0
    for ind_u=1:numel(u) % create array of derivatives of the inputs
        input_der(ind_u,:) = [u(ind_u),sym(strcat(char(u(ind_u)),sprintf('_d')),[1 nd])];  %[1 max(opts.nnzDerU)]
        input_der(ind_u,(opts.nnzDerU(ind_u)+2):end)=0;
    end          
else
    input_der = [];
end
syms zero_input_der_dummy_name

%- Create array of unknown inputs and set certain derivatives to zero:
if numel(w)>0
    for ind_w=1:numel(w) 
        w_der(ind_w,:) = [w(ind_w),sym(strcat(char(w(ind_w)),sprintf('_d')),[1 nd+1])];  % [1 max(opts.nnzDerW)]
        w_der(ind_w,(opts.nnzDerW(ind_w)+2):end)=0;
    end 
    wlvector     = reshape(w_der(:,1:end-1),[],1);  % reshape array as a column vector
    wlvector_dot = reshape(w_der(:,2:end),[],1);    % vector of derivatives of the unknown inputs 
    %-- Include as states only nonzero inputs/derivatives:
    [nzi,nzj,nz_wlvec] = find(wlvector); %#ok<ASGLU>
    wlvector     = nz_wlvec;   
    wlvector_dot = wlvector_dot(nzi); 
else
    wlvector      = [];
    wlvector_dot  = [];
end  

%======================================================================
% Augment state vector, dynamics:
xaug = [x;p;wlvector];                      
faug = [f;zeros(numel(p),1);wlvector_dot];   

%======================================================================
% Build Oi:
ind        = 0; % Lie derivative index (k)
lasttime   = 0;     
past_Lie   = h; 
onx        = zeros(m*(1+nd),n+q+numel(wlvector)); 
onx        = sym(onx);
onx(1:m,:) = jacobian(h,xaug); % 1st block
totaltime  = toc;       
%----------------------------------------------------------------------           
while ind < nd && lasttime < opts.maxLietime % 2nd and subsequent blocks
    tic
    ind = ind+1;
    Lieh = onx(((ind-1)*m+1):ind*m,:)*faug;
    extra_term = 0;
    if numel(u) > 0
        for j=0:ind-1
            if j < size(input_der,2) 
                lo_u_der   = input_der(:,j+1);   
                hi_u_der   = input_der(:,j+2);
                lo_u_der   = subs(lo_u_der,0,zero_input_der_dummy_name);
                extra_term = extra_term + jacobian(past_Lie,lo_u_der)*hi_u_der;
            end                        
        end
    end    
    past_Lie = Lieh + extra_term;         
    onx((ind*m+1):((ind+1)*m),:) = jacobian(past_Lie,xaug);
    lasttime = toc;
    totaltime = totaltime + lasttime;        
    fprintf('%d ',ind); % index of the last derivative to be computed
end
if ind == nd
    obsidentmatrix = sprintf('obs_ident_matrix_%s_%d_Lie_deriv',modelname,nd);
    obsidentfile = strcat(nmf,filesep,'results',filesep,obsidentmatrix);
    save(obsidentfile);
    increaseLie = 1;
    while increaseLie == 1
        fprintf('\n >>> Observability-Identifiability matrix built with %d Lie derivatives',nd);
        fprintf('\n     (calculated in %d seconds)',totaltime);
        %==========================================================================
        % Check identifiability by calculating rank:
        fprintf('\n >>> Calculating rank...');
        tic               
        rango = double(rank(onx)); 
        fprintf('\n     Rank = %d (calculated in %d seconds)',rango,toc);         
        if rango == numel(xaug)%r 
            obs_states    = x;
            obs_inputs    = w;
            identifiables = p;
            increaseLie   = 0;
        else % With that number of Lie derivatives the array is not full rank.                
            %----------------------------------------------------------
            % If there are unknown inputs, we may want to check id/obs of (x,p,w) and not of dw/dt:
            if numel(w)>0
                [identifiables,nonidentif,obs_states,unobs_states,obs_inputs,unobs_inputs] = ...
                    elim_and_recalc(unmeas_x_indices,rango,onx,opts,identifiables,obs_states,obs_inputs);   
                obs_in_no_der = intersect(w,obs_inputs);
                if ( numel(identifiables)==numel(p) && numel(obs_states)+numel(meas_x)==numel(x) && numel(obs_in_no_der)==numel(w) )
                    obs_states    = x;
                    obs_inputs    = obs_in_no_der;
                    identifiables = p;
                    increaseLie   = 0; % -> with this we skip the next 'if' block and jump to the end of the algorithm 
                    isFISPO       = 1;
                end       
            end
            %----------------------------------------------------------  
            % If possible (& necessary), calculate one more Lie derivative and retry:
            if nd < numel(xaug) && lasttime < opts.maxLietime && rango ~= lastrank && increaseLie == 1
                tic
                nd = nd+1;  
                ind = nd;  
                extra_term = 0; % reset for each new Lie derivative                    
                %-  Known input derivatives: --------------------------
                if numel(u) > 0 % Extra terms of extended Lie derivatives
                    % may have to add extra input derivatives (note that 'nd' has grown):
                    clear input_der
                    for ind_u=1:numel(u) 
                        input_der(ind_u,:) = [u(ind_u),sym(strcat(char(u(ind_u)),sprintf('_d')),[1 nd])];
                        input_der(ind_u,(opts.nnzDerU(ind_u)+2):end)=0;
                    end                        
                    for j=0:ind-1
                        if j < size(input_der,2)
                            lo_u_der   = input_der(:,j+1);
                            hi_u_der   = input_der(:,j+2);
                            lo_u_der   = subs(lo_u_der,0,zero_input_der_dummy_name);
                            extra_term = extra_term + jacobian(past_Lie,lo_u_der)*hi_u_der;
                        end 
                    end
                end
                
                %- Unknown input derivatives: ------------------------- 
                % (add new derivatives, if they are not zero:
                if numel(w)>0
                    prev_size = numel(wlvector);
                    clear w_der
                    for ind_w=1:numel(w) 
                        w_der(ind_w,:) = [w(ind_w),sym(strcat(char(w(ind_w)),sprintf('_d')),[1 nd+1])];  
                        w_der(ind_w,(opts.nnzDerW(ind_w)+2):end)=0;
                    end                         
                    wlvector     = reshape(w_der(:,1:end-1),[],1);  % reshape array as a column vector
                    wlvector_dot = reshape(w_der(:,2:end),[],1);    % vector of derivatives of the unknown inputs 
                    %-- Include as states only nonzero inputs/derivatives:
                    [nzi,nzj,nz_wlvec] = find(wlvector); %#ok<ASGLU>
                    wlvector     = nz_wlvec;   
                    wlvector_dot = wlvector_dot(nzi); 
                    %-- Augment state vector & dynamics with new input derivs:
                    xaug = [x;p;wlvector];                      
                    faug = [f;zeros(numel(p),1);wlvector_dot]; 
                    %-- Augment size of the Obs-Id matrix if needed:
                    new_size = numel(wlvector);
                    onx = [onx,zeros(ind*m,new_size-prev_size)];
                end  
                
                newLie   = onx(((ind-1)*m+1):ind*m,:)*faug; 
                past_Lie = newLie + extra_term;        
                newOnx   = jacobian(past_Lie,xaug); 
                onx      = [onx; newOnx];   
                clear newLie newOnx
                lasttime  = toc;
                totaltime = totaltime + lasttime;
                obsidentmatrix = sprintf('obs_ident_matrix_%s_%d_Lie_deriv',modelname,nd);
                obsidentfile = strcat(nmf,filesep,'results',filesep,obsidentmatrix);
                save(obsidentfile);
                lastrank = rango;                                             
            % If that is not possible, there are several possible causes:
            else 
                if nd >= numel(xaug) % The maximum number of Lie derivatives has been reached
                    unidflag = 1; 
                    fprintf('\n    The model is structurally unidentifiable as a whole');
                else
                    if rango == lastrank
                        onx = onx(1:(end-m),:);
                        nd = nd - 1;
                        unidflag = 1; % note that the pars may still be identifiable (rank deficiency may be due to initial conditions)
                    else
                        if lasttime >= opts.maxLietime
                            fprintf('\n    => More Lie derivatives would be needed to see if the model is structurally unidentifiable as a whole.');
                            fprintf('\n    However, the maximum computation time allowed for calculating each of them has been reached.');
                            fprintf('\n    You can increase it by changing <<opts.maxLietime>> (currently opts.maxLietime = %d)',opts.maxLietime);
                            unidflag = 0; 
                        end
                    end
                end  
                if isFISPO == 0
                    % Eliminate columns one by one to check identifiability of the associated parameters: 
                    [identifiables,nonidentif,obs_states,unobs_states,obs_inputs,unobs_inputs] = ...
                         elim_and_recalc(unmeas_x_indices,rango,onx,opts,identifiables,obs_states,obs_inputs);                            
                     obs_in_no_der = intersect(w,obs_inputs);
                     if ( numel(identifiables)==numel(p) && numel(obs_states)+numel(meas_x)==numel(x) && numel(obs_in_no_der)==numel(w) )
                        obs_states    = x;
                        obs_inputs    = obs_in_no_der;
                        identifiables = p; 
                        isFISPO       = 1;
                     end   
                     increaseLie = 0;
                end
            end
        end
    end        
else% If the maxLietime has been reached, but the minimum of Lie derivatives has not been calculated:
    fprintf('\n    More Lie derivatives would be needed to analyse the model.');
    fprintf('\n    However, the maximum computation time allowed for calculating each of them has been reached.');
    fprintf('\n    You can increase it by changing <<opts.maxLietime>> (currently opts.maxLietime = %d)',opts.maxLietime);
    fprintf('\n >>> Calculating rank...');
    tic
    rango = double(rank(onx));
    fprintf('\n     Rank = %d (calculated in %d seconds)',rango,toc); 
    [identifiables,nonidentif,obs_states,unobs_states,obs_inputs,unobs_inputs] = ...
    elim_and_recalc(unmeas_x_indices,rango,onx,opts,identifiables,obs_states,obs_inputs);                             
end

%==========================================================================
% Build the vectors of identifiable / non identifiable parameters and
% observable / unobservable state variables (avoid possible repetitions):
p_id         = symvar(identifiables);        
p_un         = symvar(nonidentif);          
obs_states   = symvar(obs_states);       
unobs_states = symvar(unobs_states); 
obs_inputs   = symvar(obs_inputs);       
unobs_inputs = symvar(unobs_inputs); 

%==========================================================================
% Report results:
fprintf('\n\n ------------------------ \n');
fprintf(' >>> RESULTS SUMMARY:\n');
fprintf(' ------------------------ \n');
% load(modelname)
if (numel(p_id) == numel(p)) && (numel(obs_states) == numel(x)) && (numel(obs_inputs) == numel(w))
	fprintf('\n >>> The model is Fully Input-State-Parameter Observable (FISPO):');
	if numel(w)>0, fprintf('\n     All its unknown inputs are observable.'); end
	fprintf('\n     All its states are observable.');
	fprintf('\n     All its parameters are locally structurally identifiable.');	
else
	if numel(p_id) == numel(p)
		fprintf('\n >>> The model is structurally identifiable:');
		fprintf('\n     All its parameters are structurally identifiable.');
	else    
		if unidflag == 1 
			fprintf('\n >>> The model is structurally unidentifiable.');
			fprintf('\n >>> These parameters are identifiable:\n      %s ',char(p_id));
			fprintf('\n >>> These parameters are unidentifiable:\n      %s \n',char(p_un));
		else             
			fprintf('\n >>> These parameters are identifiable:\n      %s ',char(p_id));
		end  
	end
		
	if numel(obs_states)>0,    fprintf('\n >>> These states are observable (and their initial conditions, if considered unknown, are identifiable):\n      %s ',char(obs_states)); end
	if numel(unobs_states)>0,  fprintf('\n >>> These states are unobservable (and their initial conditions, if considered unknown, are unidentifiable):\n      %s ',char(unobs_states)); end
	if numel(meas_x)>0,        fprintf('\n >>> These states are directly measured:\n      %s ',char(meas_x)); end

	if numel(obs_inputs)>0,    fprintf('\n >>> These unmeasured inputs are observable:\n      %s ',char(obs_inputs)); end
	if numel(unobs_inputs)>0,  fprintf('\n >>> These unmeasured inputs are unobservable:\n      %s ',char(unobs_inputs)); end
	if numel(u)>0,             fprintf('\n >>> These inputs are known:\n      %s ',char(u)); end
end

%==========================================================================
totaltime = toc(tStart);
fprintf('\n Total execution time: %d \n\n',totaltime);

%==========================================================================
% Save results:
resultsname = sprintf('id_results_%s_%s',modelname,date);   
fullresultsname = strcat(nmf,filesep,'results',filesep,resultsname);
save(fullresultsname);

%==========================================================================
% Delete auxiliary files:
if exist("current_options.m",'file')
    delete("current_options.m")
end

end    
