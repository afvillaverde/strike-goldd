%--------------------------------------------------------------------------
% Function that analyses the identifiability and observability of
% nonlinear models with an implementation of the method presented in: 
% Sedoglavic, A. "A probabilistic algorithm to test local algebraic 
% observability in polynomial time." J Symbol Comput 33.5 (2002): 735-755.
%--------------------------------------------------------------------------

function prob_obs_test(modelname,opts,prev_ident_pars,nmf)

%==========================================================================
% Load model:

% Convert model to multi-experiment form (Optional):
if opts.multiexp == 1
    ME_analysis(modelname,opts);
    modelname=strcat(modelname,'_',num2str(opts.multiexp_numexp),'Exp');
end

load(modelname)  %#ok<LOAD> 
fprintf(['\n Analyzing the %s model with the probabilistic' ...
    ' algorithm... \n'], modelname)

%==========================================================================
% Remove parameters that have already been classified as identifiable:
if exist('prev_ident_pars','var') == 1
    for np=1:numel(prev_ident_pars)
        [~, original_index] = ismember(prev_ident_pars(np),p); 
        p(original_index)=[];                %#ok<AGROW> 
    end
end

%==========================================================================
% Dimensions of the problem:
m    = numel(h);                  %#ok<NODEF> % number of outputs
n    = numel(x);                  %#ok<NODEF> % number of states
q    = numel(p);                  % number of unknown parameters
if size(h,2)>size(h,1),h=h.';end
if size(x,2)>size(x,1),x=x.';end
if size(p,2)>size(p,1),p=p.';end
if exist('w','var')
    nw = numel(w); %#ok<NODEF> % number of unknown inputs
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
if size(w,2)>size(w,1),w=w.';end
if exist('u','var')
    nu = numel(u); %#ok<NODEF> % number of known inputs
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
if size(u,2)>size(u,1),u=u.';end
fprintf('\n >>> The model contains:\n %d states:\n %s',n,char(x));
fprintf('\n %d outputs:\n %s',m,char(h));
fprintf('\n %d known inputs:\n %s',nu,char(u));
fprintf('\n %d unknown inputs:\n %s',nw,char(w));
fprintf('\n %d parameters:\n %s',q,char(p));

%==========================================================================
% Create augmented state vector with unknown inputs and derivatives 
nwl=0;
if nw~=0
    if length(opts.nnzDerW) == 1
        if opts.nnzDerW == inf
             fprintf('\n\n ------------------------ \n');
            fprintf(' >>> WARNING:\n');
            fprintf(' ------------------------ \n');
            fprintf([' Not able to do calculations with infinite ' ...
                'number of unknown input derivatives:\n'])
            fprintf(['     Number of derivatives for unknown input ' ...
                'fixed to 3.\n'])
            fprintf(['     To increse number of unknown input ' ...
                'derivatives change number on file options.\n'])
            fprintf(['     To try the analysis with an infinite number' ...
                ' of unknown input derivatives use FISPO algorithm.\n'])
            opts.nnzDerW = 3;
            w_der=sym('w_der',[nw*opts.nnzDerW,1]);
            for ind_w=1:nw 
                w_der(ind_w:nw:end-nw+ind_w) = ...
                    sym(strcat(char(w(ind_w)),sprintf('_d')),[1 opts.nnzDerW]);  
            end 
            wlvector     = [w;w_der];  % reshape array as a column vector
            wlvector_dot = [w_der;zeros(nw,1)];    % vector of derivatives
                                                   % of the unknown inputs
        elseif opts.nnzDerW == 0
            wlvector=w;
            wlvector_dot = zeros(nw,1);
        else
            w_der=sym('w_der',[nw*opts.nnzDerW,1]);
            for ind_w=1:nw 
                w_der(ind_w:nw:end-nw+ind_w) = ...
                    sym(strcat(char(w(ind_w)),sprintf('_d')),[1 opts.nnzDerW]);  
            end 
            wlvector     = [w;w_der];  % reshape array as a column vector
            wlvector_dot = [w_der;zeros(nw,1)];    % vector of derivatives
                                                   % of the unknown inputs 
        end
    else
        for ind_w=1:nw
            if opts.nnzDerW(ind_w) == inf
                fprintf('\n\n ------------------------ \n');
                fprintf(' >>> WARNING:\n');
                fprintf(' ------------------------ \n');
                fprintf([' Not able to do calculations with infinite' ...
                    ' number of unknown input derivatives:\n'])
                fprintf(['     Number of derivatives for unknown input' ...
                    ' fixed to 3.\n'])
                fprintf(['     To increse number of unknown input' ...
                    ' derivatives change number on file options.\n'])
                fprintf(['     To try the analysis with an infinite' ...
                    ' number of unknown input derivatives use FISPO' ...
                    ' algorithm.\n'])
                opts.nnzDerW(ind_w)=3;
            end 
        end
        wlvector=sym('w_der',[nw+sum(opts.nnzDerW(1:nw)),1]);
        wlvector_dot=sym('wlvector_dot',[nw+sum(opts.nnzDerW(1:nw)),1]);
        ind=1;
        for ind_w=1:nw  
            wlvector(ind)=w(ind_w);
            wlvector(ind+1:ind+opts.nnzDerW(ind_w)) = ...
                sym(strcat(char(w(ind_w)),sprintf('_d')), ...
                [1 opts.nnzDerW(ind_w)]); 
            wlvector_dot(ind:ind+opts.nnzDerW(ind_w)-1)=...
                wlvector(ind+1:ind+opts.nnzDerW(ind_w));
            wlvector_dot(opts.nnzDerW(ind_w)+ind)=0;
            ind=ind+opts.nnzDerW(ind_w)+1;
        end 
    end
    nwl=numel(wlvector);
    % Construct augmented state vector and state function by taking unknown
    % inputs ass state variables
    x = [x;wlvector];
    f = [f;wlvector_dot]; %#ok<NODEF> 
    n = numel(x);
end



%==========================================================================
% Assignment of the prime numer used by the method:
Myprime = nextprime(10^6);
fprintf('\n >>> Computations are done modulo: %d \n',Myprime);

%==========================================================================
% Build Oi, the observability-identifiability matrix:
tic
[onx]=build_OI_sed(Myprime,x,p,u,h,f,n,q,nu,m,opts);
timematrix=toc;

% Report results:
fprintf('\n\n ------------------------ \n');
fprintf(' >>> RESULTS SUMMARY:\n');
fprintf(' ------------------------ \n');
fprintf(['\n Observability-Identifiability matrix calculated in %d ' ...
    'seconds. \n\n'],timematrix);

%==========================================================================
% Analysis of identifiability and observability:
tic
[p_id,  p_un, obs_states, unobs_states, obs_inputs, unobs_inputs,meas_x,rango]...
    = ObservabilityAnalysis(x,p,onx,Myprime,n,q,nwl,h,u); %#ok<ASGLU> 
timeAnalysis=toc;

fprintf('\n >>> Matrix analyzed in %d seconds. \n',timeAnalysis);

totaltime=timematrix+timeAnalysis;
fprintf('\n Total execution time: %d seconds. \n',totaltime);

xaug=[x;p]; %#ok<NASGU> 

resultsname = sprintf('id_results_%s_%s',modelname,date);
fullresultsname = strcat(nmf,filesep,'results',filesep,resultsname);
save(fullresultsname);

warning('on','all')
