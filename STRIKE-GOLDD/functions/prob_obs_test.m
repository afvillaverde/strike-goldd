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
fprintf('\n Analyzing the %s model with the probabilistic algorithm... \n',modelname)

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
% Assignment of the prime numer used by the method:
Myprime = nextprime(10^6);
fprintf('\n >>> Computations are done modulo: %d \n',Myprime);

%==========================================================================
% Build Oi, the observability-identifiability matrix:
tic
[onx]=build_OI_sed(Myprime,x,p,u,w,h,f,n,q,nu,nw,m,opts);
timematrix=toc;
fprintf('\n >>> Observability-Identifiability matrix calculated in %d seconds. \n',timematrix);

%==========================================================================
    % Analysis of identifiability and observability:
tic
 [observable,unobservable] = ObservabilityAnalysis(x,p,w,onx,Myprime,n,q,nw);
timeAnalysis=toc;
fprintf('\n >>> Matrix analyzed in %d seconds. \n',timeAnalysis);

totaltime=timematrix+timeAnalysis;
fprintf('\n >>> Total time: %d seconds. \n',totaltime);

resultsname = sprintf('id_results_%s_%s',modelname,date);
fullresultsname = strcat(nmf,filesep,'results',filesep,resultsname);
save(fullresultsname);
