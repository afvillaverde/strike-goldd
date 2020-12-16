%=========================================================================%
%====ORC-DF: An implementation of the Observability Rank Condition for====%
%======================systems with Direct Feedthrough====================%

function ORC_DF(modelname,opts)

tStart=tic;

%============================Load model===================================%

%(Optional) Create multi-experiment model:
if opts.multiexp == 1
    ME_analysis(modelname,opts);
    modelname=strcat(modelname,'_',num2str(opts.multiexp_numexp),'Exp');
end

affine_modelname = sprintf('affine_%s',modelname);
path=strcat(pwd,filesep,'models',filesep,strcat(affine_modelname,'.mat'));
exist_affine_model=exist(path,'file');

%If the model has already been converted into affine in inputs form:
if exist_affine_model==2
    modelname=affine_modelname;
end

load(modelname) %#ok<*LOAD>

fprintf('\n >>> Analyzing observability of %s with affine in control system algorithm\n',modelname);

%======================Dimensions of the problem==========================%

%Number of states:
ns=numel(x);   %#ok<*USENS>
%Number of outputs:
m=numel(h);
%Number of unknown parameters:
if exist('p','var')
    np=numel(p);
else
    p=[];
    np=0;
end
%Number of unknown inputs:
if exist('w','var')&&numel(w)>0 %#ok<*NODEF>
    nw=numel(w);
    if opts.multiexp == 1
        opts.nnzDerW=repmat(opts.nnzDerW,1,opts.multiexp_numexp);
    end
    if exist_affine_model~=2
    %Chek that system is affine in unknown inputs:    
        syms coeff_w coefh_w
        for j=1:nw
            for i=1:ns
                try
                    coeff_w=coeffs(f(i),w(j),'all');
                    if numel(coeff_w)>2
                        warning('Unable to convert model to affine in inputs form. System dynamics is not affine in unknown inputs.')
                        return
                    end
                catch
                    warning('Unable to convert model to affine in inputs form. System dynamics is not affine in unknown inputs.')
                    return
                end
                variables=symvar(coeff_w);
                nz_variables=simplify(variables-w(j)*ones(1,numel(variables)));
                if numel(find(nz_variables==0))~=0
                    warning('Unable to convert system to affine in inputs form. System dynamics is not affine in unknown inputs.')
                    return
                end
            end
            for i=1:m
                try 
                    coefh_w=coeffs(h(i),w(j),'all');
                    if numel(coefh_w)>2
                        warning('Unable to convert model to affine in inputs form. Output dynamics is not affine in unknown inputs.')
                        return
                    end
                catch
                    warning('Unable to convert model to affine in inputs form. Output dynamics is not affine in unknown inputs.')
                    return
                end
                variables=symvar(coefh_w);
                nz_variables=simplify(variables-w(j)*ones(1,numel(variables)));
                if numel(find(nz_variables==0))~=0
                    warning('Unable to convert system to affine in inputs form. Output dynamics is not affine in unknown inputs.')
                    return
                end
            end
        end
    end
     
    %Create array of unknown inputs derivatives and set certain derivatives to zero:
    w_der=sym(zeros(nw,opts.affine_kmax+1));
    if numel(opts.nnzDerW)==nw
      for ind_w=1:nw
          if opts.nnzDerW(ind_w)<opts.affine_kmax+1
              w_der(ind_w,1:opts.nnzDerW(ind_w))=sym(strcat(char(w(ind_w)),sprintf('d')),[1 opts.nnzDerW(ind_w)]);
          else
              w_der(ind_w,:)=sym(strcat(char(w(ind_w)),sprintf('d')),[1 opts.affine_kmax+1]);
          end
      end   
    elseif numel(opts.nnzDerW)>0
      warning('The size of vector that contains the number of nonzero unknown inputs derivatives must be %d',nw)
      return
    else
    %if opts.nnzDerw=[],it is assumed that opts.nnzDerw(i)=inf, i=1,...,nw:
      for ind_w=1:nw, w_der(ind_w,:)=sym(strcat(char(w(ind_w)),sprintf('d_')),[1 opts.affine_kmax+1]);end
    end   
else
    w=[];
    nw=0;
    w_der=sym(zeros(nw,opts.affine_kmax+1));
end

%Number of known inputs:
if exist('u','var')&& numel(u)>0
    nu=numel(u);
    %Convert system into affine in inputs form:
    if exist_affine_model~=2
        fprintf('\n >>> Converting control system into affine in inputs form...\n');
        syms coeff_u coefh_u
        %Initialize known inputs coefficients of system dynamics: 
        fu=sym(zeros(ns,nu));
        %Initialize known inputs coefficients of outputs dynamics:
        hu=sym(zeros(m,nu));
        %Chek that system is affine in known inputs and almacenate coefficients:
        for j=1:nu
            for i=1:ns
                try
                    coeff_u=coeffs(f(i),u(j),'all');
                    if numel(coeff_u)>2
                        warning('Unable to convert system to affine in inputs form. System dynamics is not affine in known inputs.')
                        return
                    end
                catch
                    warning('Unable to convert system to affine in inputs form. System dynamics is not affine in known inputs.')
                    return
                end
                variables=symvar(coeff_u);
                nz_variables=simplify(variables-u(j)*ones(1,numel(variables)));
                if numel(find(nz_variables==0))~=0
                    warning('Unable to convert system to affine in inputs form. System dynamics is not affine in known inputs.')
                    return
                end
                %Almacenate uk coefficients of system dynamics: 
                if numel(coeff_u)==2, fu(i,j)=coeff_u(1);end
            end
            for i=1:m
                try
                    coefh_u=coeffs(h(i),u(j),'all');
                    if numel(coefh_u)>2
                        warning('Unable to convert system to affine in inputs form. Output dynamics is not affine in known inputs.')
                        return
                    end
                catch
                    warning('Unable to convert system to affine in inputs form. Output dynamics is not affine in known inputs.')
                    return
                end
                variables=symvar(coefh_u);
                nz_variables=simplify(variables-u(j)*ones(1,numel(variables)));
                if numel(find(nz_variables==0))~=0
                    warning('Unable to convert system to affine in inputs form. Output dynamics is not affine in known inputs.')
                    return
                end
                %Almacenate uk coefficients of output dynamics:
                if numel(coefh_u)==2, hu(i,j)=coefh_u(1);end
            end
        end
        %State-unknown inputs 0-augmented system dynamics:
        fxw=simplify(f-fu*u);
        %State-unknown inputs output dynamics:
        hxw=simplify(h-hu*u);
        %Remove dependent columns of fu, hu:
        zc_fu=[];
        zc_hu=[];
        xaug=[x;u;w;f;h];
        allvariables = symvar(xaug);
        numbers = vpa(0.1+rand(size(allvariables)));
        num_fu=subs(fu,allvariables,numbers);
        num_hu=subs(hu,allvariables,numbers);
        rank_fu=rank(num_fu);
        rank_hu=rank(num_hu);
        for j=1:nu
            aux_fu=num_fu;
            aux_hu=num_hu;
            aux_fu(:,j)=[];
            aux_hu(:,j)=[];
            if rank(aux_fu)==rank_fu
                zc_fu=[zc_fu j];
            end
            if rank(aux_hu)==rank_hu
                zc_hu=[zc_hu j];
            end
        end
        fu(:,zc_fu)=[];
        hu(:,zc_hu)=[];
        % Save affine model:
        fullaffinename = strcat(pwd,filesep,'models',filesep,affine_modelname);
        if exist('ics','var')==0,ics=[];end
        if exist('known_ics','var')==0,known_ics=[];end
        save(fullaffinename,'x','u','w','p','f','h','fu','fxw','hu','hxw','ics','known_ics');
    end
else
    u=[];
    nu=0;
    fu=[];
    hu=[];
    fxw=f;
    hxw=h;
end

%Remove parameters that have already been classified as identificable:
if exist('prev_ident_pars','var')
    for i=1:numel(prev_ident_pars)
        [~, original_index] = ismember(prev_ident_pars(i),p);
        p(original_index)=[]; %#ok<*AGROW>
    end
    np=numel(p);
end

fprintf('\n >>> The model contains:\n %d states:\n %s',ns,char(x));
fprintf('\n %d outputs:\n %s',m,char(h));
fprintf('\n %d known inputs:\n %s',nu,char(u));
fprintf('\n %d unknown inputs:\n %s',nw,char(w));
fprintf('\n %d parameters:\n %s',np,char(p));

%====================(Optional)Parallel preferences=======================%
if opts.affine_parallel_rank == 1 || opts.affine_parallel_Lie == 1
    fprintf('\n >>> Initializating parallel preferences...')
    pool=gcp('nocreate');
    mycluster=parcluster('local');
    max_workers=mycluster.NumWorkers;
    if isempty(pool)
        if opts.affine_workers<=max_workers
            parpool(opts.affine_workers);
        else
            warning('The maximum number of workers is %d. Introduce a number of workers up to %d.',max_workers,max_workers)
            return
        end
    elseif pool.NumWorkers~=opts.affine_workers
        if opts.affine_workers<=max_workers
            delete(pool);
            parpool(opts.affine_workers);
        else
            warning('The maximum number of workers is %d. Introduce a number of workers up to %d.',max_workers,max_workers)
            return
        end
    end
end

%%===========================Initialization==============================%%

%Initialize computation time of each stage:
stage_time=zeros(opts.affine_kmax+1,1);
%Stage_time can be separated in:
%Time spent calculating the differential of Omega:
dif_time=zeros(opts.affine_kmax+1,1);
%Time spent calculating the rank of Omega:
rank_time=zeros(opts.affine_kmax+1,1);
%Time spent calculating partial ranks of Omega:
partial_rank_time=zeros(opts.affine_kmax+1,1);
%Initialize number of states in each stage:
nxau=zeros(1,opts.affine_kmax+1);
%Number of states of 0-augmented system:
nxau(1)=ns+np+nw;
%Maximum number of states:
max_ns=nxau(1)+numel(find(w_der));
%Initialize rank of Omega for each stage:
rank_dif_omega=zeros(1,opts.affine_kmax);
%Initialize matrix whose components are the values of partial ranks for each state and iteration:
partial_ranks=zeros(max_ns,opts.affine_kmax);

tStage=tic;
%Stage counter:
k=0;
%0-augmented system state:
xau=[x;p;w];
%Reshape hu as a column vector:
hu=reshape(hu,[],1); 
%System dynamics of 0-augmented system:
fxw=[fxw;zeros(np+nw,1)]+[zeros(ns+np,1);w_der(:,1)]; 
%Known inputs contribution of 0-augmented system dynamics:
if numel(fu)~=0
    fu=[fu;zeros(np+nw,numel(fu(1,:)))]; 
end
%Initialize matrix which (i,k)-entry is =0 if the ith state is obs. at the kth stage and =1 if it is not:
unobs_states_ind=ones(max_ns,opts.affine_kmax);
%First codistribution vector:
Delta_Omega=[hxw;hu];

tDer=tic;
%First derivative of codistribution vector (using parallel toolbox):
if opts.affine_parallel_Lie==1
    parfor i=1:nxau(1)
        dif_Delta_omega(:,i)=jacobian(Delta_Omega,xau(i));
    end
    dif_time(1)=toc(tDer);
%First derivative of codistribution vector (without parallel toolbox):
else
    dif_Delta_omega=jacobian(Delta_Omega,xau);
    dif_time(1)=toc(tDer);
end

dif_Omega=dif_Delta_omega;

%Index actualization:
k=k+1;
%0-stage computation time:
stage_time(1)=toc(tStage);
new_stage_time(1)=0;

%(Optional) Define numerical equivalents of the symbolic variables:
if opts.numeric == 1
    [~,~,nz_w_der]=find(w_der);
    xaug=[xau;[f;zeros(np,1);reshape(nz_w_der,[],1)];h;u];
    allvariables = symvar(xaug);
    numbers = vpa(0.1+rand(size(allvariables)));
end  

%==============================Kth stage==================================%
while k<opts.affine_kmax+1 && new_stage_time(k)<opts.affine_tStage

    fprintf('\n >>> Building observability matrix of %d-augmented system...',k)
    
    tStage=tic;
    
    %Actualization of delta_omega (Lie derivative of Omega): 
    Delta_Omega=dif_Delta_omega*[fxw fu]; 
    %Reshape delta_omega as a column vector:
    Delta_Omega=reshape(Delta_Omega,[],1);
    %Including as states only nonzero kth derivatives of unknown inputs:
    [nz_r,~,nz_w_der] = find(w_der(:,k)); 
    xau=[xau;nz_w_der];
    %Number of states of k-augmented system:
    nxau(k+1)=numel(xau);
    %Calculate partial derivatives of Delta_Omega (using parallel toolbox):
    tDer=tic;
    if opts.affine_parallel_Lie==1
        dif_Delta_omega=sym(zeros(numel(Delta_Omega),nxau(k+1)));
        parfor i=1:nxau(k+1)
            dif_Delta_omega(:,i)=jacobian(Delta_Omega,xau(i));
        end
    %Calculate partial derivatives of Delta_Omega (without parallel toolbox):
    else
        dif_Delta_omega=jacobian(Delta_Omega,xau);
    end
    %Actualization of differential of Omega:
    dif_Omega=[dif_Omega zeros(numel(dif_Omega(:,1)),nxau(k+1)-nxau(k)); dif_Delta_omega];
    %Computation time of differential of Omega:
    dif_time(k+1)=dif_time(k)+toc(tDer);
    
%==============Investigate observability of k-augmented system============%

    fprintf('\n >>> Calculating rank of %d-augmented system observability matrix...',k)
    
    %(Optional) Replace known initial conditions:
    if opts.replaceICs == 1
        xind = find(known_ics);
        if size(ics) ~= size(x)
            ics = transpose(ics);
        end
        dif_Omega= subs(dif_Omega,x(xind),ics(xind));
    end 
    
    %(Optional) Replace numerical equivalents of symbolic variables:
    if opts.numeric == 1
        num_dif_omega = subs(dif_Omega,allvariables,numbers); 
    else
        num_dif_omega = dif_Omega;
    end 
    
    %Calculate rank of k-observability matrix:
    tRank=tic;
    rank_dif_omega(k)=rank(num_dif_omega);
    %Actualization of rank computation time:
    tRank=toc(tRank);
    rank_time(k)=rank_time(k)+tRank;
    rank_time(k+1)=rank_time(k);
    
    fprintf('\n     Rank = %d (calculated in %d seconds)',rank_dif_omega(k),tRank)
    
    %Actualization of partial ranks matrix:
    partial_ranks(1:nxau(k+1),k)=rank_dif_omega(k)-1;
    %Almacenate indexes of unobservable states:
    [nf_unobs,~]=find(unobs_states_ind(1:nxau(k+1),k));
    
    %Comprobate if the model is FISPO:
    if rank_dif_omega(k)==nxau(k+1)
        %Print results if the model is FISPO:
        new_stage_time(k+1)=toc(tStage);
        stage_time(k+1)=stage_time(k)+new_stage_time(k+1);
        
        fprintf('\n\n ------------------------ \n');
        fprintf(' >>> RESULTS SUMMARY:\n');
        fprintf(' ------------------------ \n');
        fprintf('\n >>> The model is k-row observable for k = %d  \n',k)
        fprintf('\n >>> The model is Fully Input-State-Parameter Observable (FISPO):');
        if nw>0, fprintf('\n All its unknown inputs are observable.'); end
        fprintf('\n All its states are locally structurally observable.');
        if np>0, fprintf('\n All its parameters are locally structurally identifiable.'); end
        
        totaltime = toc(tStart);
        fprintf('\n Total execution time: %d \n\n',totaltime);
        if opts.affine_graphics == 1
        %All states are k-row observable:
        unobs_states_ind(:,k)=zeros(max_ns,1);
        %Ticks including system states for labelling y-axis:
        xau_ticks=flip(arrayfun(@char, xau, 'uniform',0));                                     
        for i=1:k
            %Unobservable states at k=i:
            unobs_i=find(unobs_states_ind(1:nxau(i+1),i));
            %Number of unobservable states at k=i:
            n_unobs=numel(unobs_i);
            %Plot unobservable and observable states for {1,...,k}:            
            figure(1)
            %Print x-tick at unobservable states:
            scatter(i*ones(n_unobs,1),(nxau(k+1)+1)*ones(n_unobs,1)-unobs_i,50,'x','b','LineWidth',1)
            grid on
            hold on
            %Observable states at k=i:
            obs_i=find(unobs_states_ind(1:nxau(i+1),i)==0);  
            %Number of observable states at k=i:
            n_obs=numel(obs_i);
            %Print o-tick at observable states:
            scatter(i*ones(n_obs,1),(nxau(k+1)+1)*ones(n_obs,1)-obs_i,50,'o','b','LineWidth',1)
            %Undefined states at k=i:
            non_def_states=(nxau(i+1)+1):nxau(k+1);
            %Print black point at undefined states:
            scatter(i*ones(1,nxau(k+1)-nxau(i+1)),(nxau(k+1)+1)*ones(1,nxau(k+1)-nxau(i+1))-non_def_states,70,'.','k','LineWidth',1.5)
        end
        %Plot settings:
        lgd=legend('Unobservable states','Observable states','Non-defined states','location','northeast');
        lgd.FontSize=7;
        axis([0 k+4 0 nxau(k+1)+1]);                        
        xticks(1:k);
        xlabel('Stage');
        yticks(1:nxau(k+1));
        yticklabels(xau_ticks);
        title('Classification of states')
        hold off
        dif_time=dif_time(1:k+1);
        stage_time=stage_time(1:k+1);
        rank_time=rank_time(1:k);
        partial_rank_time=partial_rank_time(1:k);
        %Plot Omega_k partial derivatives computation time for each iteration:
        figure(2)
        subplot(1,4,1)
        semilogy(0:k,dif_time,'LineStyle',':','Linewidth',1.5,'Color',[0 0.75 0.5],'Marker','*','MarkerSize',6)
        xlabel('k')
        ylabel('t[s]')
        xticks(0:k)
        yticks(sort(dif_time))
        title('Partial derivatives computation time.')
   
        %Plot Omega_k rank computation time for each iteration:
        subplot(1,4,2)
        semilogy(1:k,rank_time,'LineStyle',':','Linewidth',1.5,'Color',[0 0.75 0.5],'Marker','*','MarkerSize',6)
        xlabel('k')
        ylabel('t[s]')
        xticks(0:k)
        yticks(sort(rank_time))
        title('Rank computation time')
        
        %Plot partial rank computation time for each iteration:
        subplot(1,4,3)
        semilogy(1:k,partial_rank_time,'LineStyle',':','Linewidth',1.5,'Color',[0 0.75 0.5],'Marker','*','MarkerSize',6)
        xlabel('k')
        ylabel('t[s]')
        xticks(1:k)
        yticks(sort(partial_rank_time(1:k-1)))
        xlim([1 k])
        title('Partial rank computation time')
    
        %Plot total computation time for each iteration:
        subplot(1,4,4)
        semilogy(0:k,stage_time,'LineStyle',':','Linewidth',1.5,'Color',[0 0.75 0.5],'Marker','*','MarkerSize',6)
        xlabel('k')
        ylabel('t[s]')
        xticks(0:k)
        yticks(sort(stage_time))
        title('Stage computation time')
        end
        %Save results if model is FISPO:
        resultsname = sprintf('id_results_%s',modelname);
        fullresultsname = strcat(pwd,filesep,'results',filesep,resultsname,'_',date);
        save(fullresultsname);
        return 
        
    elseif nw==0 && k>1 && rank_dif_omega(k)==rank_dif_omega(k-1)
        
        %If rank(Omega_k)=rank(Omega_k-1), the rank will not increase anymore and system is not observable:
        new_stage_time(k+1)=toc(tStage);
        stage_time(k+1)=stage_time(k)+new_stage_time(k+1);
        fprintf('\n >>> Observability matrix dimension will not increase including further derivatives.')
        fprintf('\n >>> The control system is not k-row observable for k higher or equal to 1.')
        totaltime = toc(tStart);
        
        if opts.affine_graphics == 1
        %Ticks including system states for labelling y-axis:
        xau_ticks=flip(arrayfun(@char, xau, 'uniform',0)); 
        
        for i=1:k
            %Unobservable states at k=i:
            unobs_i=find(unobs_states_ind(1:nxau(i+1),i));
            %Number of unobservable states at k=i:
            n_unobs=numel(unobs_i);
            %Plot unobservable and observable states for {1,...,k}:            
            figure(1)
            %Print x-tick at unobservable states:
            scatter(i*ones(n_unobs,1),(nxau(k+1)+1)*ones(n_unobs,1)-unobs_i,50,'x','b','LineWidth',1)
            %Plot settings:
            axis([0 k+4 0 nxau(k+1)+1]);                        
            xticks(1:k);
            xlabel('k');
            yticks(1:nxau(k+1));
            yticklabels(xau_ticks);
            grid on
            hold on
            %Observable states at k=i:
            obs_i=find(unobs_states_ind(1:nxau(i+1),i)==0);  
            %Number of observable states at k=i:
            n_obs=numel(obs_i);
            %Print o-tick at observable states:
            scatter(i*ones(n_obs,1),(nxau(k+1)+1)*ones(n_obs,1)-obs_i,50,'o','b','LineWidth',1)
            %Undefined states at k=i:
            non_def_states=(nxau(i+1)+1):nxau(k+1);
            %Print black point at undefined states:
            scatter(i*ones(1,nxau(k+1)-nxau(i+1)),(nxau(k+1)+1)*ones(1,nxau(k+1)-nxau(i+1))-non_def_states,70,'.','k','LineWidth',1.5)
        end
        legend('Unobservable states.','Observable states.','Non-defined states.')
        hold off
        title('Observable and unobservable states VS stage number')
        
        dif_time=dif_time(1:k+1);
        stage_time=stage_time(1:k+1);
        rank_time=rank_time(1:k);
        partial_rank_time=partial_rank_time(1:k);
        
        %Plot Omega_k partial derivatives computation time for each iteration:
        figure(2)
        subplot(1,4,1)
        semilogy(0:k,dif_time(1:k+1),'LineStyle',':','Linewidth',1.5,'Color',[0 0.75 0.5],'Marker','*','MarkerSize',6)
        xlabel('k')
        ylabel('t[s]')
        xticks(0:k)
        yticks(sort(dif_time))
        title('Partial derivatives computation time VS stage number')
   
        %Plot Omega_k rank computation time for each iteration:
        %figure(3)
        subplot(1,4,2)
        semilogy(1:k,rank_time,'LineStyle',':','Linewidth',1.5,'Color',[0 0.75 0.5],'Marker','*','MarkerSize',6)
        xlabel('k')
        ylabel('t[s]')
        xticks(0:k)
        yticks(sort(rank_time))
        title('Rank computation time VS stage number.')
        
        %Plot partial rank computation time for each iteration:
        %figure(4)
        subplot(1,4,3)
        semilogy(1:k,partial_rank_time,'LineStyle',':','Linewidth',1.5,'Color',[0 0.75 0.5],'Marker','*','MarkerSize',6)
        xlabel('k')
        ylabel('t[s]')
        xticks(1:k)
        yticks(sort(partial_rank_time(1:k-1)))
        xlim([1 k])
        title('Partial ranks computation time VS stage number.')
    
        %Plot total computation time for each iteration:
        %figure(5)
        subplot(1,4,4)
        semilogy(0:k,stage_time,'LineStyle',':','Linewidth',1.5,'Color',[0 0.75 0.5],'Marker','*','MarkerSize',6)
        xlabel('k')
        ylabel('t[s]')
        xticks(0:k)
        yticks(sort(stage_time))
        title('Stage computation time VS stage number.')
        end    
        %Print observability results:
        tic;
        [nf_obs_states,~]=find(unobs_states_ind(1:ns,k)==0);
        %Observable states:
        obs_states=xau(nf_obs_states);
        unobs_states_ind(1:ns,k)=1;
        %Observable parameters:
        [nf_obs_par,~]=find(unobs_states_ind(1:(ns+np),k)==0);
        obs_par=xau(nf_obs_par);
        
        fprintf('\n\n ------------------------ \n');
        fprintf(' >>> RESULTS SUMMARY:\n');
        fprintf(' ------------------------ \n');
        if numel(obs_states)>0, fprintf('\n >>> The original observable states are: \n   %s', char(obs_states));end
        if numel(obs_par)>0, fprintf('\n >>> The observable parameters are: \n    %s',char(obs_par)); end
        
        totaltime=totaltime+toc;
        fprintf('\n Total execution time: %d \n\n',totaltime);
        resultsname = sprintf('id_results_%s',modelname);
        fullresultsname = strcat(pwd,filesep,'results',filesep,resultsname,'_',date);
        save(fullresultsname);
        return  
    end
    
    fprintf('\n >>> The %d-augmented system is not FISPO.',k)
     
%=========Investigate partial observability of k-augmented system=========%
 
    fprintf('\n >>> Investigating partial observability of %d -augmented system...',k)
    tPartial=tic;
    %Partial rank computation (using parallel toolbox):
    if opts.affine_parallel_rank
        %For the (k-1)-row unobs.states, remove its column and recalculate rank: 
        parfor i=1:numel(nf_unobs)                  
            elim_dif_omega=num_dif_omega;
            elim_dif_omega(:,nf_unobs(i))=[];     
            new_rank(i)=rank(elim_dif_omega);
        end
        %Actualization of obs-unobs states array:
        for i=1:numel(nf_unobs)
            if new_rank(i)<rank_dif_omega(k)
                unobs_states_ind(nf_unobs(i),k)=0;
                partial_ranks(nf_unobs(i),k)=new_rank(i);
            end
        end
        partial_rank_time(k)=partial_rank_time(k)+toc(tPartial);
        partial_rank_time(k+1)=partial_rank_time(k);
        %Partial rank calculation (without parallel toolbox):    
    else
        for i=1:numel(nf_unobs)
            elim_dif_omega=num_dif_omega; 
            elim_dif_omega(:,nf_unobs(i))=[];      
            partial_ranks(nf_unobs(i),k)=rank(elim_dif_omega);
            if partial_ranks(nf_unobs(i),k)<rank_dif_omega(k), unobs_states_ind(nf_unobs(i),k)=0;end
        end
        partial_rank_time(k)=partial_rank_time(k)+toc(tPartial);
        partial_rank_time(k+1)=partial_rank_time(k);
    end
%==================Actualization of system dynamics=======================%

    %Almacenate unobs-obs. indexes at the kth stage.
    unobs_states_ind(:,k+1)=unobs_states_ind(:,k); 
    %Actualization of state-unknown inputs system dynamics:
    fxw=[fxw;w_der(nz_r,k+1)];
    %Actualization of known inputs contribution to system dynamics:
    fu=[fu;zeros(nxau(k+1)-nxau(k),nu)];
    %Actualization of stage computation time:
    new_stage_time(k+1)=toc(tStage);
    stage_time(k+1)=stage_time(k)+new_stage_time(k+1);
    %Index actualization:
    k=k+1;
end 

%================Observability results if k=opts.kmax+1===================%

totaltime = toc(tStart);

if opts.affine_graphics==1
%Ticks with states for labelling y-axis:
xau_ticks=flip(arrayfun(@char, xau, 'uniform',0));
%Plot unobservable and observable states for k=1,...,opts.kmax:
figure(1)
for i=1:k-1
    %Unobservable states at the ith stage:
    unobs_i=find(unobs_states_ind(1:nxau(i+1),i));     
    %Number of unobservable states at the ith stage:
    n_unobs=numel(unobs_i);
    %Print x-tick at unobservable states:
    scatter(i*ones(n_unobs,1),(nxau(k)+1)*ones(n_unobs,1)-unobs_i,50,'x','b','LineWidth',1)
    grid on
    hold on
    %Observable states at the ith stage:
    obs_i=find(unobs_states_ind(1:nxau(i+1),i)==0);  
    %Number of observable states at the ith stage:
    n_obs=numel(obs_i);
    %Print o-tick at observable states:
    scatter(i*ones(n_obs,1),(nxau(k)+1)*ones(n_obs,1)-obs_i,50,'o','b','LineWidth',1)
    %Undefined states at the ith stage:
    non_def_states=(nxau(i+1)+1):nxau(k);
    %Print black dot at undefined states:
    scatter(i*ones(1,nxau(k)-nxau(i+1)),(nxau(k)+1)*ones(1,nxau(k)-nxau(i+1))-non_def_states,70,'.','k','LineWidth',1.5)
end
legend('Unobservable states.','Observable states.','Non-defined states.')
hold off
%Plot settings:
title('Observable and unobservable states VS stage number')
xlabel('Stage')
axis([0 k+4 0 nxau(k)+1]);                        
xticks(1:k-1);
yticks(1:nxau(k));
yticklabels(xau_ticks);

%Plot Omega_k partial derivatives computation time:
figure(2)
subplot(1,4,1)
semilogy(0:k-1,dif_time,'LineStyle',':','Linewidth',1.5,'Color',[0 0.75 0.5],'Marker','*','MarkerSize',6)
xlabel('Stage')
ylabel('t[s]')
xticks(0:k-1)
yticks(dif_time(1:k))
xlim([0 k-1])
title('Partial derivatives computation time')

%Plot omega_k rank computation time:
%figure(3)
subplot(1,4,2)
semilogy(1:k-1,rank_time(1:k-1),'LineStyle',':','Linewidth',1.5,'Color',[0 0.75 0.5],'Marker','*','MarkerSize',6)
xlabel('Stage')
ylabel('t[s]')
xticks(1:k-1)
xlim([1 k-1])
yticks(rank_time(1:k-1))
title('Rank computation time')

%Plot partial rank computation time:
%figure(4)
subplot(1,4,3)
semilogy(1:k-1,partial_rank_time(1:k-1),'LineStyle',':','Linewidth',1.5,'Color',[0 0.75 0.5],'Marker','*','MarkerSize',6)
xlabel('Stage')
ylabel('t[s]')
xticks(1:k-1)
yticks(partial_rank_time(1:k-1))
xlim([1 k-1])
title('Partial ranks computation time')

%Plot total execution time:
subplot(1,4,4)
semilogy(0:k-1,stage_time(1:k),'LineStyle',':','Linewidth',1.5,'Color',[0 0.75 0.5],'Marker','*','MarkerSize',6)
xlabel('Stage')
xticks(0:k-1)
xlim([0 k-1])
ylabel('t[s]')
yticks(stage_time(1:k))
title('Stage computation time')
end

tic;
%Print results at k=opts.kmax:
[nf_obs_states,~]=find(unobs_states_ind(1:ns,k)==0);
obs_states=xau(nf_obs_states);

unobs_states_ind(1:ns,k)=1;
[nf_obs_par,~]=find(unobs_states_ind(1:(ns+np),k)==0);
obs_par=xau(nf_obs_par);

unobs_states_ind((ns+1):(ns+np),k)=1;
[nf_obs_inputs,~]=find(unobs_states_ind(:,k)==0);
obs_inputs=xau(nf_obs_inputs);
if k>opts.affine_kmax
    fprintf('\n >>> Maximum number of iterations or rebased.');
    fprintf('\n >>> The model is not k-row observable for k < %d  \n',k);
else
    fprintf('\n >>> Maximum computation time rebased.');
end
fprintf('\n\n ------------------------ \n');
fprintf(' >>> RESULTS SUMMARY:\n');
fprintf(' ------------------------ \n');
if numel(obs_states)>0, fprintf('\n >>> The original observable states are: \n   %s', char(obs_states));end
if numel(obs_par)>0, fprintf('\n >>> The observable parameters are: \n    %s',char(obs_par)); end
if numel(obs_inputs)>0, fprintf('\n >>> The observable unknown inputs are: \n   %s',char(obs_inputs));end

totaltime=totaltime+toc;
fprintf('\n Total execution time: %d \n\n',totaltime);
resultsname = sprintf('id_results_%s',modelname);

fullresultsname = strcat(pwd,filesep,'results',filesep,resultsname,'_',date);
save(fullresultsname);

end
