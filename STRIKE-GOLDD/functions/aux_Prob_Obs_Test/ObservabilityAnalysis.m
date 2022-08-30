%--------------------------------------------------------------------------
% Function that runs the analysis of the observability-identifiability
% matrix that was previously calculated.
%--------------------------------------------------------------------------

function [p_id, p_un,obs_states,unobs_states,obs_inputs,unobs_inputs, ...
    meas_x,r] = ObservabilityAnalysis(X,Theta,ObservabilityMatrix, ...
    p,n,l,nwl,h,u)

p_id         = [];        
p_un         = [];          
obs_states   = [];       
unobs_states = []; 
obs_inputs   = [];       
unobs_inputs = []; 
xtw          = [Theta;X];

r=gfrank2(ObservabilityMatrix,p);
fprintf(' >>> The rank of the observability matrix is %d \n', r);
trnsdeg=n+l-r;
fprintf(['     transcendence degree of k(U,Y) --> k(U,Y,x,p)' ...
    ' is %d \n'], trnsdeg);

if trnsdeg==0
    fprintf(['\n >>> The model is structurally identifiable, ' ...
        'observable and reconstrutible.']);
    if nwl>0, fprintf('\n     All its unknown inputs are observable.'); end
	fprintf('\n     All its states are observable.');
	fprintf(['\n     All its parameters are locally structurally' ...
        ' identifiable.']);
    p_id=xtw(1:l);
    obs_states=xtw(1+l:l+n-nwl);  
    obs_inputs=xtw(1+l+n-nwl:l+n);
    %=======================================================================
    % Check which states are directly measured, if any: 
    saidas     = cell(numel(h),1);
    ismeasured = zeros(1,numel(obs_states));
    for i=1:numel(h)
        saidas{i} = char(h(i));
    end
    for i=1:numel(obs_states)
        ismeasured(i) = sum(strcmp(char(obs_states(i)),saidas));
    end
    meas_x = obs_states(ismeasured~=0);     % names of the measured states
    obs_states=obs_states(ismeasured==0);
else
    fprintf(['\n >>> The model is not observable, identifiable or' ...
        ' reconstructible:']);
    for i=1:l
            O=ObservabilityMatrix;
            O(:,i)=[];
        if gfrank2(O,p)==r
            p_un=[p_un,xtw(i)]; %#ok<AGROW>
        else
            p_id=[p_id,xtw(i)]; %#ok<AGROW>
        end
    end
    for i=1+l:l+n-nwl
            O=ObservabilityMatrix;
            O(:,i)=[];
        if gfrank2(O,p)==r
            unobs_states=[unobs_states,xtw(i)]; %#ok<AGROW> 
        else
            obs_states=[obs_states,xtw(i)]; %#ok<AGROW> 
        end
    end
    for i=1+l+n-nwl:l+n
            O=ObservabilityMatrix;
            O(:,i)=[];
        if gfrank2(O,p)==r
            unobs_inputs=[unobs_inputs,xtw(i)]; %#ok<AGROW> 
        else
            obs_inputs=[obs_inputs,xtw(i)]; %#ok<AGROW> 
        end
    end

    %======================================================================
    % Check which states are directly measured, if any: 
    saidas     = cell(numel(h),1);
    ismeasured = zeros(1,numel(obs_states));
    for i=1:numel(h)
        saidas{i} = char(h(i));
    end
    for i=1:numel(obs_states)
        ismeasured(i) = sum(strcmp(char(obs_states(i)),saidas));
    end
    meas_x = obs_states(ismeasured~=0);     % names of the measured states
    obs_states=obs_states(ismeasured==0);

    if numel(p_un)>0 
        if numel(p_id)>0          
            fprintf(['\n >>> These parameters are identifiable:\n' ...
                '      %s '],char(p_id)); 
        end          
        fprintf(['\n >>> These parameters are unidentifiable:\n' ...
            '      %s '],char(p_un));
    else
        fprintf('\n >>> The model is structurally identifiable.');
	    fprintf(['\n     All its parameters are locally structurally' ...
            ' identifiable.']);
    end

    if numel(unobs_states)>0  	
        if numel(obs_states)>0    
            fprintf(['\n >>> These states are observable (and their' ...
                ' initial conditions, if considered unknown, are' ...
                ' identifiable):\n      %s '],char(obs_states)); 
        end 
        fprintf(['\n >>> These states are unobservable (and their ' ...
            'initial conditions, if considered unknown, are ' ...
            'unidentifiable):\n      %s '],char(unobs_states));
        if numel(meas_x)>0        
            fprintf(['\n >>> These states are directly measured:\n' ...
                '      %s '],char(meas_x)); 
        end
    else
        fprintf('\n >>> The model is observable.');
	    fprintf('\n     All its states are observable.');
	    
    end
    
    if nwl>0 
        if numel(unobs_inputs)>0  
        if numel(obs_inputs)>0    
                fprintf(['\n >>> These unmeasured inputs are ' ...
                    'observable:\n      %s '],char(obs_inputs)); 
        end
            fprintf(['\n >>> These unmeasured inputs are ' ...
                'unobservable:\n      %s '],char(unobs_inputs));
        else
            fprintf('\n >>> The model is reconstrutible.');
            fprintf('\n     All its unknown inputs are observable.')
        end
    end
    if numel(u)>0             
        fprintf('\n >>> These inputs are known:\n      %s \n',char(u)); 
    end

end
    