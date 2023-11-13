function ctrl_LARC_GSC(modelname,opts,n,nu,f1,f2,jac_f1,jac_f2,x,x0)
% ======================================================================= %
% =========================Function that checks========================== %
% ===================the Lie Algebraic Rank Condition==================== % 
% ============= and Sussmann's General Sufficient Condition ============= %
% ======================================================================= %

tLARC = tic;

fprintf(['\n Computing the LIE ALGEBRAIC RANK CONDITION ' ...
    '(LARC)...\n Computing Lie bracket number']);

% ==============================PREPROCESS=============================== %
Di = [f1,f2];                       % Initial distribution
Di_x0 = subs(Di,x,x0);              % Distribution at the point
Di_rank = rank(Di);                 % Rank of the initial distribution
Di_x0_rank = rank(Di_x0);           % Rank of the distribution at the point
Di_cols = size(Di,2);               % Number of colums of the distribution
Lie_index = [];                     % Index for GSC
last_brackets_index = eye(Di_cols); % Index for GSC of last brackets
last_brackets_index(:,1) = []; 
last_brackets = f2;               % First vectors to construct the brackets
last_brackets_numb = size(f2,2);  % Number of vectors to construct brackets
                         % ============================================== %
                         % ====Check if the following test properly ===== %
                         % ======determines if the point is regular====== %
                         % ============================================== %
if Di_rank == Di_x0_rank % Regularity of x0 for the initial distribution
    regular = 1;
else 
    regular = 0;
end

% =======================DISTRIBUTION COMPUTATION======================== %
% Initialization
increaseLie=1;
iLie = 1;
lasttime = 0;
fprintf(' %d... ',iLie)
% Loop
while increaseLie == 1 && regular == 1 && lasttime < opts.maxtime % Iterations while not involutive 
                                       % and regular
    % ======================NEW BRACKETS COMPUTATION===================== %                                
    new_brackets = [Lie_bracket(f1,last_brackets,x,n,jac_f1), ... % Construction 
                    Lie_bracket(f2,last_brackets,x,n,jac_f2)];    % of brackets
    % Indexes for the new brackets
    new_brackets_index=zeros(nu+1,last_brackets_numb*(nu+1)); 
    for i=0:nu
        aux_brackets_index=last_brackets_index;
        aux_brackets_index(i+1,:)=aux_brackets_index(i+1,:)+eye(1,last_brackets_numb);
        new_brackets_index(:,i*last_brackets_numb+1:(i+1)*last_brackets_numb) = aux_brackets_index;
    end
    % Number of new brackets
    new_brackets_numb = last_brackets_numb*(nu+1);
    % =============Remove null columns from the new brackets============= %
    zero_cols_ind = [];
    for icol=1:size(new_brackets,2) 
        if new_brackets(:,icol) == zeros(n,1) % Check for vectors of zeros
            zero_cols_ind = [zero_cols_ind,icol]; %#ok<AGROW> 
        end
    end
    % Updates
    new_brackets(:,zero_cols_ind) = [];
    new_brackets_index(:,zero_cols_ind) = [];
    new_brackets_numb = new_brackets_numb - length(zero_cols_ind);
    % =====================COMPUTE NEW DISTRIBUTION====================== %
    Dc = [Di,new_brackets];
    Dc_x0 = [Di_x0,subs(new_brackets,x,x0)];
    Dc_rank = rank(Dc);
    Dc_x0_rank = rank(Dc_x0);
    Dc_cols = Di_cols + new_brackets_numb;
    % ==============VERIFICATION OF DIFFERENT REQUIREMENTS=============== %
    if new_brackets_numb == 0 % If no new non zero brackets STOP
        increaseLie = 0;
    elseif Dc_rank ~= Dc_x0_rank % If non regular point change loop 
        regular = 0;
    else
        if Dc_rank == n % If full rank reached STOP 
            increaseLie = 0;
        elseif Dc_rank == Di_rank % If rank not increased STOP
            increaseLie = 0;
            % updates
            Dc = Di; 
            Dc_x0 = Di_x0; 
            Dc_cols = Di_cols;
            new_brackets_numb = 0;
            new_brackets_index = [];
            new_brackets = [];
        end
        if Dc_cols ~= Dc_rank % If colums dependent remove them
            colum_depend_index = [];
            Dc_check_old = Dc;
            for icol=Di_cols+1:Dc_cols
                Dc_check = Dc;
                Dc_check(:,[colum_depend_index,icol]) = [];
                if rank(Dc_check_old) == rank(Dc_check)
                    colum_depend_index = [colum_depend_index,icol]; %#ok<AGROW> 
                    Dc_check_old = Dc;
                    Dc_check_old(:,colum_depend_index) = [];
                end
            end
            % Updates
            Dc = Dc_check_old;
            Dc_x0(:,colum_depend_index) = [];
            new_brackets(:,colum_depend_index-Di_cols) = [];
            new_brackets_index(:,colum_depend_index-Di_cols) = [];
            new_brackets_numb = new_brackets_numb-length(colum_depend_index);
        end
    end   
    % ==============================UPDATE=============================== %
    Di = Dc;
    Di_x0 = Dc_x0;
    Di_rank = Dc_rank;
    Di_cols = Di_cols+new_brackets_numb;
    
    last_brackets = new_brackets;
    Lie_index = [Lie_index,new_brackets_index]; %#ok<AGROW> 
    last_brackets_index = new_brackets_index;
    last_brackets_numb = new_brackets_numb;
    iLie = iLie + 1;
    fprintf(' %d... ',iLie)  
    lasttime = toc(tLARC);
end

base = Di;
base_index = (1:size(base,2));

while increaseLie == 1 && lasttime < opts.maxtime %construct distribution i when x0 not regular
    % ======================NEW BRACKETS COMPUTATION===================== %                                
    new_brackets = [Lie_bracket(f1,last_brackets,x,n,jac_f1),Lie_bracket(f2,last_brackets,x,n,jac_f2)];
    % Indexes for the new brackets
    new_brackets_index=zeros(nu+1,last_brackets_numb*(nu+1));
    for i=0:nu
        aux_brackets_index=last_brackets_index;
        aux_brackets_index(i+1,:)=aux_brackets_index(i+1,:)+ones(1,last_brackets_numb);
        new_brackets_index(:,i*last_brackets_numb+1:(i+1)*last_brackets_numb) = aux_brackets_index;
    end
    new_brackets_numb = last_brackets_numb*(nu+1);
    % =============Remove null columns from the new brackets============= %
    zero_cols_ind = [];
    for icol=1:size(new_brackets,2) 
        if new_brackets(:,icol) == zeros(n,1) % Check for vectors of zeros
            zero_cols_ind = [zero_cols_ind,icol]; %#ok<AGROW> 
        end
    end
    new_brackets(:,zero_cols_ind) = [];
    new_brackets_index(:,zero_cols_ind) = [];
    % Number of new brackets
    new_brackets_numb = new_brackets_numb - length(zero_cols_ind);
    % ==============Remove some of the dependent brackets================ %
    % Initialization
    last_brackets = simplify(new_brackets);
    new_brackets = [];
    last_brackets_numb = new_brackets_numb;
    new_brackets_numb = 0;
    last_brackets_index = new_brackets_index;
    new_brackets_index = [];
    Dc = Di;
    % Find dependencies
    for i=1:last_brackets_numb % Check each new bracket
        % Initializations
        logic=0;
        base_cols = size(base,2);
        base_index_dell = [];
        for j=1:base_cols % Increase gradually the base colums
            aux = [base(:,end-j+1:end),last_brackets(:,i)];
            check = null(aux); % Compute matrix of possible dependance
            if ~isempty(check) % Check if there exists any dependance
                check_end = check(:,end);
                try % Check for some smooth dependance
                    subs(check_end,x,x0);
                    logic=1; % Smooth dependance 
                    break
                catch % Check for some smooth dependance to reduce base
                    index = find(check_end);
                    if length(index)==2
                        auxx = aux;
                        auxx(:,index(1)) = aux(:,index(2));
                        auxx(:,index(2)) = aux(:,index(1));
                        check_aux = null(auxx);
                        try % Reduction of the base
                            subs(check_aux(:,end),x,x0);
                            base_index_dell = [base_index_dell,base_cols-j+index(1)]; %#ok<AGROW> 
                        catch

                        end
                    end
                end
            end
        end
            if logic == 0 % Dependance not found
                % Updates
                Dc = [Dc,last_brackets(:,i)]; %#ok<AGROW> Update distribution 
                base(:,base_index_dell) = [];
                base_index(base_index_dell) = [];
                base_index = [base_index,size(Dc,2)]; %#ok<AGROW> 
                base = [base,last_brackets(:,i)]; %#ok<AGROW> Update base
                new_brackets_numb=new_brackets_numb+1;
                new_brackets = [new_brackets,last_brackets(:,i)]; %#ok<AGROW> 
                new_brackets_index = [new_brackets_index,last_brackets_index(:,i)]; %#ok<AGROW> 
            end
    end
    % =====================COMPUTE NEW DISTRIBUTION====================== %
    Dc = [Di,new_brackets];
    Dc_x0 = [Di_x0,subs(new_brackets,x,x0)];
    Dc_x0_rank = rank(Dc_x0);
    Dc_cols = Di_cols + new_brackets_numb;
    % ==============VERIFICATION OF DIFFERENT REQUIREMENTS=============== %
    if new_brackets_numb == 0 || Dc_x0_rank == n % If no new non zero brackets
                                                 % or full rank reached STOP 
        increaseLie = 0;
    end
    % ==============================UPDATE=============================== %
    Di = Dc;
    Di_x0 = Dc_x0;
    Di_cols = Dc_cols;
    
    last_brackets = new_brackets;
    Lie_index = [Lie_index,new_brackets_index]; %#ok<AGROW> 
    last_brackets_index = new_brackets_index;
    last_brackets_numb = new_brackets_numb;
    iLie = iLie + 1;
    fprintf(' %d... ',iLie)
    lasttime = toc(tLARC);
end

% =============================LARC display============================== %
if opts.LARC == 1
    fprintf('\n dim(Delta_c(x)) = %d',Dc_x0_rank);
    if Dc_x0_rank == n % Reachable model
        fprintf(['   => The LARC is fulfilled => ' ...
            'Model %s is accessible.\n'], modelname);    
    else % Unreachable model
        if lasttime >= opts.maxtime
            fprintf('\n    => More Lie brackets would be needed to see if the model is accessible.');
            fprintf('\n    However, the maximum computation time allowed is reached.');
            fprintf('\n    You can increase it by changing <<opts.maxtime>> (currently opts.maxtime = %d)\n',opts.maxtime);
            fprintf('\n    GSC can not be computed since LARC stoped before getting a conclusion');
            opts.GSC = 0;
        else
            fprintf(['   => The LARC is not fulfilled => ' ...
                'Model %s is inaccessible.\n'], modelname);  
        end
    end
    fprintf(' (LARC computation completed in %d seconds.)\n ',toc(tLARC));
end

% ==============================GSC display============================== %
if opts.GSC == 1 
    if Dc_x0_rank == n
        tGSC = tic;
        Dc_x0_brackets = Dc_x0(:,nu+2:end);
        if isempty(Lie_index) % GSC fulfilled
            fprintf('\n The GSC is fulfilled => model %s is STLC.', modelname);
        else % Bad brackets exist
            % Finding bad brackets
            bad_brackets_ind_aux = sum(abs(Dc_x0_brackets));
            bad_brackets_ind_aux = find(bad_brackets_ind_aux);
            bad_brackets_ind = Lie_index(:,bad_brackets_ind_aux);
            bad_brackets_ind(1,:) = mod(bad_brackets_ind(1,:),2);
            bad_brackets_ind(2:end,:) = mod(1.+bad_brackets_ind(2:end,:),2);
            bad_brackets_ind = sum(bad_brackets_ind);
            bad_brackets_ind = nu+1.-bad_brackets_ind;
            bad_brackets_ind = bad_brackets_ind_aux(bad_brackets_ind == 0);
            if isempty(bad_brackets_ind) % GSC fulfilled
                fprintf('\n The GSC is fulfilled => model %s is STLC.', ...
                    modelname);
            else % Bad brackets found
                logic2 = 0;
                for i=1:length(bad_brackets_ind)
                    check = null([base(:,base_index < 1+nu+bad_brackets_ind(i)),Dc(:,1+nu+bad_brackets_ind(i))]);
                    if size(check,2) < 2 % Proportional comb. already checked
                        logic = 0;
                    else % Check for possible combinations
                        check_end = check(:,end);
                        index = find(check_end);
                        for j=1:length(index) % Check for combinations of null vectors elements
                            try
                            subs(check_end(index(j)),x,x0); % Comb. not needed
                            catch % Check for possible combinations
                                [num,dem] = frac_elem_sym(check_end(index(j)));
                                logic = 0;
                                for k=1:size(check,2)-1 % Check for pos. colums
                                    [num2,dem2] = frac_elem_sym(check(j,end-k));
                                    try
                                        aux = num*dem2/(dem*num2);
                                        value = subs(aux,x,x0);
                                        if isinf(value) || isnan(value)
                                            break 
                                        end
                                        check_end=check_end-aux*check(:,end-k);
                                        logic = 1; % Element cobination found
                                    catch
                                    end
                                    if logic == 1 % Stop looking
                                        break
                                    end
                                end
                            end
                        end
                    end
                    if logic == 0 % Bad bracket combination not found
                        fprintf('\n Bracket [')
                        for j=1:nu+1
                            fprintf('%d',Lie_index(j,bad_brackets_ind(i)))
                            if j~=nu+1
                                fprintf(',')
                            end
                        end
                        fprintf('] does not fulfill the GSC.')
                        logic2 = 1;
                    end
                end
                if logic2 == 0 % Combinations for all bad brackets found
                    fprintf('\n Found combinations for the bad brackets.')
                    fprintf('\n The GSC is fulfilled => model %s is STLC.', modelname);
                end
            end
        end
        fprintf(' (GSC check completed in %d seconds.)\n ',toc(tGSC));
        fprintf(' (GSC total computation completed in %d seconds.)\n ',toc(tLARC));
    else
        fprintf('Model %s is not STLC since it is not accesible', modelname);
    end
end
end