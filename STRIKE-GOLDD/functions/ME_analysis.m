%=========================================================================%
%========ME_analysis: Convert a model into Multi-Experiment form==========%

function ME_analysis(modelname,opts)

%============================Load model===================================%
load(modelname); %#ok<*LOAD>

%======================Dimensions of the problem==========================%
% Number of states:
n=numel(x); %#ok<*NODEF>
% Number of outputs:
m=numel(h);
% Number of known inputs:
if exist('u','var')
    nu=numel(u);
else
    u=[];
    nu=0;
end
% Number of unknown inputs:
if exist('w','var')
    nw=numel(w);
else
    w=[];
    nw=0;
end
% Number of initial conditions:
if exist('ics','var')
    nics=numel(ics);
else
    ics=[];
    nics=0;
end
%========================Initialize variables=============================%

me_x   = sym(zeros(n*opts.numexp,1));
me_h   = sym(zeros(m*opts.numexp,1));
me_f   = sym(zeros(n*opts.numexp,1));
me_u   = sym(zeros(nu*opts.numexp,1));
me_w   = sym(zeros(nw*opts.numexp,1));
if exist('ics','var')
    me_ics = zeros(nics*opts.numexp,1);
    if exist('known_ics','var')
        me_known_ics=zeros(nics*opts.numexp,1);
    else
        me_known_ics=[];
    end
else
    me_ics=[];
end

% Replicate initial conditions:
for ind_ics=1:nics
    me_ics=repmat(ics,1,opts.numexp);
    me_known_ics=repmat(known_ics,1,opts.numexp);
end


%========================Multi-experiment model===========================%

variables=[x;w;u];

for i=1:opts.numexp
    
    num_exp=strcat('Exp',num2str(i));
    
    for ind_u=1:nu,me_u(ind_u+nu*(i-1)) = sym(strcat(char(u(ind_u)),sprintf(num_exp)));end
    for ind_w=1:nw,me_w(ind_w+nw*(i-1)) = sym(strcat(char(w(ind_w)),sprintf(num_exp)));end
    for ind_x=1:n,me_x(ind_x+n*(i-1))= sym(strcat(char(x(ind_x)),sprintf(num_exp)));end

    aug_variables=[me_x(1+n*(i-1):n*i);me_w(1+nw*(i-1):nw*i);me_u(1+nu*(i-1):nu*i)];
    
    for ind_f=1:n,me_f(ind_f+n*(i-1))= subs(f(ind_f),variables,aug_variables);end
    for ind_h=1:m,me_h(ind_h+m*(i-1))= subs(h(ind_h),variables,aug_variables);end
end

% New variables: 
u=me_u;
w=me_w;
x=me_x;
ics=me_ics;
known_ics=me_known_ics;
% Multi-experiment dynamics:
f=me_f;
% Multi-experiment output:
h=me_h;

if exist('p','var')==0
    p=[];
end

save(strcat(pwd,filesep,'models',filesep,strcat(modelname,'_',num2str(opts.numexp)),'Exp'),'x','p','u','w','ics','known_ics','f','h')  
