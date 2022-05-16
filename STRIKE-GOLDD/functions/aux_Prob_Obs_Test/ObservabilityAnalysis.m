%--------------------------------------------------------------------------
% Function that runs the analysis of the observability-identifiability
% matrix that was previously calculated.
%--------------------------------------------------------------------------

function [observable,unobservable] = ObservabilityAnalysis(X,Theta,ObservabilityMatrix,p,n,l)

observable   = [];
unobservable = [];
xtw          = [Theta;X];

r=gfrank2(ObservabilityMatrix,p);
fprintf('The rank of the observability matrix is %d \n', r);
trnsdeg=n+l-r;
fprintf('transcendence degree of k(U,Y) --> k(U,Y,x,p) is %d \n', trnsdeg);

if trnsdeg==0
    fprintf('\n >>> The model is structurally identifiable and observable.');
    observable = xtw;
else
    fprintf('\n >>> The model is not fully structurally identifiable/observable.\n');
    for i=1:n+l
            O=ObservabilityMatrix;
            O(:,i)=[];
        if gfrank2(O,p)==r
            unobservable=[unobservable,xtw(i)];
        else
            observable=[observable,xtw(i)];
        end
    end
    fprintf('The identifiable/observable variables are:');
    display(arrayfun(@char,observable,'UniformOutput',0));
    fprintf('The unidentifiable/unobservable variables are:');
    display(arrayfun(@char,unobservable,'UniformOutput',0));
end
    