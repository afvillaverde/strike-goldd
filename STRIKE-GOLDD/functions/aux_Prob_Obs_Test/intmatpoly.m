%--------------------------------------------------------------------------
% function that computes the integration of polynomials as vector
% in the colums of a matrix
%--------------------------------------------------------------------------

function sol=intmatpoly(Mat,Order)

[~,n]=size(Mat);
sol=zeros(Order,n);

% trun polynomns
for i=1:n
    sol(:,i)=trunpoly(polyint(Mat(:,i)')',Order);
end