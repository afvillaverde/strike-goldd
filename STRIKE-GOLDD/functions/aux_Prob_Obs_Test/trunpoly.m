%--------------------------------------------------------------------------
% restrict the order of polynomials as vector on the columns of a matrix
%--------------------------------------------------------------------------

function T=trunpoly(Mat,Order)

[r,c]=size(Mat);
if r==Order
    T=Mat;
elseif r<Order
    T=[zeros(Order-r,c);Mat];
else
    T=Mat(r-Order+1:r,:);
end