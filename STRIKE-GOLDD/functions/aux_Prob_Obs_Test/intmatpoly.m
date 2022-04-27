% function that computes the intgration of polynomials as vector
% in the colums of a matrix
function sol=intmatpoly(Mat,Order)

[~,n]=size(Mat);
sol=zeros(Order,n);

% truncar polinomios
for i=1:n
    sol(:,i)=trunpoly(polyint(Mat(:,i)')',Order);
end