% function that calculates the derivatives of polynomials as vectors
% on columns of a matrix
function sol=dermatpoly(Mat,Order)

[~,n]=size(Mat);
sol=zeros(Order,n);

for i=1:n
    sol(:,i)=trunpoly(polyder(Mat(:,i)')',Order);
end