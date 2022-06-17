%--------------------------------------------------------------------------
% funtion to multiply polynomials as vector on the columns of a matrix
%--------------------------------------------------------------------------

function sol=multmatpolytrun(mat1,mat2,order,nr1,nc1,nr2,nc2)

if nargin==4
    s=nr1;
    ind=1;
    
    sol=zeros(order,nr1*nr1);
    for i=1:nr1
        for j=1:nr1
            aux=zeros(order,1);
            for k=0:nr1-1
                aux=aux+trunpoly(conv(mat1(:,1+s*(i-1)+k), ...
                    mat2(:,j+s*k)),order);
            end
            sol(:,ind)=aux;
            ind=ind+1;
        end
    end
else
    if nc1~=nr2
        error('error en multmatpolytrun matriz dimension must agry')
    else
        ind=1;
    
        sol=zeros(order,nr1*nc2);
        for i=1:nr1
            for j=1:nc2
                aux=zeros(order,1);
                for k=0:nc1-1
                    aux=aux+trunpoly(conv(mat1(:,1+nc1*(i-1)+k), ...
                        mat2(:,j+nc2*k)),order);
                end
                sol(:,ind)=aux;
                ind=ind+1;
            end
        end
    end
end
