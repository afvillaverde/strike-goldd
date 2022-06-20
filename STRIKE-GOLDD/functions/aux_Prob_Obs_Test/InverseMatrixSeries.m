%--------------------------------------------------------------------------
% use of newton operator for computation of InvHomSolInvA
%--------------------------------------------------------------------------

function sol=InverseMatrixSeries(Mat,LocalOrder,MatSize,MyPrime)

if LocalOrder==1
    sol=Mat(end,:);
    aux=zeros(MatSize,MatSize);
    for i=1:MatSize
        aux(i,:)=sol((i-1)*MatSize+1:i*MatSize);
    end
    aux=inv(aux);    
    for i=1:MatSize
        sol((i-1)*MatSize+1:i*MatSize)=aux(i,:);
    end    
else
    NewOrder = 2^ceil(log2(LocalOrder)-1);
    aux=trunpoly(Mat,NewOrder);
    aux=InverseMatrixSeries(aux,NewOrder,MatSize,MyPrime);
    aux=trunpoly(aux,LocalOrder);

    sol=mod(multmatpolytrun(-aux,Mat,LocalOrder,MatSize),MyPrime);
    sol(end,(1:MatSize:MatSize*MatSize)+(0:MatSize-1))= ...
        sol(end,(1:MatSize:MatSize*MatSize)+(0:MatSize-1))+2;
    sol=mod(multmatpolytrun(sol,aux,LocalOrder,MatSize),MyPrime);
end
