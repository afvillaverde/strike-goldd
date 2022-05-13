% function that replicates the behavior of the function mod in maple
function sol=ratmod(x,p)
[row,col]=size(x);
if (row>1) || (col>1)
    sol=zeros(row,col);
    for i=1:row
        for j=1:col
            sol(i,j)=ratmod(x(i,j),p);
        end
    end
else
    xr=round(x);
        if xr==x
            sol=mod(x,p);
        else
            a=1;
            b=1;
            tol=1/p;
            ex=10;
            while a/b~=x && tol~=0
                [a,b]=rat(x,tol);
                tol=10^-ex;
                ex=ex*10;
            end
            [~,sol]=gcd(b,p);
            sol=mod(sol*a,p);
        end
end