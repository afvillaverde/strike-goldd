%--------------------------------------------------------------------------
% application of the mod function (as in maple) to the coefficients of a
% symbolic ecuation
%--------------------------------------------------------------------------

function sol=polmod(x,p,logic)

[row,col]=size(x);
if (row>1) || (col>1)
    sol=sym('sol',[row,col]);
    for i=1:row
        for j=1:col
            sol(i,j)=polmod(x(i,j),p,logic);
        end
    end
elseif isSymType(x,'number')
    xr=round(x);
        if xr==x
            sol=mod(x,p);
        else
            if logic == 1
                a=1;
                b=1;
                ex=20;
                tol=10^-ex;
                while a/b~=c(i) && ex<=10^8
                    [a,b]=rat(c(i),tol);
                    tol=10^-ex;
                    ex=ex^2;
                end
                [~,sol]=gcd(b,p);
                sol=mod(sol*a,p);
            else
                [a,b]=rat(x,1e-20);
                [~,sol]=gcd(b,p);
                sol=mod(sol*a,p);
            end
        end
else
    [c,t]=coeffs(x);
    cl=length(c);
    cmod=zeros(1,cl);
    for i=1:cl
        cr=round(c(i));
        if cr==c(i)
            cmod(i)=mod(c(i),p);
        else
            if logic == 1
                a=1;
                b=1;
                ex=20;
                tol=10^-ex;
                while a/b~=c(i) && ex<=10^8
                    [a,b]=rat(c(i),tol);
                    tol=10^-ex;
                    ex=ex^2;
                end
                [~,sol]=gcd(b,p);
                cmod(i)=mod(sol*a,p);
            else
                [a,b]=rat(c(i),1e-20);
                [~,sol]=gcd(b,p);
                cmod(i)=mod(sol*a,p);
            end
        end
    end
    sol=sum(cmod.*t);
end