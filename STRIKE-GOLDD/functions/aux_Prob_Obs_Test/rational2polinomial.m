%--------------------------------------------------------------------------
% Function that changes a rational equation into a polynomial one
%--------------------------------------------------------------------------

function [a,b]=rational2polinomial(exp)

lexp=length(exp);
if lexp==1
    if isSymType(exp,'number') || isSymType(exp,'variable')
        a=exp;
        b=1;
    elseif isSymType(exp,'power')
        c=children(exp);
        if iscell(c)==1
            c=[c{:}];
        end
        if isnumeric(eval(c(2)))
            if round(c(2))~=c(2)
                c(2)=round(c(2));
                warning('MATLAB:exponetial',['Found noninteger ' ...
                    'exponents in the equations replaced with their' ...
                    ' rounded value'])
                warning('off','MATLAB:exponetial')
            end
            if c(2)<0
                [aux1,aux2]=rational2polinomial(c(1));
                a=aux2^(-c(2));
                b=aux1^(-c(2));
            elseif c(2)==0
                a=1;
                b=1;
            else
                [aux1,aux2]=rational2polinomial(c(1));
                a=aux1^(c(2));
                b=aux2^(c(2));
            end
        else
            error(['The model is not rational: it cannot be analysed' ...
                ' with this algorithm. To use the FISPO algorithm' ...
                ' instead, set in options.m: opts.algorithm = 1'])
        end
    elseif isSymType(exp,'times')
        c=children(exp);
        if iscell(c)==1
            c=[c{:}];
        end
        a=1;
        b=1;
        for j=1:length(c)
            [aux1,aux2]=rational2polinomial(c(j));
            a=a*aux1;
            b=b*aux2;
        end
    elseif isSymType(exp,'plus')
        c=children(exp);
        if iscell(c)==1
            c=[c{:}];
        end
        a=0;
        b=1;
        for j=1:length(c)
            [aux1,aux2]=rational2polinomial(c(j));
            a=a*aux2+aux1*b;
            b=b*aux2;
        end
    else
        try
            c=children(exp);
            if iscell(c)==1
                c=[c{:}];
            end
            a=taylor(exp,c,0,'Order',100);
            b=1;
            warning('MATLAB:Taylor',['Found nonrational functions in' ...
                ' the equations and replaced them with their' ...
                ' Taylor expansion'])
            warning("off",'MATLAB:Taylor')
        catch
            try
                c=children(exp);
                if iscell(c)==1
                    c=[c{:}];
                 end
                a=taylor(exp,c,1,'Order',100);
                b=1;
                warning('MATLAB:Taylor',['Found nonrational functions' ...
                    ' in the equations and replaced them with their' ...
                    ' Taylor expansion'])
                warning("off",'MATLAB:Taylor')
            catch
                try
                    c=children(exp);
                    if iscell(c)==1
                        c=[c{:}];
                    end
                    syms cs
                    exp=subs(exp,c,cs);
                    a=taylor(exp,cs,0,'Order',100);
                    a=subs(a,cs,c);
                    b=1;
                    warning('MATLAB:Taylor',['Found nonrational functions' ...
                        ' in the equations and replaced them with their' ...
                        ' Taylor expansion'])
                    warning("off",'MATLAB:Taylor')
                catch
                    try
                        c=children(exp);
                        if iscell(c)==1
                            c=[c{:}];
                        end
                        syms cs
                        exp=subs(exp,c,cs);
                        a=taylor(exp,cs,1,'Order',100);
                        a=subs(a,cs,c);
                        b=1;
                        warning('MATLAB:Taylor',['Found nonrational functions' ...
                        ' in the equations and replaced them with their' ...
                        ' Taylor expansion'])
                        warning("off",'MATLAB:Taylor')
                    catch
                        error(['The model is not rational: it cannot be' ...
                            ' analysed with this algorithm. To use the FISPO' ...
                            ' algorithm instead, set in' ...
                            ' options.m: opts.algorithm = 1'])
                    end
                end
            end
        end
    end
else
    error('Argument must be just one equation in rational2polinomial')
end
