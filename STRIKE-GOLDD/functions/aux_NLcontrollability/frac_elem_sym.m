%--------------------------------------------------------------------------
% Auxiliar function that changes a rational equation into its numerator and
% denominator
%--------------------------------------------------------------------------

function [num,den]=frac_elem_sym(exp)

lexp=length(exp);
if lexp==1
    if isSymType(exp,'number') || isSymType(exp,'variable')
        num=exp;
        den=1;
    elseif isSymType(exp,'power')
        c=children(exp);
        if iscell(c)==1
            c=[c{:}];
        end
        try 
            if c(2)<0
                [aux1,aux2]=frac_elem_sym(c(1));
                num=aux2^(-c(2));
                den=aux1^(-c(2));
            elseif c(2)==0
                num=1;
                den=1;
            else
                [aux1,aux2]=frac_elem_sym(c(1));
                num=aux1^(c(2));
                den=aux2^(c(2));
            end
        catch
            error(['The model is not rational: it cannot be analysed' ...
                ' with this algorithm.'])
        end
    elseif isSymType(exp,'times')
        c=children(exp);
        if iscell(c)==1
            c=[c{:}];
        end
        num=1;
        den=1;
        for j=1:length(c)
            [aux1,aux2]=frac_elem_sym(c(j));
            num=num*aux1;
            den=den*aux2;
        end
    elseif isSymType(exp,'plus')
        c=children(exp);
        if iscell(c)==1
            c=[c{:}];
        end
        num=0;
        den=1;
        for j=1:length(c)
            [aux1,aux2]=frac_elem_sym(c(j));
            num=num*aux2+aux1*den;
            den=den*aux2;
        end
    else
%         try
%             c=children(exp);
%             if iscell(c)==1
%                 c=[c{:}];
%             end
%             num=taylor(exp,c,0,'Order',100);
%             den=1;
%             warning('MATLAB:Taylor',['Found nonrational functions in' ...
%                 ' the equations and replaced them with their' ...
%                 ' Taylor expansion'])
%             warning("off",'MATLAB:Taylor')
%         catch
%             try
%                 c=children(exp);
%                 if iscell(c)==1
%                     c=[c{:}];
%                  end
%                 num=taylor(exp,c,1,'Order',100);
%                 den=1;
%                 warning('MATLAB:Taylor',['Found nonrational functions' ...
%                     ' in the equations and replaced them with their' ...
%                     ' Taylor expansion'])
%                 warning("off",'MATLAB:Taylor')
%             catch
%                 try
%                     c=children(exp);
%                     if iscell(c)==1
%                         c=[c{:}];
%                     end
%                     syms cs
%                     exp=subs(exp,c,cs);
%                     num=taylor(exp,cs,0,'Order',100);
%                     num=subs(num,cs,c);
%                     den=1;
%                     warning('MATLAB:Taylor',['Found nonrational functions' ...
%                         ' in the equations and replaced them with their' ...
%                         ' Taylor expansion'])
%                     warning("off",'MATLAB:Taylor')
%                 catch
%                     try
%                         c=children(exp);
%                         if iscell(c)==1
%                             c=[c{:}];
%                         end
%                         syms cs
%                         exp=subs(exp,c,cs);
%                         num=taylor(exp,cs,1,'Order',100);
%                         num=subs(num,cs,c);
%                         den=1;
%                         warning('MATLAB:Taylor',['Found nonrational functions' ...
%                         ' in the equations and replaced them with their' ...
%                         ' Taylor expansion'])
%                         warning("off",'MATLAB:Taylor')
%                     catch
%                         error(['The model is not rational: it cannot be' ...
%                             ' analysed with this algorithm. To use the FISPO' ...
%                             ' algorithm instead, set in' ...
%                             ' options.m: opts.algorithm = 1'])
%                     end
%                 end
%             end
%         end
    end
else
    error('Argument must be just one equation in frac_elem')
end