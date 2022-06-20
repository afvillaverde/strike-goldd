%--------------------------------------------------------------------------
% function that fixes the form of SLPOut to obtain the observability -
% identifiability matrix
%--------------------------------------------------------------------------

function onx=constructonx(SLPOut,OldOrder)

[nr,nc]=size(SLPOut);
onx=zeros(nr*nc,OldOrder);
 for i=1:nr
    for j=1:nc
        if SLPOut(i,j) == 0
            aux=zeros(1,OldOrder);
        else
            aux=coeffs(SLPOut(i,j),'t','All');
        end
        onx(nr*length(aux)-nr+i:-nr:i,j)=aux;
    end
 end