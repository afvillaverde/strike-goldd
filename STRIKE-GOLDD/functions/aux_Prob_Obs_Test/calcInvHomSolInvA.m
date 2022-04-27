% computation of calcInvHomSolInvA
function InvHomSolInvA=calcInvHomSolInvA(HomSol,invsep_esp,Order,n,Myprime,t)

InvHomSolInvA=InverseMatrixSeries(HomSol,Order,n,Myprime);

InvHomSolInvA=mod(multmatpolytrun(InvHomSolInvA,invsep_esp,Order,n),Myprime);