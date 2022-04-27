% computation of the homogeneous resolution
function HomSol=calcHomSol(logder_esp,Order,n,Myprime)

% integration of polynomials in a matrix
HomSol=-intmatpoly(logder_esp,Order);


HomSol=ExpMatrixSeries(HomSol,Order,n,Myprime);