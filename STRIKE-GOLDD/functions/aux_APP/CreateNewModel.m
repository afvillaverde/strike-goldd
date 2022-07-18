function CreateNewModel()
clear all
syms Tu Ti V
x = transpose([Tu Ti V]);
syms lambda rho N delta c
p = transpose([lambda rho N delta c]);
h = transpose([V, Tu+Ti]);
syms eta
u = [eta];
 
w = [];
f = [ lambda-rho*Tu-eta*Tu*V;
 eta*Tu*V-delta*Ti;
N*delta*Ti-c*V];
ics  = [[000]];
known_ics = [1 ];
save("C:\Users\xabor\Downloads\APP\strike-goldd-dev\STRIKE-GOLDD\models\HIV","x","p","u","w","h","f","ics","known_ics");
end