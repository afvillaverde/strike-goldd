function CreateNewModel()
clear all
syms  x1 x2 x3
x = transpose([x1 x2 x3]);
syms p1 p2 p3 p4 p5
p = transpose([p1 p2 p3 p4 p5]);
h = transpose([x3, x2+x1]);
syms u1
u = [u1];
w = [];f = [ p1 - p2*x1 - u1*x1*x3;
 u1*x1*x3 - p4*x2;
p3*p4*x2 - p5*x3];
ics  = [0 0 0];
known_ics = [1 1 1 ];
save("C:\Users\xabor\OneDrive - Universidade de Vigo\Escritorio\STRIKE-GOLDD-Ultima-version\strike-goldd-master\STRIKE-GOLDD\models\EXAMPLE","x","p","u","w","h","f","ics","known_ics");
end