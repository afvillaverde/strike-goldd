%--------------------------------------------------------------------------
% File that defines the CHO model.
% It stores the variables in a mat-file called CHO_115p.mat 
% The model is a modification of the B4 model from the BioPreDyn-bench
% collection. It is taken from:
%--------------------------------------------------------------------------
% Villaverde AF et al (2015) BioPreDyn-bench: a suite of benchmark problems 
% for dynamic modelling in systems biology. BMC systems biology, 9(1), 1.
%--------------------------------------------------------------------------
% Changes w.r.t. B4: the following parameters are grouped:
% par47 - par48 --> (new) par47
% -par55 + par57 --> (new) par55
%--------------------------------------------------------------------------
clear all;

% 34 states, x1, x2, ..., x(34)
x = sym('x', [34 1]);

% 13 outputs:
h = x([5 4 3 2 1 29 27 21 15 13 30 32 11]);

% 117 parameters, par1, par2, ..., par117
p = sym('par', [117 1]);

% inputs:
u    = [];

% constants:
rho = 141.4710605D0;
vol_1 = 0.9D0;
vol_2 = 0.1D0;
rss6 = 2413.6562849217303D0*vol_1*1.0/vol_2;
rss7 = 603.414071180172D0;
rss8 = 40133.38664791273D0*vol_1*1.0/vol_2;
rss9 = 74836.04665510054D0*vol_1*1.0/vol_2;
rss10 = 180754.26889005644D0*vol_1*1.0/vol_2;
rss11 = 164390.2745505474D0;
rss12 = 36512.902220891134D0;
rss13 = 166803.93083553354D0;
rss14 = 39529.972576829896D0;
rss15 = 1810.242213701644D0;
rss16 = 127877.37232994902D0;
rss17 = 166200.51676427448D0;
rss18 = 216060.34296802033D0*vol_1*1.0/vol_2;
rss19 = 5.028450593613567D0;
rss20 = 35306.0740779408D0*vol_1*1.0/vol_2;
rss21 = 0.0D0*vol_1*1.0/vol_2;
rss22 = 532531.1670427525D0;
rss23 = 38926.55850561767D0;
rss24 = 40133.38664808034D0;
rss25 = 689982.1797742102D0;
rss26 = 495414.85075147304D0*vol_1*1.0/vol_2;
rss27 = 2413.6562849329107D0;
rss28 = 83401.9654177705D0;
rss29 = 603.4140712332287D0;
rss30 = 603.4140712332287D0;
rss31 = 2413.656284931604D0;
rss32 = 533134.5811139438D0;
css1=100000;
css2=1;
css3=1000;
css4=1000;
css5=30;
css6=1000;
css7=1000;
css8=1000;
css9=1000;
css10=1000;
css11=1000;
css12=1000;
css13=3000;
css14=1000;
css15=1000;
css16=1000;
css17=100;
css18=100;
css19=100;
css20=100;
css21=1000;
css22=1000;
css23=1000;
css24=1000;
css25=1000;
css26=1000;
css27=1000;
css28=1000;
css29=1000;
css30=3000;
css31=100000;
css32=1000;
css33=1000;
css34=1000;

% elasticities:
e1=p(1);
e2=p(2);
e3=p(3);
e4=p(4);
e5=p(5);
e6=p(6);
e7=p(7);
e8=p(8);
e14=p(9);
e15=p(10);
e16 = -p(11);
e17 = -p(12);
e18=p(13);
e19=p(14);
e20 = -p(15);
e21 = -p(16);
e22=p(17);
e23=p(18);
e24 = -p(19);
e25 = -p(20);
e26=p(21);
e27=p(22);
e28=p(23);
e29 = -p(24);
e30 = -p(25);
e31=p(26);
e32=p(27);
e33=p(28);
e34 = -p(29);
e35 = -p(30);
e36 = -p(31);
e37=p(32);
e38=p(33);
e39 = -p(34);
e40 = -p(35);
e41 = 0.0D0;
e42=p(36);
e43=p(37);
e44=p(38);
e45 = -p(39);
e46 = -p(40);
e47=p(41);
e48 = 0.0D0;
e49 = -p(42);
e50=p(43);
e51=p(44);
e52 = -p(45);
e53 = -p(46);
e54=p(47);
e55 = -p(48);
e56 = -p(49);
e57 = -p(50);
e58=p(51);
e59=p(52);
e60=p(53);
e61 = -p(54);
e62 = -p(55);
e63 = -p(56);
e64=p(57);
e65 = -p(58);
e66=p(59);
e67=p(60);
e68=p(61);
e69 = -p(62);
e70 = -p(63);
e71 = -p(64);
e72=p(65);
e73=p(66);
e74 = -p(67);
e75 = -p(68);
e76=p(69);
e77 = -p(70);
e78=p(71);
e79=p(72);
e80 = -p(73);
e81 = -p(74);
e82=p(75);
e83=p(76);
e84=p(77);
e85=p(78);
e86=p(79);
e87 = -p(80);
e88 = -p(81);
e89 = -p(82);
e90 = -p(83);
e91 = 0.0D0;
e92 = 0.0D0;
e93 = -p(84);
e94 = 0.7D0;
e95 = -2.0D-1;
e96=p(85);
e97=p(86);
e98 = -p(87);
e99 = -p(88);
e100=p(89);
e101=p(90);
e102 = -p(91);
e103 = -p(92);
e104=p(93);
e105=p(94);
e106 = -p(95);
e107 = -p(96);
e108=p(97);
e109 = -p(98);
e110=p(99);
e111=p(100);
e112=p(101);
e113 = -p(102);
e114 = -p(103);
e115=p(104);
e116=p(105);
e117 = -p(106);
e118=p(107);
e119 = -p(108);
e120=p(109);
e121 = -p(110);
e122=p(111);
e123 = -p(112);
e124=p(113);
e125 = -p(114);
e126=p(115);
e127 = -p(116);
e128 = -p(117);

% auxiliary variables:
clogk1=log(x(1));
clogk2=log(x(2));
clogk3=log(x(3));
clogk4=log(x(4));
clogk5=log(x(5));
clogk6=log(x(6));
clogk7=log(x(7));
clogk8=log(x(8));
clogk9=log(x(9));
clogk10=log(x(10));
clogk11=log(x(11));
clogk12=log(x(12));
clogk13=log(x(13));
clogk14=log(x(14));
clogk15=log(x(15));
clogk16=log(x(16));
clogk17=log(x(17));
clogk18=log(x(18));
clogk19=log(x(19));
clogk20=log(x(20));
clogk21=log(x(21));
clogk22=log(x(22));
clogk23=log(x(23));
clogk24=log(x(24));
clogk25=log(x(25));
clogk26=log(x(26));
clogk27=log(x(27));
clogk28=log(x(28));
clogk29=log(x(29));
clogk30=log(x(30));
clogk31=log(x(31));
clogk32=log(x(32));
clogk33=log(x(33));
clogk34=log(x(34));

% rates:
r1=0.0D0;
r2=0.0D0;
r3=0.0D0;
r4=0.0D0;
r5=0.0D0;
r6=((rss6*(1.0D0+(((e14*clogk8)+(e15*clogk9))+((e16*clogk6)+(e17*clogk7))))));
r7=((rss7*(1.0D0+(((e18*clogk12)+(e19*clogk13))+((e20*clogk10)+(e21*clogk11))))));
r8=((rss8*(1.0D0+((((e22*clogk15)+(e23*clogk6))+((e24*clogk8)+(e25*clogk14)))+(e26*clogk16)))));
r9=((rss9*(1.0D0+(((e27*clogk16)+(e28*clogk7))+((e29*clogk15)+(e30*clogk9))))));
r10=((rss10*(1.0D0+((((e31*clogk9)+(e32*clogk19))+(e33*clogk20))+(((e34*clogk7)+(e35*clogk17))+(e36*clogk18))))));
r11=((rss11*(1.0D0+((((e37*clogk22)+(e38*clogk11))+((e39*clogk21)+(e40*clogk13)))+(e41*clogk13)))));
r12=((rss12*(1.0D0+(((((((e42*clogk15)+(e43*clogk7))+(e44*clogk21))+((e45*clogk8)+(e46*clogk9)))+(e47*clogk32))+(e48*clogk9))+(e49*clogk30)))));
r13=((rss13*(1.0D0+(((((((e50*clogk25)+(e51*clogk26))+((e52*clogk23)+(e53*clogk24)))+(e54*clogk11))+(e55*clogk11))+(e56*clogk13))+(e57*clogk22)))));
r14=((rss14*(1.0D0+((((((e58*clogk28)+(e59*clogk29))+(e60*clogk23))+(((e61*clogk12)+(e62*clogk27))+(e63*clogk25)))+(e64*clogk27))+(e65*clogk10)))));
r15=((rss15*(1.0D0+((((e66*clogk32)+(e67*clogk22))+(e68*clogk16))+(((e69*clogk15)+(e70*clogk30))+(e71*clogk27))))));
r16=((rss16*(1.0D0+(((e72*clogk21)+(e73*clogk23))+((e74*clogk2)+(e75*clogk25))))));
r17=((rss17*(1.0D0+((e76*clogk24)+(e77*clogk22)))));
r18=((rss18*(1.0D0+(((e78*clogk20)+(e79*clogk17))+((e80*clogk18)+(e81*clogk19))))));
r19=((rss19/(((((((((css24/e1)*(css25/e2))*(css12/e3))*(css33/e4))*(css34/e5))*(css29/e6))*(css10/e7))*(css13/e8))/((((((((1.0D0+(css24/e1))*(1.0D0+(css25/e2)))*(1.0D0+(css12/e3)))*(1.0D0+(css33/e4)))*(1.0D0+(css34/e5)))*(1.0D0+(css29/e6)))*(1.0D0+(css10/e7)))*(1.0D0+(css13/e8)))))*((((((((((x(24)*css24)/e1)*((x(25)*css25)/e2))*((x(12)*css12)/e3))*((x(33)*css33)/e4))*((x(34)*css34)/e5))*((x(29)*css29)/e6))*((x(10)*css10)/e7))*((x(13)*css13)/e8))/((((((((1.0D0+((x(24)*css24)/e1))*(1.0D0+((x(25)*css25)/e2)))*(1.0D0+((x(12)*css12)/e3)))*(1.0D0+((x(33)*css33)/e4)))*(1.0D0+((x(34)*css34)/e5)))*(1.0D0+((x(29)*css29)/e6)))*(1.0D0+((x(10)*css10)/e7)))*(1.0D0+((x(13)*css13)/e8)))));
r20=((rss20*(1.0D0+(((((((((e82*clogk8)+(e83*clogk7))+(e84*clogk19))+(e85*clogk31))+(e86*clogk32))+((((e87*clogk16)+(e88*clogk9))+(e89*clogk17))+(e90*clogk30)))+(e91*clogk31))+(e92*clogk9))+(e93*clogk15)))));
r21=((rss21*(1.0D0+((e94*clogk18)+(e95*clogk20)))));
r22=((rss22*(1.0D0+(((e96*clogk11)+(e97*clogk30))+((e98*clogk32)+(e99*clogk13))))));
r23=((rss23*(1.0D0+(((e100*clogk8)+(e101*clogk27))+((e102*clogk28)+(e103*clogk16))))));
r24=((rss24*(1.0D0+(((e104*clogk14)+(e105*clogk12))+((e106*clogk29)+(e107*clogk6))))));
r25=((rss25*(1.0D0+((e108*clogk13)+(e109*clogk11)))));
r26=((rss26*(1.0D0+((((e110*clogk32)+(e111*clogk31))+(e112*clogk18))+((e113*clogk30)+(e114*clogk20))))));
r27=((rss27*(1.0D0+(((e115*clogk27)+(e116*clogk31))+(e117*clogk16)))));
r28=((rss28*(1.0D0+((e118*clogk1)+(e119*clogk26)))));
r29=((rss29*(1.0D0+((e120*clogk3)+(e121*clogk33)))));
r30=((rss30*(1.0D0+((e122*clogk4)+(e123*clogk34)))));
r31=((rss31*(1.0D0+((e124*clogk6)+(e125*clogk12)))));
r32=((rss32*(1.0D0+((e126*clogk18)+((e127*clogk31)+(e128*clogk20))))));

% dynamic equations:
f = sym('f');
f(1)=(r1 -r28*1.0D0/rho)/css1;
f(2)=(-r2+r16*1.0D0/rho)/css2;
f(3)=(r3 -r29*1.0D0/rho)/css3;
f(4)=(r4 -r30*1.0D0/rho)/css4;
f(5)=(-r5+r19*1.0D0/rho)/css5;
f(6)=(r6-r8+r24*vol_1/vol_2-r31*vol_1/vol_2)/css6;
f(7)=(r6-r9+r10-2.0D0*r12*vol_1/vol_2-r20)/css7;
f(8)=(-r6+r8+r12*vol_1/vol_2-r20-r23*vol_1/vol_2)/css8;
f(9)=(-r6+r9-r10+2.0D0*r12*vol_1/vol_2+r20)/css9;
f(10)=(r7-120.0D0*r19)/css10;
f(11)=(r7-r11+1260.0D0*r19-r22+r25)/css11;
f(12)=(-r7+r14-240.0D0*r19-r24+r31)/css12;
f(13)=(-r7+r11-1260.0D0*r19+r22-r25)/css13;
f(14)=(r8-r24*vol_1/vol_2)/css14;
f(15)=(-r8+r9-r12*vol_1/vol_2+r15*vol_1/vol_2)/css15;
f(16)=(-r9-r15*vol_1/vol_2+r20+r23*vol_1/vol_2+r27*vol_1/vol_2)/css16;
f(17)=(2.0D0*r10-2.0D0*r18+2.0D0*r20)/css17;
f(18)=(4.0D0*r10+6.0D0*r18-r21-3.0D0*r26-r32*vol_1/vol_2)/css18;
f(19)=(-2.0D0*r10+2.0D0*r18-2.0D0*r20)/css19;
f(20)=(-4.0D0*r10-6.0D0*r18+r21+3.0D0*r26+r32*vol_1/vol_2)/css20;
f(21)=(r11-r12-r16)/css21;
f(22)=(-r11-r15+r17)/css22;
f(23)=(r13-r14-r16+120.0D0*r19)/css23;
f(24)=(r13-r17-120.0D0*r19)/css24;
f(25)=(-r13+r14+r16-120.0D0*r19)/css25;
f(26)=(-0.5D0*r13+r28)/css26;
f(27)=(r14+r15-r23-r27)/css27;
f(28)=(-r14+120.0D0*r19+r23)/css28;
f(29)=(-r14-120.0D0*r19+r24)/css29;
f(30)=(r15*vol_1/vol_2+r20-r22*vol_1/vol_2+r26)/css30;
f(31)=(-r20-r26-r27*vol_1/vol_2+r32*vol_1/vol_2)/css31;
f(32)=(-r15*vol_1/vol_2-r20+r22*vol_1/vol_2-r26)/css32;
f(33)=(-120.0D0*r19+r29)/css33;
f(34)=(-120.0D0*r19+r30)/css34; 
f = f.';

%%% Group the unidentifiable parameters:
subs(f,p(47)-p(48),p(47));
subs(f,-p(55)+p(57),p(55));
p = p([1:47,49:56,58:end]);
%%%

% initial conditions:
ics = [];

% which initial conditions are known:
known_ics = zeros(1,numel(x));  

save('CHO_115p','x','h','u','p','f','ics','known_ics');
