function [T4]=solveInlet(M4,T0,M0,alpha)
MW=[2.01588 31.99880 28.01340 17.00734 18.01528 1.00794 15.99940];
Ru=8.314510;
diff=1;
H4tot=0;
xi0(1)=0.21;
xi0(2)=0.79;
MW0=xi0(2)*MW(3)+MW(2)*xi0(1);
xi4=[0 0 0];
xi4(1)=1/(1+alpha+79/21*alpha);
xi4(2)=alpha/(1+alpha+79/21*alpha);
xi4(3)=79/21*alpha/(1+alpha+79/21*alpha);
a = importfile1('termodata.dat', [2,4,6,8,10,12,14],[2,4,6,8,10,12,14]);
MW4=xi4(1)*MW(1)+xi4(2)*MW(2)+xi4(3)*MW(3);
R4=Ru/MW4*1000;

u0=M0*sqrt(computeGamma0(T0,xi0,a)*287*T0);
T4=10;
H0(1)=Ru*T0*computeEnthalpy(T0, a, 2);
H0(2)=Ru*T0*computeEnthalpy(T0, a, 3);
H0tot=xi0(1)*H0(1)+xi0(2)*H0(2);
H0tot/MW0;%*1000;
Hto=H0tot+0.5*u0^2;

while abs(diff) >0.01
u4=M4*sqrt(computeGamma4(T4,xi4,a)*R4*T4);
%u4=sqrt(0.97)*u0;
H4=Hto-0.5*u4^2;

for i=1:3
H4aux(i)=Ru*T4*computeEnthalpy(T4, a, i);
H4totaux=H4tot+H4aux(i)*xi4(i);
end
H4totaux=H4totaux/MW4;%*1000;
diff=abs(H4totaux-H4);
if H4totaux>H4
    T4=T4-diff/100000;
else
    T4=T4+diff/100000;
end
end

end