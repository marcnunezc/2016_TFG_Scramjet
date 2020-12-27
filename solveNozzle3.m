function [M9 T9 u9 A9]=solveNozzle3(T5,xi5,p0,p5,M5,T9_1,A5,alpha,precision)
S5tot=0;
H9tot=0;
H5tot=0;
diff=1;
dif=1;
a = importfile1('termodata.dat', [2,4,6,8,10,12,14],[2,4,6,8,10,12,14]);

MW=[2.01588 31.99880 28.01340 17.00734 18.01528 1.00794 15.99940];
MW5=xi5(1)*MW(4)+xi5(2)*MW(5)+xi5(3)*MW(1)+xi5(4)*MW(2)+xi5(5)*MW(6)+xi5(6)*MW(7)+xi5(7)*MW(3);
Ru=8.314510;
R5=Ru/MW5*1000;
gamma5=computeGamma5(T5,xi5,a);
u5=M5*sqrt(gamma5*R5*T5);


S5(1)= Ru*computeEntropy(T5, a, 4); %H20
S5(2)= Ru*computeEntropy(T5, a, 5); %OH
S5(3)= Ru*computeEntropy(T5, a, 1); %H2
S5(4)= Ru*computeEntropy(T5, a, 2); %O2
S5(5)= Ru*computeEntropy(T5, a, 6); %H
S5(6)= Ru*computeEntropy(T5, a, 7); %O
S5(7)= Ru*computeEntropy(T5, a, 3); %N2

H5(1)= Ru*T5*computeEnthalpy(T5, a, 4); %H20
H5(2)= Ru*T5*computeEnthalpy(T5, a, 5); %OH
H5(3)= Ru*T5*computeEnthalpy(T5, a, 1); %H2
H5(4)= Ru*T5*computeEnthalpy(T5, a, 2); %O2
H5(5)= Ru*T5*computeEnthalpy(T5, a, 6); %H
H5(6)= Ru*T5*computeEnthalpy(T5, a, 7); %O
H5(7)= Ru*T5*computeEnthalpy(T5, a, 3); %N2   

 for i=1:7
     H5tot=H5tot+H5(i)*xi5(i);
     S5tot=S5tot+S5(i)*xi5(i);
 end
 H5tot=H5tot/MW5*1000;

T9old=1;
xi9=[0 0 0 0 0 0 0];
while diff >0.0001
Kp(1)=-(computeGibbs(T9old,a,4,xi9(1),p5)-computeGibbs(T9old,a,1,xi9(3),p5)-0.5*computeGibbs(T9old,a,2,xi9(4),p5));
Kp(2)=-(2*computeGibbs(T9old,a,5,xi9(2),p5)-computeGibbs(T9old,a,1,xi9(3),p5)-computeGibbs(T9old,a,2,xi9(4),p5));
Kp(3)=-(2*computeGibbs(T9old,a,6,xi9(5),p5)-computeGibbs(T9old,a,1,xi9(3),p5));
Kp(4)=-(2*computeGibbs(T9old,a,7,xi9(6),p5)-computeGibbs(T9old,a,2,xi9(4),p5));
    
xi9=solveXiBetaEpsilon(p0,Kp, precision,alpha);
MW9=xi9(1)*MW(4)+xi9(2)*MW(5)+xi9(3)*MW(1)+xi9(4)*MW(2)+xi9(5)*MW(6)+xi9(6)*MW(7)+xi9(7)*MW(3); 
S9tot=S5tot-Ru*log(p5)+Ru*log(p0);
x=T9old;

while abs(dif)>0.0001
dif=S9tot-xi9(1)*Ru*(a(4,1)*log(x)+a(4,2)*x+a(4,3)*x^2/2+a(4,4)*x^3/3+a(4,5)*x^4/4+a(4,7))...
    -xi9(2)*Ru*(a(5,1)*log(x)+a(5,2)*x+a(5,3)*x^2/2+a(5,4)*x^3/3+a(5,5)*x^4/4+a(5,7))...
    -xi9(3)*Ru*(a(1,1)*log(x)+a(1,2)*x+a(1,3)*x^2/2+a(1,4)*x^3/3+a(1,5)*x^4/4+a(1,7))...
    -xi9(4)*Ru*(a(2,1)*log(x)+a(2,2)*x+a(2,3)*x^2/2+a(2,4)*x^3/3+a(2,5)*x^4/4+a(2,7))...
    -xi9(5)*Ru*(a(6,1)*log(x)+a(6,2)*x+a(6,3)*x^2/2+a(6,4)*x^3/3+a(6,5)*x^4/4+a(6,7))...
    -xi9(6)*Ru*(a(7,1)*log(x)+a(7,2)*x+a(7,3)*x^2/2+a(7,4)*x^3/3+a(7,5)*x^4/4+a(7,7))...
    -xi9(7)*Ru*(a(3,1)*log(x)+a(3,2)*x+a(3,3)*x^2/2+a(3,4)*x^3/3+a(3,5)*x^4/4+a(3,7));
if dif>0
    x=x+0.1*x*dif/S9tot;
else
    x=x-0.1*x*dif/S9tot;
end
end
T9=x;
 diff=abs(T9-T9old);
 T9old=T9;
end
H9(1)= Ru*T9*computeEnthalpy(T9, a, 4); %H20
H9(2)= Ru*T9*computeEnthalpy(T9, a, 5); %OH
H9(3)= Ru*T9*computeEnthalpy(T9, a, 1); %H2
H9(4)= Ru*T9*computeEnthalpy(T9, a, 2); %O2
H9(5)= Ru*T9*computeEnthalpy(T9, a, 6); %H
H9(6)= Ru*T9*computeEnthalpy(T9, a, 7); %O
H9(7)= Ru*T9*computeEnthalpy(T9, a, 3); %N2          
 for i=1:7
     H9tot=H9tot+H9(i)*xi5(i);
 end

H9tot=H9tot/MW9*1000;
R9=Ru/MW9*1000;
rho5=p5/R5/T5;
rho9=p0/R9/T9;
u9=sqrt(u5^2+2*(abs(H9tot-H5tot)));
A9=rho5*u5*A5/rho9/u9;

M9=u9/sqrt(computeGamma5(T9,xi9,a)*R9*T9);
end
