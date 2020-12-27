 %i=1 Hydrogen2, i=2 Oxygen2, , i=3 Nitrogen2,i=4 H20,i=5 OH,i=6 hydrogen, i=7 oxygen     

 function  [p4 T5 xi5 M5 u5 m5]=solveExit(T4,M4,p5,A5,alpha,precision)

MW=[2.01588 31.99880 28.01340 17.00734 18.01528 1.00794 15.99940];

Ru=8.314510;
H4=[0 0 0];
xi4=[0 0 0];
Kp=[0 0 0 0];
%calculate xi4
xi4(1)=1/(1+alpha+79/21*alpha);
xi4(2)=alpha/(1+alpha+79/21*alpha);
xi4(3)=79/21*alpha/(1+alpha+79/21*alpha);

MW4=xi4(1)*MW(1)+xi4(2)*MW(2)+xi4(3)*MW(3);
R4=Ru/MW4*1000;


a = importfile1('termodata.dat', [2,4,6,8,10,12,14],[2,4,6,8,10,12,14]);

% /*Previous Calculations*/
gamma4=computeGamma4(T4,xi4,a);
u4=M4*sqrt(gamma4*R4*T4);


%calculate enthalpy 4
for i=1:3
H4(i)=Ru*T4*computeEnthalpy(T4, a, i);
end
%A. Assume T5
T5=T4;

diff=1;
xi5=[0 0 0 0 0 0 0];

 while diff>precision

Kp(1)=-(computeGibbs(T5,a,4,xi5(1),p5)-computeGibbs(T5,a,1,xi5(3),p5)-0.5*computeGibbs(T5,a,2,xi5(4),p5));
Kp(2)=-(2*computeGibbs(T5,a,5,xi5(2),p5)-computeGibbs(T5,a,1,xi5(3),p5)-computeGibbs(T5,a,2,xi5(4),p5));
Kp(3)=-(2*computeGibbs(T5,a,6,xi5(5),p5)-computeGibbs(T5,a,1,xi5(3),p5));
Kp(4)=-(2*computeGibbs(T5,a,7,xi5(6),p5)-computeGibbs(T5,a,2,xi5(4),p5));


      
%calculate xi5
xi5=solveXiBetaEpsilon(p5,Kp, precision,alpha);
 MW5=xi5(1)*MW(4)+xi5(2)*MW(5)+xi5(3)*MW(1)+xi5(4)*MW(2)+xi5(5)*MW(6)+xi5(6)*MW(7)+xi5(7)*MW(3);
 R5=Ru/MW5*1000;
%calculate enthalpy 5
H5(1)= Ru*T5*computeEnthalpy(T5, a, 4); %H20
H5(2)= Ru*T5*computeEnthalpy(T5, a, 5); %OH
H5(3)= Ru*T5*computeEnthalpy(T5, a, 1); %H2
H5(4)= Ru*T5*computeEnthalpy(T5, a, 2); %O2
H5(5)= Ru*T5*computeEnthalpy(T5, a, 6); %H
H5(6)= Ru*T5*computeEnthalpy(T5, a, 7); %O
H5(7)= Ru*T5*computeEnthalpy(T5, a, 3); %N2
           
 %B.calculate u5 eq 10.13
 H4tot=0;
 H5tot=0;
 for i=1:3
     H4tot=H4tot+H4(i)*xi4(i);
 end
 for i=1:7
     H5tot=H5tot+H5(i)*xi5(i);
 end
 H5tot=H5tot/MW5;%*1000;
 H4tot=H4tot/MW4;%*1000;
 u5=sqrt(2*(H4tot+u4^2/2-H5tot));
 gamma5=computeGamma5(T5,xi5,a);
 M5=u5/sqrt(gamma5*R5*T5);
 %C. Find rho4 and p4

rho5=p5/R5/T5;
m5=rho5*u5*A5;
rho4=rho5*u5/u4;
p4=rho4*R4*T4;
%D. Calcualte p4* until p4=p4* iterating T5
p4aux=rho5*u5*(u5-u4)+p5;
 
diff=abs(p4aux-p4);

if p4aux<p4
  T5=T5+T5*diff/p4;%diff*10;%
else
  T5=T5-T5*diff/p4;%diff*10;%
end


 end


