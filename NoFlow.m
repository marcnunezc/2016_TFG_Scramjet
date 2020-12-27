%i=1 Hydrogen2, i=2 Oxygen2, , i=3 Nitrogen2,i=4 H20,i=5 OH,i=6 hydrogen, i=7 oxygen     
 clc
clear all

MW=[2.01588 31.99880 28.01340 17.00734 18.01528 1.00794 15.99940];
N=20;
p5=1;
for alpha=1:10
    if alpha==0
        alpha=0.1
    end
    alpha=alpha/2
Ru=8.314510;
H4=[0 0 0];
xi4=[0 0 0];
Kp=[0 0 0 0];
m=1;

a = importfile1('termodata.dat', [2,4,6,8,10,12,14],[2,4,6,8,10,12,14]);
precision=0.00001;
% /*Previous Calculations*/
for T4=linspace(0,2500,N)
 T4vec(m)=T4;
u4=0;
%u4=3*sqrt(1.4*287*T4);
%calculate xi4
    xi4(1)=1/(1+alpha+79/21*alpha);
xi4(2)=alpha/(1+alpha+79/21*alpha);
xi4(3)=79/21*alpha/(1+alpha+79/21*alpha);
xi5=[0 0 0 0 0 0 0];
        %calculate enthalpy 4
      
        for i=1:3
        H4(i)=Ru*T4*computeEnthalpy(T4, a, i);
        end
%A. Assume T5
T5=819; %kelvin
diff=1;
%compute KP
while diff>0.0001
diff=1;




Kp(1)=-(computeGibbs(T5,a,4,xi5(1),p5)-computeGibbs(T5,a,1,xi5(3),p5)-0.5*computeGibbs(T5,a,2,xi5(4),p5));
Kp(2)=-(2*computeGibbs(T5,a,5,xi5(2),p5)-computeGibbs(T5,a,1,xi5(3),p5)-computeGibbs(T5,a,2,xi5(4),p5));
Kp(3)=-(2*computeGibbs(T5,a,6,xi5(5),p5)-computeGibbs(T5,a,1,xi5(3),p5));
Kp(4)=-(2*computeGibbs(T5,a,7,xi5(6),p5)-computeGibbs(T5,a,2,xi5(4),p5));

p5=1; %atm

      
%calculate xi5
xi5=solveXiBetaEpsilon(p5,Kp, precision,alpha);
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

 %C. Find rho4 and p4
 MW5=xi5(1)*MW(4)+xi5(2)*MW(5)+xi5(3)*MW(1)+xi5(4)*MW(2)+xi5(5)*MW(6)+xi5(6)*MW(7)+xi5(7)*MW(3);
 MW4=xi4(1)*MW(1)+xi4(2)*MW(2)+xi4(3)*MW(3);
 R4=Ru/MW4;
 R5=Ru/MW5;

  diff=H4tot-H5tot;
if alpha<=0.5
 if H4tot<H5tot
      T5=T5-abs(diff)/10;
    else
      T5=T5+abs(diff)/10;
    end
else
    if H4tot<H5tot
      T5=T5-abs(diff)/100;
    else
      T5=T5+abs(diff)/100;
    end
end

end
T5vec(m)=T5;
m=m+1;
end

hold on
plot(T5vec,T5vec-T4vec)
end 
m=1;
for T4=linspace(500,2500,5)
    
 for   alpha=linspace(1,10,N)
    alpha=alpha/2;
Ru=8.314510;
H4=[0 0 0];
xi4=[0 0 0];
Kp=[0 0 0 0];


a = importfile1('termodata.dat', [2,4,6,8,10,12,14],[2,4,6,8,10,12,14]);
precision=0.00001;
% /*Previous Calculations*/


u4=0;

%calculate xi4
    xi4(1)=1/(1+alpha+79/21*alpha);
xi4(2)=alpha/(1+alpha+79/21*alpha);
xi4(3)=79/21*alpha/(1+alpha+79/21*alpha);
xi5=[0 0 0 0 0 0 0];
        %calculate enthalpy 4
      
        for i=1:3
        H4(i)=Ru*T4*computeEnthalpy(T4, a, i);
        end
%A. Assume T5
T5=819; %kelvin
diff=1;
%compute KP
while diff>0.001
diff=1;




Kp(1)=-(computeGibbs(T5,a,4,xi5(1),p5)-computeGibbs(T5,a,1,xi5(3),p5)-0.5*computeGibbs(T5,a,2,xi5(4),p5));
Kp(2)=-(2*computeGibbs(T5,a,5,xi5(2),p5)-computeGibbs(T5,a,1,xi5(3),p5)-computeGibbs(T5,a,2,xi5(4),p5));
Kp(3)=-(2*computeGibbs(T5,a,6,xi5(5),p5)-computeGibbs(T5,a,1,xi5(3),p5));
Kp(4)=-(2*computeGibbs(T5,a,7,xi5(6),p5)-computeGibbs(T5,a,2,xi5(4),p5));

p5=1; %atm

      
%calculate xi5
xi5=solveXiBetaEpsilon(p5,Kp, precision,alpha);
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

 %C. Find rho4 and p4
 MW5=xi5(1)*MW(4)+xi5(2)*MW(5)+xi5(3)*MW(1)+xi5(4)*MW(2)+xi5(5)*MW(6)+xi5(6)*MW(7)+xi5(7)*MW(3);
 MW4=xi4(1)*MW(1)+xi4(2)*MW(2)+xi4(3)*MW(3);
 R4=Ru/MW4;
 R5=Ru/MW5;

  diff=H4tot-H5tot;
if alpha<=0.5
 if H4tot<H5tot
      T5=T5-abs(diff)/10
    else
      T5=T5+abs(diff)/10
    end
else
    if H4tot<H5tot
      T5=T5-abs(diff)/1000
    else
      T5=T5+abs(diff)/1000
    end
end
end
T5vec1(m)=T5;
m=m+1;
end

plot(T5vec1,T5vec1-T4)
end 
ylim([0 1800])
xlim([500 3000])