function [P4]=solveInletPerfect(M4,T0,T4,M0,P0,alpha)
MW=[2.01588 31.99880 28.01340 17.00734 18.01528 1.00794 15.99940];
Ru=8.314510;
eta_d=0.97;
Cptot=0;
a = importfile1('termodata.dat', [2,4,6,8,10,12,14],[2,4,6,8,10,12,14]);
diff=1;
xi0(1)=0.21;
xi0(2)=0.79;
xi4=[0 0 0];
xi4(1)=1/(1+alpha+79/21*alpha);
xi4(2)=alpha/(1+alpha+79/21*alpha);
xi4(3)=79/21*alpha/(1+alpha+79/21*alpha);
a = importfile1('termodata.dat', [2,4,6,8,10,12,14],[2,4,6,8,10,12,14]);
MW4=xi4(1)*MW(1)+xi4(2)*MW(2)+xi4(3)*MW(3);
R4=Ru/MW4*1000;

gamma0=computeGamma0(T0,xi0,a);
u0=M0*sqrt(computeGamma0(T0,xi0,a)*287*T0);
pi_d=(1+(1-eta_d)*0.5*(gamma0-1)*M0^2)^(-gamma0/(gamma0-1));
Pt0=P0*((1+(gamma0-1)/2*M0^2)^(gamma0/(gamma0-1)));
Pt4=Pt0*pi_d;

% Tt0=T0*(1+(gamma0-1)/2*M0^2);
% T4old=1000;
% while diff>0.0001
% gamma4=computeGamma4(T4old,xi4,a);
% T4=Tt0/(1+(gamma4-1)/2*M4^2);
% diff=abs(T4-T4old);
% T4old=T4*0.4+T4old*0.6;
% end
% T4old=2000;
% while diff >0.001
%         Cp(1)= Ru*computeCp(T4old, a, 1); %H2
%         Cp(2)= Ru*computeCp(T4old, a, 2); %O2
%         Cp(3)= Ru*computeCp(T4old, a, 3); %N2
%         for i=1:length(xi4)
%             Cptot=Cptot+Cp(i)*xi4(i);
%         end
%         Cptot=Cptot/MW4*1000;
%         T4=T0+u0^2/(2*Cptot)*(1-eta_d);
%         diff=abs(T4-T4old);
%         T4old=T4;
% end
gamma4=computeGamma4(T4,xi4,a);
P4=Pt4/((1+(gamma4-1)/2*M4^2)^(gamma4/(gamma4-1)));
end