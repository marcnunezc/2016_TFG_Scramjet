clc 
clear all
precision=0.0001;
N=30;
T5=zeros(1,N);
f=0.0293;
g0=9.81;
p0=1197/101325;
T0=273+15-46.64;

MAX=0;
    
 M0=linspace(6,15,N);
% 
A5=1;
% 

 for alpha=2
     if alpha==1
         alpha=1.2;
     end
 alpha=alpha/2
 M4=linspace(4,4,N);
 for i=1:N 
     
    u0(i)=M0(i)*sqrt(1.4301*287*T0);
   T4(i)=solveInlet(M4(i),T0,M0(i),alpha);
%   [T4(i) p4(i)]=solveInletPerfect(M4(i),T0,M0(i),p0,alpha)
 end
% p5=p4; %atm
p5=linspace(1,1,N);
% 
%T4=linspace(500,2500,N);
 for i=1:N
  i
[p4(i), T5(i), xi5(i,:), M5(i), u5(i) m5(i)]=solveExit(T4(i),M4(i),p5(i),A5,alpha,precision);
[M9_1(i), T9_1(i), u9_1(i)]=solveNozzle1(T5(i),xi5(i,:),p0,p5(i), M5(i));
[M9_2(i), T9_2(i), u9_2(i), A9_2(i) m9_2(i)]=solveNozzle2(T5(i),xi5(i,:),p0,p5(i), M5(i), T9_1(i),A5);
[M9_3(i), T9_3(i), u9_3(i), A9_3(i)]=solveNozzle3(T5(i),xi5(i,:),p0,p5(i), M5(i), T9_1(i),A5,alpha,precision);
% 
 F1m(i)=(1+f)*u9_1(i)-u0(i);
 F2m(i)=(1+f)*u9_2(i)-u0(i);
 F3m(i)=(1+f)*u9_3(i)-u0(i);
 Isp1(i)=F1m(i)/f/g0;
%Isp2(i)=F2m(i)/f/g0;
% Isp3(i)=F3m(i)/f/g0;
%
% T1(i)=Isp1(i)*m5(i)*f*g0;
% T2(i)=Isp2(i)*m5(i)*f*g0;
% T3(i)=Isp3(i)*m5(i)*f*g0;
if Isp1(i)>MAX
    Mmax=M0(i);
    MAX=Isp1(i);
end
 end
%  hold on
%  grid on
%  plot(M0,M9_1,'b')
%  plot(M0,M9_2,'r')
%  plot(M0,M9_3,'k')
%  hold off
% plot(T4,M5,'g')

% figure
% hold on

  plot(M0,Isp1)
%  plot(M0,Isp2,'r')
%  plot(M0,Isp3,'k')
% for i=1:N
%     div(i)=M5(i)/M4(i);
% end
%  plot(M4,div,'b')
end
% ylim([6 9])
% xlim([5 15])