clc 
clear all
precision=0.0001;
N=30;
T5=zeros(1,N);
g0=9.81;
p0=1197/101325;
T0=241.3600;
for M4=6:9
M0=linspace(22,22,1);

A5=1;

hold on
for alpha=2
alpha=alpha/2;
f=0.0293;
M4=linspace(M4,M4,N);
p5=linspace(1,1,N);
p4old=linspace(1,1,N);

T4=linspace(500,2500,N);
for i=1
u0(i)=M0(i)*sqrt(1.4301*287*T0);
[T4(i)]=solveInlet(M4(i),T0,M0(i),alpha);
[p4old(i)]=solveInletPerfect(M4(i),T0,T4(i),M0(i),p0,alpha);

[p4(i), T5(i), xi5(i,:), M5(i), u5(i) m5(i),p5(i)]=iterateExit(p4old(i),T4(i),M4(i),p5(i),A5,alpha,precision);

[M9_1(i), T9_1(i), u9_1(i)]=solveNozzle1(T5(i),xi5(i,:),p0,p5(i), M5(i));
[M9_2(i), T9_2(i), u9_2(i), A9_2(i) m9_2(i)]=solveNozzle2(T5(i),xi5(i,:),p0,p5(i), M5(i), T9_1(i),A5);
[M9_3(i), T9_3(i), u9_3(i), A9_3(i)]=solveNozzle3(T5(i),xi5(i,:),p0,p5(i), M5(i), T9_1(i),A5,alpha,precision);

 F1m(i)=(1+f)*u9_1(i)-u0(i);
 F2m(i)=(1+f)*u9_2(i)-u0(i);
 F3m(i)=(1+f)*u9_3(i)-u0(i);
 Isp1(i)=F1m(i)/f/g0;
 Isp2(i)=F2m(i)/f/g0;
 Isp3(i)=F3m(i)/f/g0;

 T1(i)=Isp1(i)*m5(i)*f*g0;
 T2(i)=Isp2(i)*m5(i)*f*g0;
 T3(i)=Isp3(i)*m5(i)*f*g0;
end
% plot(M0,M9_1,'b')
% plot(M0,M9_2,'r')
% plot(M0,M9_3,'k')
% plot(T4,M5,'g')
% for i=1:N
%     div(i)=M9_2(i)/M9_3(i);
% end
% plot(M0,div)
% figure
% hold on

 plot(M0,Isp1)
% plot(M0,Isp2,'r')
 %plot(M0,Isp3,'k')


end
matrix(M4-5,:)=Isp1
end

