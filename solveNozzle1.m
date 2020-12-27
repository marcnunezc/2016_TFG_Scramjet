function [M9 T9 u9]=solveNozzle1(T5,xi5,p0,p5,M5)

T9old=0;

MW=[2.01588 31.99880 28.01340 17.00734 18.01528 1.00794 15.99940];
 MW5=xi5(1)*MW(4)+xi5(2)*MW(5)+xi5(3)*MW(1)+xi5(4)*MW(2)+xi5(5)*MW(6)+xi5(6)*MW(7)+xi5(7)*MW(3);
hf=[8468.102 8680.104 8670.104 8.813E3 9.904E3 6.197E3 6.725E3];
Ru=8.314510;
R5=Ru/MW5*1000;

a = importfile1('termodata.dat', [2,4,6,8,10,12,14],[2,4,6,8,10,12,14]);
diff=1;
while diff>0.001
Tm=(T9old+T5)/2;
gamma=computeGamma5(Tm,xi5,a);
pt5=p5*(1+0.5*(gamma-1)*M5^2)^(gamma/(gamma-1));
Tt5=T5*(1+0.5*(gamma-1)*M5^2);
M9=sqrt(2/(gamma-1)*((pt5/p0)^((gamma-1)/gamma)-1));
T9=Tt5/(1+0.5*(gamma-1)*M9^2);
diff=abs(T9-T9old);
T9old=T9;
end 

u9=0.985*M9*sqrt(computeGamma5(T9,xi5,a)*R5*T9);
end