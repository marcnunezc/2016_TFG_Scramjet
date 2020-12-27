function h=computeEnthalpy(T,a,i)
h=a(i,1)+a(i,2)*T/2+a(i,3)*T^2/3+a(i,4)*T^3/4+a(i,5)*T^4/5+a(i,6)/T;
end