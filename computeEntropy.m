function s=computeEntropy(T,a,i)
s=a(i,1)*log(T)+a(i,2)*T+a(i,3)*T^2/2+a(i,4)*T^3/3+a(i,5)*T^4/4+a(i,7);
end