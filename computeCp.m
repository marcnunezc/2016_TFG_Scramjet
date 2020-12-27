function cp=computeCp(T,a,i)
cp=a(i,1)+a(i,2)*T+a(i,3)*T.^2+a(i,4)*T.^3/4+a(i,5)*T.^4;
end