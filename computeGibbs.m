function G=computeGibbs(T,a,i,xi5,p5)
G=computeEnthalpy(T,a,i)-computeEntropy(T,a,i)-log(p5)*xi5;
end