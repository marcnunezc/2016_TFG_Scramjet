function gamma=computeGamma5(T,xi,a)
Ru=8.314510;
Cptot=0;
Cp(1)= Ru*computeCp(T, a, 4); %H20
Cp(2)= Ru*computeCp(T, a, 5); %OH
Cp(3)= Ru*computeCp(T, a, 1); %H2
Cp(4)= Ru*computeCp(T, a, 2); %O2
Cp(5)= Ru*computeCp(T, a, 6); %H
Cp(6)= Ru*computeCp(T, a, 7); %O
Cp(7)= Ru*computeCp(T, a, 3); %N2
for i=1:length(xi)
    Cptot=Cptot+Cp(i)*xi(i);
end
gamma=Cptot/(Cptot-Ru);
end