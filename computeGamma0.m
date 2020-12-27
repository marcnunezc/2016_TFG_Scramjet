function gamma=computeGamma0(T,xi,a)
Ru=8.314510;
Cptot=0;
Cp=[0 0];

Cp(1)= Ru*computeCp(T, a, 2); %O2
Cp(2)= Ru*computeCp(T, a, 3); %N2


for i=1:length(xi)
    Cptot=Cptot+Cp(i)*xi(i);
end
gamma=Cptot/(Cptot-Ru);
end