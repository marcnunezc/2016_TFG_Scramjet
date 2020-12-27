

 function xi5=solveXi(p5,Kp,precision,alpha)
 MW=[2.01588 31.99880 28.01340 17.00734 18.01528 1.00794 15.99940];
 j=0;
 sigma0=5;
 beta0=0.1;
 epsilon0=0.1;
 x0(1)=beta0;
 x0(2)=epsilon0;
 x0(3)=sigma0;
 X=[0;0;0;0];
 B=[0;0;0;0];
 B1=[0;0;0];
 A=[0 -1 0 0;
     2 -1 0 0;
     0 -1 2 0;
     0 0 0 2];
 A1=[1 2 0;
     2 0 0;
     -1 -1 1];
 error=1;
   while error>precision

            
            B(1)=Kp(1)-log(x0(1))+0.5*log(x0(2))-0.5*(log(x0(3))-log(p5));
            B(2)=Kp(2)+log(x0(2));
            B(3)=Kp(3)-log(p5)+log(x0(3));
            B(4)=Kp(4)+log(x0(2))-log(p5)+log(x0(3));
            X=A\B;
            gamma=exp(X(1));
            delta=exp(X(2));
            zeta=exp(X(3));
            eta=exp(X(4));
            B1(1)=2*alpha-gamma-eta;
            B1(2)=2-2*delta-gamma-zeta;
            B1(3)=delta+gamma+eta+zeta+alpha*79/21;
            X1=A1\B1;
            beta=X1(1);
            epsilon=X1(2);
            sigma=X1(3);
            %--------------------------------------------------------------
            xi(1)=beta;
            xi(2)=epsilon;
            xi(3)=sigma;
            error=0;
                for i=1:3
                    diff=abs(xi(i)-x0(i));
                    if error<diff
                        error=diff;
                    end
                    x0(i)=xi(i);     
                end              

   end
   
  xi5=1/sigma* [beta  gamma delta  epsilon    zeta   eta alpha*79/21];
        end