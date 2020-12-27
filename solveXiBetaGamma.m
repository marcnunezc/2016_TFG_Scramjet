

        function xi5=solveXi(p5,Kp,precision,alpha)
j=0;
 sigma0=8.0215;
 beta0=0.6;
 gamma0=0.3;
 x0(1)=beta0;
 x0(2)=gamma0;
 x0(3)=sigma0;
 X=[0;0;0;0];
 B=[0;0;0;0];
 A=[-1 -0.5 0 0;
     -1 -1 0 0;
     -1 0 2 0;
     0 -1 0 2];
 error=1;
   while error>precision
   %while j<7
            
            B(1)=Kp(1)-log(x0(1))-0.5*(log(x0(3))-log(p5));
            B(2)=Kp(2)+log(x0(2));
            B(3)=Kp(3)-log(p5)+log(x0(3));
            B(4)=Kp(4)-log(p5)+log(x0(3));
            X=A\B;
            delta=exp(X(1));
            epsilon=exp(X(2));
            zeta=exp(X(3));
            eta=exp(X(4));
            %--------------------------------------------------------------
            beta=2-2*alpha-2*delta+2*epsilon-zeta+eta;
            gamma=2*alpha-beta-2*epsilon-eta;
            sigma=beta+delta+gamma+epsilon+zeta+eta+alpha*79/21;
            %--------------------------------------------------------------
            xi(1)=beta;
            xi(2)=gamma;
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
   
  xi5=  [beta  gamma delta  epsilon    zeta   eta alpha*79/21]
        end