function [p4, T5, xi5, M5, u5, m5, p5]=iterateExit(p4old,T4,M4,p5,A5,alpha,precision)
diff=1;
    while diff>precision
    [p4, T5, xi5, M5, u5, m5]=solveExit(T4,M4,p5,A5,alpha,precision);
    diff=abs(p4old-p4);
        if p4old>p4
            p5=p5+p5*diff/p4;
        else
            p5=p5-p5*diff/p4;
        end
    end
end
