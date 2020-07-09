N = 50;
M = 50;
C = 0.030;
Z1 = 30;
X1 = 10;
delxbar = 1/N;
delzbar = 1/M;
const1 = X1*X1/(Z1*Z1);
ITER = 10000;
for I = 1: N+1
    for J = 1: M+1
        p(I,J) = 0.0;
    end
end
sum(1) = 0.0;
for K =1: ITER
    sumij = 0.0;
    for I = 2:N
        X(I) = 1/N*(I-1);
        h = 2/3*(2-X(I));
        hm = 2/3*(2-X(I)-0.5*delxbar);
        hp = 2/3*(2-X(I)+0.5*delxbar);
        hm1 = 2/3*(2-X(I)- delxbar);
        hp1 = 2/3*(2-X(I)+delxbar);
        cubh = h*h*h;
        cubhm = hm*hm*hm;
        cubhp = hp*hp*hp;
        const2 = (cubhp+cubhm+2*const1*cubh);
        A = const1*cubh/const2;
        C = cubhp/const2;
        D = cubhm/const2;
        E = -0.5*delxbar*(hp1-hm1)*C/(const2*X1);
        for J=2:M
            Z(J) = 1/M*(J-1);
            p(I,J) = A* p(I,J+1) + A* p(I,J-1) + C* p(I+1,J) + D* p(I-1,J) - E;
            sumij = sumij +  p(I,J);
        end
    end
    sum(K+1) = sumij;
    percentage = abs(sum(K+1)- sum(K))/abs(sum(K+1));
    if percentage < 0.0001
        break
    end
end
y = K
surf(p);
title('Pressure profile of journal bearing');
xlabel('along the z-direction');
ylabel('along the x-direction');
zlabel('pressure variation');
        
        
        
        
        
        
        
      