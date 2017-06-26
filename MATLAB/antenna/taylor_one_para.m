function I = taylor_one_para(N,SLL0)

%SLL0= -SLL;
B1= 0.1;
B2 = 2;
SLL1= SLL_CALC(B1);
SLL2= SLL_CALC(B2);
while(abs(SLL1+SLL0) > 0.0001 || abs(SLL2+SLL0)> 0.0001)
    B =(B1+B2)/2;
    SLL = SLL_CALC(B);
    if SLL + SLL0 > 0
        B2 =B;
    else
        B1=B;
    end
    SLL1= SLL_CALC(B1);
    SLL2= SLL_CALC(B2);
end
B =round((B1+B2)/2*10000)/10000;
SLL= SLL_CALC(B);

lambda = 1;
d = lambda/2;
%N = 19;
L = d*(N-1);
i = zeros(1,N);

for m = 1:1:N
    i(m) = besselj(0,1i*pi*B*sqrt(1-(2*(m-(N+1)/2)*d/L)^2));
end
I = i/max(i);