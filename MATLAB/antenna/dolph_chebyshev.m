%%
%dolph_chebyshev array

%%
% Author : Leon
% Email: yangli0534@gmail.com
% blog: www.cnblogs.com/hiramlee0534
% dolph-chebyshev ????????
% ???????????? ?????????? ??????????????????????

% clc;
% clear all;
% close all;
%format long;
function [BeamWidth, D,x] = dolph_chebyshev(N,SLL,onoff)
%%
% desired parameters
% N = 30; %elemetns number
% SLL = -40;%Sidelobe level

%%
% ????????????
R0 = 10^(-SLL/20);
x = zeros(1,N-1);
u = zeros(1,N-1);
s = zeros(1,N);
I = zeros(1,N);
x0 = 1/2*(power(R0+sqrt(R0^2-1),1/(N-1))+power(R0-sqrt(R0^2-1),1/(N-1)));
% for i = 1:1:round((N-1)/2)
%     x(i) = cos((2*i-1)*pi/(2*(N-1)));
%     x(N-i) =  cos((2*i-1)*pi/(2*(N-1)));
% end
% u = 2*acos(x/x0);
% for i = 1:1:round((N-1)/2)
%     %x(i) = cos((2*i-1)*pi/(2*(N-1)));
%     u(N-i) =  -u(i)
% end
% theta = linspace(0,pi,M);
% S = zeros(1,M);
% S = abs(exp(1i*pi*cos(theta))-exp(1i*u(1))).*abs(exp(1i*pi*cos(theta))-exp(1i*u(2))).* ...
%     abs(exp(1i*pi*cos(theta))-exp(1i*u(3))).*abs(exp(1i*pi*cos(theta))-exp(1i*u(4)));
% % for i = 1:1:N
% %     s = s+exp(j*(i-1)*u(i))
% % end
% plot(theta,10*log10(S/max(S)))
if mod(N,2) == 1
    K = round((N-1)/2);
    
    for n =2:1:K + 1
        for p = n:1:K+1
            I(n) = I(n)+(-1)^(K-p+1)*K*factorial(p+K-2)/factorial(p-n)/factorial(p+n-2)/factorial(K-p+1)*power(x0,2*(p-1))   ;
        end
        
    end
    for p = 1:1:K+1
        I(1) = I(1)+(-1)^(K-p+1)*K*factorial(p+K-2)/factorial(p-1)/factorial(p-1)/factorial(K-p+1)*power(x0,2*(p-1))   ;
    end
    
    
    %x(1) = 1;
    %x(N) = 1;
    for i = 1:1:K+1
        x(i) = I(K+2-i)/I(K+1);
        x(N+1-i) = x(i);
    end
    
else
    K = round(N/2);
    
    for n =1:1:K
        for p = n:1:K
            I(n) = I(n)+(-1)^(K-p)*(2*K-1)*factorial(p+K-2)/2/factorial(p-n)/factorial(p+n-1)/factorial(K-p)*power(x0,2*p-1);
        end
        I(N+1-n) = I(n);
    end
    
    for i = 1:1:K
        x(i) = I(K+1-i)/I(K);
        x(N+1-i) = x(i);
    end
end
x= x/max(x);
if onoff ==1
    figure;
    plot(x);
    title('normalized amplitude distribution ');
end
%%
% AF
M = 2000; % theta sample points
E0 = zeros(1,M);%
H0 = zeros(1,M);%
U0= zeros(1,M);%
f0 = 1e9;
lambda = 3e8/f0;
d = lambda/2;% elements spacing
k =2*pi/lambda;% wave constant
imp = 120*pi;%wave impdance
I0= 1;
r = 100;
l = 0.01;

af = zeros(1,M);% array factor buffer

theta = linspace(0,pi,M)+pi/10e10;
phi = linspace(0,2*pi,M);
af = zeros(1,M);
a = 0;
b = 0;
for j = 1:1:N
    a= a+x(j);
    b = b + power(x(j),2);
    af = af+ exp(1i*(j-1)*k*d*cos(theta))*x(j);
end
af = abs(af/N );
E0 = af *imp;
Wav = real(E0 .* E0 )/2/imp;
Prad = imp*pi*3*(I0*l/lambda)^2;
theta = theta/pi*180;
af= 10*log10(af/max(af));
U= real(E0 .*E0)*r^2/2/imp;
U= 10*log10(U/max(U));

error = 0.01;
af=U;
pos1_3dB = [];
pos_max = find(max(af)==af);
while(isempty(pos1_3dB))
    pos1_3dB = find(abs(((af(1:pos_max)-af(pos_max)))+3) < error);
    error = error + 0.005;
end
error = 0.01;
pos2_3dB = [];
while(isempty(pos2_3dB))
    pos2_3dB = find(abs(((af(pos_max:end)-af(pos_max)))+3) < error);
    error = error + 0.005;
end
BeamWidth= (theta(pos2_3dB(1)+pos_max)-theta(pos1_3dB(end)));
D = 10*log10(a^2/b);
tau = a^2/b/N;
if onoff == 1
    figure;
    plot(theta,U);
    hold on;
    axis([0 180 -120 0]);
    str = strcat('N=', num2str(N),', SLL = ',num2str(SLL) ,'dB by Dolph-Chebyshev Synthesis');
    bw = strcat('main lobe beamwidth=',num2str(BeamWidth),'degree, gain =  ',num2str(D),'dB')
    title(str)
    text(60,-10,bw,'fontsize',10);
end