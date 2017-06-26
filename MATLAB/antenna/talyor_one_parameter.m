
clc;
clear all;
close all;
format long;
%
% N= 1000;
% B = linspace(0.1,1.8,N)';
% SLL = 20*log10(4.603*sinh(pi*B)./(pi*B));
%
% plot(B,SLL)
% xlabel('B');
% ylabel('-SLL:dB');
SLL0= -30;
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
B =round((B1+B2)/2*10000)/10000
SLL= SLL_CALC(B)

lambda = 1;
d = lambda/2;
N = 19;
L = d*(N-1);
i = zeros(1,N);

for m = 1:1:N
    i(m) = besselj(0,1i*pi*B*sqrt(1-(2*(m-(N+1)/2)*d/L)^2))
end
i = i/max(i);
plot(i)

x = i;
onoff = 1;
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
    str = strcat('N=', num2str(N),', SLL = ',num2str(SLL) ,'dB by Taylor Synthesis');
    bw = strcat('main lobe beamwidth=',num2str(BeamWidth),'degree, gain =  ',num2str(D),'dB')
    title(str)
    text(60,-10,bw,'fontsize',10);
    ylim([-60 0])
end
radiation_pattern(i);
% dolph_chebyshev(N,SLL0,1)