function [af,bw,gain ] = radiation_pattern(amp)
%%
% calc radiation pattern of given current distribution for a line array
% antenna
% parameters: 
% amp : current distribution
M = 100000;
af = zeros(1,M);% array factor buffer
f0 = 3e9;
lambda = 3e8/f0;
d = lambda/2;% elements spacing
k =2*pi/lambda;% wave constant
imp = 120*pi;%wave impdance
I0= 1;
r = 100;
theta = linspace(-pi/2,pi/2,M)+pi/10e10;
N = length(amp);
af = zeros(1,M);
a = 0;
b = 0;
for j = 1:1:N
    a= a+amp(j);
    b = b + power(amp(j),2);
    af = af+ exp(1i*(j-1)*k*d*sin(theta))*amp(j);
end
D = 10*log10(a^2/b);
tau = a^2/b/N;
af = abs(af );
af = 20*log10(af/max(af));

error = 0.0001;
% af=U;
pos1_3dB = [];
pos_max = find(max(af)==af);
tmp = abs(af+3);
% while(isempty(pos1_3dB))
%     pos1_3dB = find(abs(((af(1:pos_max)-af(pos_max)))+3) < error);
%     error = error + 0.0001;
% end
pos1_3dB=find(tmp(1:pos_max)==min(tmp(1:pos_max)));
error = 0.0001;
pos2_3dB = [];
% while(isempty(pos2_3dB))
%     pos2_3dB = find(abs(((af(pos_max:end)-af(pos_max)))+3) < error);
%     error = error + 0.0001;
% end
% HP = 2*asind(sigma/L/pi*sqrt((acosh(10^(R/20)))^2-(acosh(10^(R/20)/sqrt(2)))^2));
pos2_3dB=find(tmp(pos_max:end)==min(tmp(pos_max:end)));
BeamWidth= (theta(pos2_3dB(1)+pos_max)-theta(pos1_3dB(end)))/pi*180;
% figure;
% plot(theta/pi*180,af);
% grid on;
% hold on;
%str = strcat('L=', num2str(L),'\lambda, SLL = -',num2str(R) ,'dB ');
bw = strcat(' beamwidth=',num2str(BeamWidth),'\circ');
gain = strcat('directivity =  ',num2str(D),'dB');

%text(pi*2/3,-5,str,'fontsize',12);
% text(pi*2/3,-8,bw,'fontsize',12);
% text(pi*2/3,-10,gain,'fontsize',12);
% title('Radiation Pattern');
% xlabel('Phase');
% ylabel('Amplitude');
% ylim([-60 0]);

end