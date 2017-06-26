%% waveguide slot equvalient circuit
% Leon yangli0534@yahoo.com
clc;
clear all;
close all;
format long;
%% constants
N = 18; %slot numbers
SLL = -35; % sidelobe level
u0 = 4*pi*1e-7;% permeability
e0 = 8.854187817e-12;% permittivity in free space
c = 1/sqrt(u0*e0);% light velocity in free space
%% parameters
f0 = 35e9;% frequency
l = c/f0*1e3; % lambda in free space
% a = 247.65; % width length in mm
% b = a/2 ;% height length in mm
%WR-284
% a = 20;
% b = 10;
 a = 5.69;
 b = 2.845;
%a = 7.112 ; %waveguide width
%b = 3.556 ; % waveguide height

t = 1 ;% thickness in mm
w = a*0.0625/0.9 ;% width of the slot in mm l/200 < w < l/10
lc = 2*a ; % cutoff lambda
lg = l/sqrt(1-(l/lc)^2); % waveguide wavelength

G_2_slot=1.0/N;
New_G1=2.09*(lg/l)*(a/b)*(cos(0.464*pi*l/lg)-cos(0.464*pi))^2;
New_Y=G_2_slot/New_G1;
Soff=(a/pi)*sqrt(abs(asin(New_Y)));% offset distance

Slot_wl=0.210324*G_2_slot^4-0.338065*G_2_slot^3+0.12712*G_2_slot^2+0.034433*G_2_slot+0.48253;
l_slot=l*Slot_wl;%slot length
w_slot=a*0.0625/0.9;%slot width
s_slot=lg/2 ; % space betwenn slots


d = linspace(0,a/2,1000)+1e-10;
%% stevenson 利用等效传输线理论和格林函数，推导了归一化电导
g = 2.09*a*lg/b/l*power(cos(pi/2*l/lg),2)*power(sin(pi/a*d),2);

%% H.Y.Yee计算方法解决了谐振长度随偏移距离变化的问题
beta = 2*pi/lg; % TE10模传播常数 
Y0 = beta/f0/2/pi/u0;%波导特征导纳
%考虑缝隙影响，把缝隙看成一段分支波导



% figure(1);
% plot(d,g,'k');
% grid on;
% grid minor;
% title('宽边纵缝的归一化谐振电导与偏置量的变化 ');
% xlabel('偏置距离:mm');
% ylabel('归一化电导')
% 



[BeamWidth, D, I]= dolph_chebyshev(N,SLL,0);
figure;
subplot(221);
plot(I,'-o');

title('Dolph-Chebyshev 综合');
hold on;
[af,bw,gain] = radiation_pattern(I);
len = length(af);
angle = linspace(-90,90,len);
subplot(222);
plot(angle,af);
ylim([-60 0]);
grid on;
title('方向图')
text(20,5,bw);
g = I.^2/sum(I.^2);
g = I/sum(I);
% 
%figure;
subplot(223)
plot(g,'k*');
title('等效电导')
hold on;
% % g = a/b*2.09*lg/l*power(cos(pi*l/lg/2),2)*power(sin(pi*x/a),2)
 x = a/pi*asin(sqrt(g/(2.09*a/b*lg/l*power(cos(pi/2*l/lg),2))));
%figure;
subplot(224)
plot(x,'r*');
title('缝隙偏移')

%%
%Taylor one parameter
I= taylor_one_para(N,SLL);
figure;
subplot(221);
plot(I,'-o');

title('taylor单变量综合');
hold on;
[af,bw,gain] =  radiation_pattern(I);
len = length(af);
angle = linspace(-90,90,len);
subplot(222);
plot(angle,af);
ylim([-60 0]);
grid on;
title('方向图')
text(20,5,bw);
g = I.^2/sum(I.^2);
g = I/sum(I);
% 
%figure;
subplot(223)
plot(g,'k*');
title('等效电导')
hold on;
% % g = a/b*2.09*lg/l*power(cos(pi*l/lg/2),2)*power(sin(pi*x/a),2)
 x = a/pi*asin(sqrt(g/(2.09*a/b*lg/l*power(cos(pi/2*l/lg),2))));
%figure;
subplot(224)
plot(x,'r*');
title('缝隙偏移')

%%
%taylor line source
I= taylor_line(N,SLL);
figure;
subplot(221);
plot(I,'-o');

title('taylor线源综合');
hold on;
[af,bw,gain] =  radiation_pattern(I);
len = length(af);
angle = linspace(-90,90,len);
subplot(222);
plot(angle,af);
ylim([-60 0]);
grid on;
title('方向图')
text(20,5,bw);
g = I.^2/sum(I.^2);
g = I/sum(I);
% 
%figure;
subplot(223)
plot(g,'k*');
title('等效电导')
hold on;
% % g = a/b*2.09*lg/l*power(cos(pi*l/lg/2),2)*power(sin(pi*x/a),2)
 x = a/pi*asin(sqrt(g/(2.09*a/b*lg/l*power(cos(pi/2*l/lg),2))));
%figure;
subplot(224)
plot(x,'r*');
title('缝隙偏移')