%%
% paraboloidal reflector antenna design
% Leon ,yangli0534@yahoo.com
clc;
close all;
clear all;
%% 
% design parameters
G = 46; %dB
gain = 10^(G/10)%gain decimal

fc = 5400e6;%center frequency in Hz
BW = 100e6; % bandwidth in Hz

SLL = -25; % sidelobe level must be <= -25dB
HalfBeamWidth = 1.2;

K = 80;
ea = 0.5;
%%
%parabolodial reflector design
lc = 3e8/fc*1000 % wavelength in mm
ll= 3e8/(fc-BW/2)*1000;

d = K*ll/HalfBeamWidth;% diameter of aperture
d1 = ll/pi*sqrt(gain/ea);
if d1>d
    d = d1;
end
focus= 0.4*d %focus
theta = 2*atand(d/4/focus)

%%
%feed antenna design, paramidal antenna
SA = 40*log10(cosd(theta/2))% attenuation in max theta
EI = -20;%±ﬂ‘µ’’…‰
FT = EI -SA%  ‘⁄max theta∑ΩœÚ¿°‘¥À•ºı
theta1 = theta *sqrt(10/abs(FT))
ah = lc*79/(2*theta1-31) % flare width
bh = lc*88/2/theta1 % flare height
lg = lc/sqrt(1-power(lc/2/ah,2))%waveguide wavelength

Re = 2*bh^2/lg

% standard waveguide BJ-48
a = 47.55;
b = 22.15;
R = Re*(1-b/bh);
Rh= R/(1-a/ah);
ata_H = pi/4*ah^2/lg/Rh/pi

%%efficiency
n= 2;
a1 = 0.1;
ei = 10*log10((n*a1+1)^2*(2*n+1)/(n+1)/(2*n^2*a1^2+2*n*a1+n+1))
delta =1.5;
ed = -686*(delta/lc)^2
eb= 10*log10((1-(n+1)/(n*a1+1)*a*b/pi/d^2/4)^2)
ea= 0.5;
G = 10*log10((pi*d/lc)^2*ea)




