%%
%horn antenna calc
% Leon ,yangli0534@yahoo.com
clc;
close all;
clear all;

%% freq
fc = 92e9;
lambda = 3e8/fc*1000;

% standard waveguide WR10
a = 2.54;%mm
b = 1.27;%mm

%design parameters
G = 25;
gain = 10^(G/10);

%ah= 0.45*lambda*sqrt(gain);
en = 0.51;

syms ah
eqn = power(ah,4)-a*power(ah,3)+3*gain*power(lambda,2)*b/8/pi/en*ah == 3*power(gain,2)*power(lambda,4)/32/power(pi,2)/power(en,2);
solx = solve(eqn,ah);
syms ah clear;
ah = vpa(solx(2))%aperture  width 
bh = gain*power(lambda,2)/2.04/pi/ah%aperture height 
Rh=  power(ah,2)/3/lambda
R = Rh*(1-a/ah) %distance between  waveguide to aperture
Re = R/(1-b/bh)

HPh = 78*lambda/ah % 3dB beamwidth in H plane
HPe = 54*lambda/bh % 3dB beamwidth in E plane

 