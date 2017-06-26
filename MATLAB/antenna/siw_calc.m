%%
%SIW design 
clc;
clear all;
close all;

%%
% frequency
f0 = 27.5e9;% cnter frequency :GHz
varepsilon_0  = 8.854187817e-12;
mu0 = 4*pi*1e-7;
c = 1/sqrt(varepsilon_0*mu0);
lambda = c/f0*1000 ;% wavelength : mm

%%
% dielectric parameters 
varepsilon_r = 2.2; % Rogers 5880 
h = 0.254; % dielectirc thickness : 0.5 mm 

%%
% equivalent waveguide 
% standard waveguide WR34 : 22GHz-- 33GHz
a = 8.636; %  width in mm
b = 4.318; % height in mm

%%
% calc
lambda_g = lambda/sqrt(varepsilon_r - (lambda/2/a)^2)
dmin = lambda_g /20 %
dmax = lambda_g /18 % 
d = 0.40 % via diameter
p = 2*d % via space
ar = a+(d)^2/0.95/p
