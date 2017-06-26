%%
% siw array calc
clc;
clear all;
close all;


gain_element = 3; %antenna element gain
gain_desired = 27 ; %desired gain

element_num = 2^(ceil((gain_desired-gain_element)/3.0+0.1));% total element num
A= 2^(floor((gain_desired-gain_element)/3.0/2+0.1));%row
B= 2^(ceil((gain_desired-gain_element)/3.0/2+0.1))% 