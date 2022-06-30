
% This code was written by Yexian Ren 
% E-mail: renyexian@foxmail.com, QQ: 2538715345

% draw CRLB curve

clear;clc;

%% Airborne Array-InSAR Emei Data CRLB simulation

% The airborne data incidence angle changes greatly!

lam = 0.031;
r0 = 2.543272808657184e+03;
Bv0 = [0;0.164;0.330;0.495;0.641;0.829;0.993;1.158;1.324;1.469;1.658];
M = length(Bv0);
ind = [1:M];

figure
scatter(ind(:),Bv0(:),40,'black','o','linewidth',1);
xlabel(['m']);
ylabel('b_{m}');
title(['Baselines ','(M = ',num2str(11), '; \lambda = ',num2str(0.031),' [m]; r = ',num2str(2543),' [m])']);
grid on
box on

%% Spaceborne TerraSAR Data CRLB simulation
% 
% r0 = 5.5932e5;
% c = 3e8; fre = 9.65e9; lam = c/fre;   %wavelength
% Bv0 = [185.9077,  30.3082,  47.8664,  121.2385,  -13.738,...
%    -105.257,   -115.4251,   -171.9019,  0,     -2.7988, ...
%    96.8326,   70.472,   212.0481];
% M = length(Bv0);
% ind = [1:M];
% 
% figure
% scatter(ind(:),Bv0(:),80,'black','*');
% xlabel(['m']);
% ylabel('b_{m}');
% title(['Baselines ','(M = ',num2str(13), '; \lambda = ',num2str(0.031),' [m]; r = ',num2str(559320),' [m])']);
% grid on
% box on

SNR = 5; % dB
Scale = 1.5;

% Numerical simulation (smoothing strategy)
[s_crlb,h0] = draw_Two_Scatterer_CRLB(lam,r0,Bv0,SNR,Scale);

% Approximation
[s_crlb,h1] = draw_Two_Scatterer_CRLB_2(lam,r0,Bv0,SNR,Scale);

