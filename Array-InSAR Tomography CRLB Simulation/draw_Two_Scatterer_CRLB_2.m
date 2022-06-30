function [s_crlb,h] = draw_Two_Scatterer_CRLB_2(lam,r0,B,SNR,Scale)

% Author: Yexian Ren，renyexian@foxmail.com

%% Rayleigh resolution
N = length(B);

rho_s = lam*r0/2/(max(B)-min(B));

disp(['Rayleigh resolution: ',num2str(rho_s)]);

%% Elevation
Ns = Scale * 1000;

x = linspace(0,Scale,Ns);

% elevation difference between two scatterers (normalized by rou_s)
normalized_delt_s = linspace(0,Scale,Ns);
s1 = zeros(1,Ns);
s2 = normalized_delt_s .* rho_s;

%% CRLB of single scatterer
snr = 10^(SNR/10);
crlb1 = (lam*r0)./(4*pi)./sqrt(2*N*snr)./std(B,1);

%% CRLB of two scatterers
% empirical function
alpha = (s2 - s1)./rho_s;
c0 = sqrt(2.57.*(alpha.^(-1.5)-0.11).^2+0.62);
c0(c0<1)=1;

crlb2 = c0.*crlb1;

% CRLB
s1_crlb_low = s1 - crlb2;
s1_crlb_high = s1 + crlb2;

s2_crlb_low = s2 - crlb2;
s2_crlb_high = s2 + crlb2;

s_crlb = [s1_crlb_low(:),s1_crlb_high(:),s2_crlb_low(:),s2_crlb_high(:)];

%% draw CRLB
h=figure();
h11 = plot(x(:),s1(:)/rho_s,'LineStyle','-','LineWidth',1,'Color',[0,0,0]);
hold on
h21 = plot(x(:),s2(:)/rho_s,'LineStyle','-','LineWidth',1,'Color',[0,0,0]);
hold on
h12 = plot(x(:),s1_crlb_low(:)/rho_s,'LineStyle','--','LineWidth',1,'Color',[0,0,0]);
hold on
h13 = plot(x(:),s1_crlb_high(:)/rho_s,'LineStyle','--','LineWidth',1,'Color',[0,0,0]);
hold on
h22 = plot(x(:),s2_crlb_low(:)/rho_s,'LineStyle','--','LineWidth',1,'Color',[0,0,0]);
hold on
h23 = plot(x(:),s2_crlb_high(:)/rho_s,'LineStyle','--','LineWidth',1,'Color',[0,0,0]);
hold on
%axis([0.0 1.0 -0.2 1.2])
set(h,'position',[20,20,500,700]) %[left bottom width height]
set(gca,'xlim',[0,Scale],'xtick',[0:0.2:Scale]);
set(gca,'ylim',[-0.2,1.2*Scale],'ytick',[-0.2:0.2:1.2*Scale]);
set(gca,'XDir','reverse')
xlabel(['Normalized true elevation distance (','{\delta}_{s}/{\rho}_{s}',')']);
ylabel(['Estimated elevation position(','/{\rho}_{s}',')']);

title(['CRLB',' (M = ',num2str(N),'; SNR = ',num2str(SNR),' [dB];',' {\rho}_{s} = ',num2str(rho_s),' [m])']);
% legend([h11,h21],'Truth')
% legend([h12,h13,h22,h23],['Truth','\pm','CRLB (N=',num2str(N),')'])
legend([h11,h12],'Truth',['Truth','\pm','CRLB'])
%axis equal;
%pbaspect([1 1 1]) 
grid on;box on;

end