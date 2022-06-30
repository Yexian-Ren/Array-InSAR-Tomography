function [s_crlb,h] = draw_Two_Scatterer_CRLB(lam,r0,B,SNR,Scale)

% Author: Yexian Ren，renyexian@foxmail.com
% Note: there is a problem of singularity of Fisher matrix in numerical simulation

N = length(B);
rho_s = lam*r0/2/(max(B)-min(B));
disp(['Rayleigh resolution: ',num2str(rho_s)]);

Ksi = -2*B(:)/lam/r0;

%% Simulate two scatterers
% Amplitude (set to a similar amplitude)
a1 = 1;
a2 = 1;
% Phase
fai1 = 0;
fai2 = pi/2;

%% Slanted range
% sampling rate
Ns = 1000;

x = linspace(0,Scale,Ns);

% elevation difference between two scatterers (normalized by rou_s)
normalized_delt_s = linspace(0,Scale,Ns);
s1 = zeros(1,Ns);
s2 = normalized_delt_s .* rho_s;
%s2 = s1 + normalized_delt_s;

crlb_1 = zeros(1,Ns);
crlb_2 = zeros(1,Ns);

%% Calculate CRLB of two scatterers

flag = waitbar(0,'Calculation...');

Ms = 1000; % 1000 phases are simulated randomly for smoothing CRLB curves

for ii = 1:1:Ns  
    waitbar(ii/Ns)
    
   % Simulate scatterers with the same intensity and random phase, 
   % and obtain the averaged CRLB boundary curve for Ms times
   
   CRLB = 0;
   CRLB_s1 = zeros(1,Ms);
   CRLB_s2 = zeros(1,Ms);
    for jj = 1:Ms
        a1 = 1;
        a2 = 1;
        fai1 = -pi+2*pi*rand(1,1);
        fai2 = -pi+2*pi*rand(1,1);
        Scatterers = [a1,fai1,s1(ii);a2,fai2,s2(ii)];
        CRLB = TomSAR_CRLB_Numerical_Sim(Ksi,Scatterers,SNR);
        CRLB_s1(1,jj) = CRLB(3,3);
        CRLB_s2(1,jj) = CRLB(6,6);
    end
    
    % take average
    crlb_1(ii) = sqrt(mean(CRLB_s1)); 
    crlb_2(ii) = sqrt(mean(CRLB_s2));

end
close(flag)

s1_crlb_low = s1 - crlb_1;
s1_crlb_high = s1 + crlb_1;

s2_crlb_low = s2 - crlb_2;
s2_crlb_high = s2 + crlb_2;

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

set(h,'position',[20,20,500,700]) 
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


