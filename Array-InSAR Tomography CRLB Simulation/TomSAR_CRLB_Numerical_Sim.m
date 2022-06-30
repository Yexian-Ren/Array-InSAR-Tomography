function CRLB = TomSAR_CRLB_Numerical_Sim(Ksi,Scatterers,SNR)

% Author: github.com/Yexian-Ren, renyexian@foxmail.com

% Numerical simulation of CRLB
% Main problem: Fisher matrix is likely to be a singular matrix. Use pinv?

% input:
% Ksi = - (2 * Bn)/(lam * r)
% Scatterers [a,fai,s]
% SNR signal-to-noise ratio (dB)

% output:
% CRLB (3*K)*(3*K) martrix

N = length(Ksi);
K = size(Scatterers,1);

A = Scatterers(:,1);
Fai = Scatterers(:,2);
S = Scatterers(:,3);

meanPs = mean(A.^2);
sigma2 = meanPs * 10^(-SNR/10);

J = zeros(3*K,3*K);

for ii = 1:K
    for jj = 1:K
        a_i = A(ii);
        a_j = A(jj);
        fai_i = Fai(ii);
        fai_j = Fai(jj);
        s_i = S(ii);
        s_j = S(jj);
        J_sub = GetFisherInform(Ksi,a_i,a_j,fai_i,fai_j,s_i,s_j,sigma2);
        J(1+(ii-1)*3:ii*3,1+(jj-1)*3:jj*3) = J_sub;
    end
end

CRLB = pinv(J);
%CRLB = inv(J); % singular

end

function J_sub = GetFisherInform(Ksi,a_i,a_j,fai_i,fai_j,s_i,s_j,sigma2)

delt_ij = 2*pi*Ksi(:)*(s_i-s_j)+(fai_i-fai_j);
% delt_ij = 2*pi*Ksi(:)*(s_j-s_i)+(fai_j-fai_i);

cc = sum(cos(delt_ij));

ss = sum(sin(delt_ij));

cc_ksi =  sum(Ksi(:).*cos(delt_ij));

ss_ksi =  sum(Ksi(:).*sin(delt_ij));

cc_ksi2 = sum(Ksi(:).*Ksi(:).*cos(delt_ij));

% J_sub =  (2/sigma2).*[cc, -a_j*ss, -2*pi*a_j*ss_ksi;...
%                       a_i*ss, a_i*a_j*cc, 2*pi*a_i*a_j*cc_ksi;...
%                       2*pi*a_i*ss_ksi, 2*pi*a_i*a_j*cc_ksi, 4*pi*pi*a_i*a_j*cc_ksi2;];

J_sub =  (2/sigma2).*[cc, a_j*ss, 2*pi*a_j*ss_ksi;...
                      -a_i*ss, a_i*a_j*cc, 2*pi*a_i*a_j*cc_ksi;...
                      -2*pi*a_i*ss_ksi, 2*pi*a_i*a_j*cc_ksi, 4*pi*pi*a_i*a_j*cc_ksi2;];

end