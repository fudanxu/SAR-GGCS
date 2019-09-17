clc
clear
close all
load('.\Data\sea4.mat')
% *************************
addpath(genpath(pwd));
hh = sea4;
[row,col] = size(hh);
z = hh.*conj(hh);
za = sqrt(z);
zr = real(hh);
zi = imag(hh);
c1 = 20; c2 = 60;
c11 = -150; c22 = 150;
N_mean = 20000;
[R_simu_sum,I_simu_sum]  = Simulation_RiIG_clutter_terrasar(hh,c1,c2,c11,c22,N_mean);
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
theta = atan2(I_simu_sum,R_simu_sum);
II_simu = A_simu.^2;
A = linspace(min(sqrt(hh(:).*conj(hh(:)))), max(sqrt(hh(:).*conj(hh(:)))), 500);
[alpha,beta,delta] = Para_RiIG(hh);
P_A = RiIG_dis(alpha,beta,delta,A);
counts = hist(za(:),A);
P_A_real = counts./numel(za)./(A(end)-A(end-1));
[counts,centers] = hist(atan2(zi(:),zr(:)),100);
P_Theta_real = counts./numel(zi)./(centers(end)-centers(end-1));
counts = hist(A_simu(:),A);
P_A_simu = counts./numel(A_simu)./(A(end)-A(end-1));
[counts] = hist(theta(:),centers);
P_Theta_simu = counts./numel(theta)./(centers(end)-centers(end-1));
A_corr = autocorr2d(A_simu);
rho_sa = autocorr2d(za);
A_corr_w = Select_center_win(A_corr,20);
rho_sa_w = Select_center_win(rho_sa,20);
%% Plot figures
%% £¨a£©real data
figure; imagesc(20.*log10(za(1:row,1:col))); axis equal tight off; colorbar; colormap('gray'); caxis([c1,c2]);
tightfig;
print('-dtiff','-r300',['Results\','Sea4_real'])
%% £¨b£©simulated data
figure; imagesc(20.*log10(A_simu)); axis equal tight off; colorbar; colormap('gray'); caxis([c1,c2]);
tightfig;
print('-dtiff','-r300',['Results\','Sea4_simulated'])

