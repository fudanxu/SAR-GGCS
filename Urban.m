clc
clear;
close all;
load('.\Data\urban.mat')
hh = urban;
c1 = -40; c2 = 10;
c11 = -5; c22 = 3;
N_mean = 100;
% data processing
hh = double(hh);
sz = size(hh);
hh = hh(1:sz(1)-mod(sz(1),2),1:sz(2)-mod(sz(2),2));
[row,col] = size(hh);
num = row*col;
zr = real(hh);
zi = imag(hh);
z = hh.*conj(hh);
za = sqrt(z);
%*****************parameter estimate*****************************************************
[alpha,beta] = Para_K(hh);
fprintf('The eatimated parameter is alpha = %.4f, beta = %.4f\n',alpha,beta);
A = linspace(min(sqrt(hh(:).*conj(hh(:)))), max(sqrt(hh(:).*conj(hh(:)))), 500);
n = 1;
P_A = K_dis(alpha,beta,n,A);
counts = hist(za(:),A);
P_A_real = counts./numel(za)./(A(end)-A(end-1));
%************estimate the correlation of real SAR image***************************
rho_sr = autocorr2d(zr);
rho_si = autocorr2d(zi);
rho_s = (rho_sr+rho_si)./2;
rho_sa = autocorr2d(za);
rho_sI = autocorr2d(za.^2);
rho_c = autocorr2d(hh);
% ******compute the Convolution kernel of real/imag image********************************
rho_z_mean = rho_s;
rho_z = rho_z_mean;
Rho_z = ifftshift(fft2(fftshift(rho_z)));
Rho_z(imag(Rho_z)<1e-10) = real(Rho_z(imag(Rho_z)<1e-10));
Hr = sqrt(abs(Rho_z));
hr = ifftshift(ifft2(fftshift(Hr)));
hr = Select_center_win(hr,20);
% ********************compute the Convolution kernel of Amplitude image ********************************
rho_za = rho_sa;
Rho_A = ifftshift(fft2(fftshift(rho_za)));
Rho_A(imag(Rho_A)<1e-10) = real(Rho_A(imag(Rho_A)<1e-10));
Ha = sqrt(Rho_A);
ha = ifftshift(ifft2(fftshift(Ha)));
ha(imag(ha)<1e-10) = real(ha(imag(ha)<1e-10));
% ********************compute the Convolution kernel of Intensity image ********************************
rho_zI = rho_sI;
Rho_I = ifftshift(fft2(fftshift(rho_zI)));
Rho_I(imag(Rho_I)<1e-10) = real(Rho_I(imag(Rho_I)<1e-10));
HI = sqrt(Rho_I);
hI = ifftshift(ifft2(fftshift(HI)));
hI(imag(hI)<1e-10) = real(hI(imag(hI)<1e-10));
%% **********************single-point distribution*******************************
% K·Ö²¼ pdf
P_A = K_dis(alpha,beta,n,A);
%% *********************Coherent scatterer model*********************************
sigma = alpha./beta;
s = sqrt(sigma./2);
r = alpha;
p = alpha./(alpha+N_mean);
%*******************************Generate scatterer number N******************************************
%***********compute the Convolution kernel of Scatterer number ******************************
% (1)******************************Obtain the rho_N according to rho_I and rho_x*********************************
N0 = nbinrnd(r,p,row,col);
mean_N = mean(N0(:));
sigma_N = sqrt(var(N0(:)));
rho_x = autocorr2d(zr);
rho_I = autocorr2d(za.^2);
rho_A = autocorr2d(za);
% £¨2£©*******************Compute the rho_N according to rho_I and rho_x********
rho_N1 = (rho_I.*(2*sigma_N^2+mean_N^2)-(sigma_N^2+mean_N^2).*(rho_x.^2))./(sigma_N^2);
fprintf('The maximum value of computed rho_N is %.4f\n',max(rho_N1(:))); 
rho_N2 = (rho_I.*(2*sigma_N^2+mean_N^2))./(sigma_N^2);
[r0,c0] = find(rho_x==max(rho_x(:)));
cx = rho_x(:,c0);
kc = find(cx>=1/exp(1)./10);
rx = rho_x(r0,:);
kr = find(rx>=1/exp(1)./10);
rho_N = rho_N2;
rho_N(kc(1):kc(end),kr(1):kr(end)) = rho_N1(kc(1):kc(end),kr(1):kr(end));
%******************************************************************************
Rho_N = ifftshift(fft2(fftshift(rho_N)));
Rho_N(imag(Rho_N)<1e-10) = real(Rho_N(imag(Rho_N)<1e-10));
Hn = sqrt(Rho_N);
hn = ifftshift(ifft2(fftshift(Hn)));
% hn(imag(hn)<1e-10) = real(hn(imag(hn)<1e-10));
hn = abs(hn);
%********Generate the correlated N************
rng(1000000);
G = wgn(row,col,0);
% **************************Frequency multiple*************
H = ifftshift(fft2(fftshift(hn)));
F = ifftshift(fft2(fftshift(G)));
GF = H.*F;
% GF = Hn.*F;
G2 = ifftshift(ifft2(fftshift(GF)));
G2(imag(G2)<1e-10) = real(G2(imag(G2)<1e-10));
%**************************************************************
centers_N = 0:N_mean*3;
N_G = nbinrnd(r,p,row,col);
[counts] = hist(N_G(:),centers_N);
PN_G = counts./numel(N_G)./(centers_N(end)-centers_N(end-1));
FN_G = cdf_pr(PN_G,centers_N);
[N,~] = DisTrans_to_Fz(FN_G,centers_N,G2);
N = round(N);
% generate the gaussian scattering field
rng(1000000);
R_simu = normrnd(0,s./sqrt(N_mean),row,col,max(max(N)));
I_simu = normrnd(0,s./sqrt(N_mean),row,col,max(max(N)));
% ********************************generate the correlated gaussian scattering field********************************
h0 = waitbar(0,'Please wait...');
for k = 1:max(max(N))
     R_simu(:,:,k) = imfilter(R_simu(:,:,k),hr,'replicate'); 
    [m,sigma] = normfit(reshape(R_simu(:,:,k),1,[]));
    R_simu(:,:,k) = (R_simu(:,:,k)-m)./sigma;
    R_simu(:,:,k) = R_simu(:,:,k).*s./sqrt(N_mean);
    I_simu(:,:,k) = imfilter(I_simu(:,:,k),hr,'replicate'); 
    [m,sigma] = normfit(reshape(I_simu(:,:,k),1,[]));
    I_simu(:,:,k) = (I_simu(:,:,k)-m)./sigma;
    I_simu(:,:,k) = I_simu(:,:,k).*s./sqrt(N_mean);
    waitbar(k/max(max(N)),h0);
end
close(h0)
%******obtain the real image and imaginary image ************
R_simu_sum = zeros(row,col);
I_simu_sum = zeros(row,col);
% *********first N summation**********
for k1 = 1:row
    for k2 = 1:col
        R_simu_sum(k1,k2,:) = sum(R_simu(k1,k2,1:N(k1,k2)));
        I_simu_sum(k1,k2,:) = sum(I_simu(k1,k2,1:N(k1,k2)));
    end
end
%******obtain Amplitude image, Intensity image ************
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
theta = atan2(I_simu_sum,R_simu_sum);
II_simu = A_simu.^2;
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
A_corr_w = Select_center_win(A_corr,10);
rho_sa_w = Select_center_win(rho_sa,10);
%% £¨a£©real data
figure; imagesc(20.*log10(za(1:row,1:col))); axis equal tight off; colorbar; colormap('gray'); caxis([c1,c2]);
tightfig;
print('-dtiff','-r300',['Results\','Urban_real'])
%% £¨b£©simulated data
figure; imagesc(20.*log10(A_simu)); axis equal tight off; colorbar; colormap('gray'); caxis([c1,c2]); tightfig;
print('-dtiff','-r300',['Results\','Urban_simulated'])