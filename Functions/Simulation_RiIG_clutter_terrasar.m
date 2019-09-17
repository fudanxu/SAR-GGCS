function [R_simu_sum,I_simu_sum]  = Simulation_RiIG_clutter_terrasar(hh,c1,c2,c11,c22,N_mean)
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
%% ***************************parameter estimate*****************************************************
[alpha,beta,delta] = Para_RiIG(hh);
fprintf('The eatimated parameter is alpha = %.4f, beta = %.4f, delta = %.4f\n',alpha,beta,delta);
A = linspace(min(sqrt(hh(:).*conj(hh(:)))), max(sqrt(hh(:).*conj(hh(:)))), 500);
betax = beta;
n = 1;
counts = hist(za(:),A);
P_A_real = counts./numel(za)./(A(end)-A(end-1));
%% ****************************estimate the correlation of real SAR image***************************
rho_sr = autocorr2d(zr);
rho_si = autocorr2d(zi);
rho_s = (rho_sr+rho_si)./2;
rho_sa = autocorr2d(za);
rho_sI = autocorr2d(za.^2);
% ********************compute the Convolution kernel of real/imag image********************************
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
% RiIG分布 pdf
P_A = RiIG_dis(alpha,beta,delta,A);
%% *********************Coherent scatterer model*********************************
gamma = sqrt(alpha.^2-beta.^2);
mu = delta/gamma;
lambda = (N_mean.*gamma).^2;
pd = makedist('InverseGaussian','mu',N_mean,'lambda',lambda);
Lam = random(pd,row,col);
sigma = mu*2;
s = sqrt(sigma./2);
%% *********************************Generate scatterer number N******************************************
% 1*********************compute the Convolution kernel of Scatterer number ******************************
% (1)******************************Obtain the rho_N according to rho_I and rho_x*********************************
N0 = poissrnd(Lam,row,col);
mean_N = mean(N0(:));
sigma_N = sqrt(var(N0(:)));

rho_x = autocorr2d(zr);
rho_I = autocorr2d(za.^2);
rho_A = autocorr2d(za);
% （2）*******************Compute the rho_N according to rho_I and rho_x********
rho_N1 = (rho_I.*(2*sigma_N^2+mean_N^2)-(sigma_N^2+mean_N^2).*(rho_x.^2))./(sigma_N^2);
fprintf('The maximum value of computed rho_N is %.4f\n',max(rho_N1(:))); 
%**********find  value >1 in rho_N ***********************
temp1 = sigma_N^2/(sigma_N^2+mean_N^2);
temp2 = (rho_I-rho_x.^2)./(1-rho_I);
tempk = find(rho_I==max(rho_I(:)));
temp2(tempk) = 0;
fprintf('The value of temp1 is %.4f\nThe maximum value of temp2 is %.4f\n',temp1,max(temp2(:))); 
k = find(temp1<temp2);
if ~isempty(k)
    fprintf('Something wrong: Too small value of sigma_N or wrong data\n ');
    %**********Correct the value >1 in rho_N *****************
    rho_N1(k) = rho_N1(k)./(max(rho_N1(:))+0.1);
    rho_N1(tempk) = 1;
    %*******************************************************
end

rho_N2 = (rho_I.*(2*sigma_N^2+mean_N^2))./(sigma_N^2);
[r0,c0] = find(rho_x==max(rho_x(:)));
cx = rho_x(:,c0);
kc = find(abs(cx)>=1/exp(1)./10);
rx = rho_x(r0,:);
kr = find(abs(rx)>=1/exp(1)./10);
rho_N = rho_N2;
rho_N(kc(1):kc(end),kr(1):kr(end)) = rho_N1(kc(1):kc(end),kr(1):kr(end));
%% 2 *************if select the conv kernel of N_truth***********************************
%% ******************************************************************************
Rho_N = ifftshift(fft2(fftshift(rho_N)));
Rho_N(imag(Rho_N)<1e-10) = real(Rho_N(imag(Rho_N)<1e-10));
Hn = sqrt(Rho_N);
hn = ifftshift(ifft2(fftshift(Hn)));
% hn(imag(hn)<1e-10) = real(hn(imag(hn)<1e-10));
hn = abs(hn);
%% 3 **********Generate the correlated N************
rng(1000000);
G = wgn(row,col,0);
H = ifftshift(fft2(fftshift(hn)));
F = ifftshift(fft2(fftshift(G)));
GF = H.*F;
% GF = Hn.*F;
G2 = ifftshift(ifft2(fftshift(GF)));
G2(imag(G2)<1e-10) = real(G2(imag(G2)<1e-10));
%**************************************************************
centers_N = 0:N_mean*3;
N_G = poissrnd(Lam,row,col);
[counts] = hist(N_G(:),centers_N);
PN_G = counts./numel(N_G)./(centers_N(end)-centers_N(end-1));
FN_G = cdf_pr(PN_G,centers_N);
[N,~] = DisTrans_to_Fz(FN_G,centers_N,G2);
N = round(N);
% ***************Compare the results********************************
% generate the gaussian scattering field
rng(1000000);
betax = beta;
betay = sqrt(beta^2-betax^2);
gamma = sqrt(alpha.^2-beta.^2);
mu = delta/gamma;
sigma = mu*2;
sigma_ai = sigma./N_mean;
bi = sqrt(sigma_ai./2);
% *******************分块计算*********
len = 500;
clip = floor(max(max(N))/len);
Nk = zeros(row,col,clip+1);
for clipk = 1:clip+1
    temp = N-(clipk-1)*len;
    temp(temp<0) = 0;
    temp(temp>len) = len;
    Nk(:,:,clipk) = temp;
end
temp = sum(Nk,3);
figure;
subplot(131); imagesc(temp); axis equal tight 
subplot(132); imagesc(N); axis equal tight 
subplot(133); imagesc(temp-N); axis equal tight 
R_simu_sumk = zeros(row,col,clip+1);
I_simu_sumk = zeros(row,col,clip+1);
for clipk = 1:clip+1
    fprintf('Now is processing the %.1f of total %.1f\n',clipk,clip+1); 
    [R_simu_sumk(:,:,clipk),I_simu_sumk(:,:,clipk)] = Clip_RIN(row,col,betax,betay,mu,s,hr,N_mean,Nk(:,:,clipk));
end
R_simu_sum = sum(R_simu_sumk,3);
I_simu_sum = sum(I_simu_sumk,3);
%% **********obtain Amplitude image, Intensity image ************
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
theta = atan2(I_simu_sum,R_simu_sum);
II_simu = A_simu.^2;