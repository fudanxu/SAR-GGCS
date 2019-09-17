function [R_simu_sum,I_simu_sum]  = Simulation_correlated_Rayleigh_clutter_farmland(hh,N,row,col)
w = 10;
% data processing
hh = double(hh);
sz = size(hh);
hh = hh(1:sz(1)-mod(sz(1),2),1:sz(2)-mod(sz(2),2));
zr = real(hh);
zi = imag(hh);
z = hh.*conj(hh);
za = sqrt(z);
%  ******************
if sz(1)<=w || sz(2)<=w
    w = min(sz(1),sz(2))/2;
end
%% ***************************parameter estimate*****************************************************
m = mean(mean(sqrt(hh.*conj(hh))));
sigma = 4*m^2/pi;
A = linspace(min(sqrt(hh(:).*conj(hh(:)))), max(sqrt(hh(:).*conj(hh(:)))), 100);
fprintf('The eatimated parameter is sigma = %.4f\n',sigma);
n = 1;
% N_mean = 20;
b = sqrt(sigma./2);
P_A = raylpdf(A,b);
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
hr = Select_center_win(hr,w);

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
%% ÈðÀû·Ö²¼pdf
b = sqrt(sigma./2);
P_A = raylpdf(A,b);
%% *********************Coherent scatterer model*********************************
s = sqrt(sigma./2);
%% *********************************Generate scatterer number N******************************************
%****************************************Generate the truth Scatterer number***************************************
sigmai = s./sqrt(N);
N_truth = za.^2./(2*sigmai^2);
N_truth = ceil(mean(N_truth(:)));

fprintf('The fitted parameters of N_truth is N_truth = %.4f\n',N_truth);
fprintf('The estimated parameters of N is N = %.4f\n',N);

%% generate the gaussian scattering field
rng(1000000);
R_simu = normrnd(0,s./sqrt(N),row,col,N);
I_simu = normrnd(0,s./sqrt(N),row,col,N);
% ********************************generate the correlated gaussian scattering field********************************
h0 = waitbar(0,'Please wait...');
for k = 1:max(max(N))
     R_simu(:,:,k) = imfilter(R_simu(:,:,k),hr,'replicate'); 
    [m,sigma] = normfit(reshape(R_simu(:,:,k),1,[]));
    R_simu(:,:,k) = (R_simu(:,:,k)-m)./sigma;
    R_simu(:,:,k) = R_simu(:,:,k).*s./sqrt(N);
    I_simu(:,:,k) = imfilter(I_simu(:,:,k),hr,'replicate'); 
    [m,sigma] = normfit(reshape(I_simu(:,:,k),1,[]));
    I_simu(:,:,k) = (I_simu(:,:,k)-m)./sigma;
    I_simu(:,:,k) = I_simu(:,:,k).*s./sqrt(N);
    waitbar(k/max(max(N)),h0);
end
close(h0)
%% ********obtain the real image and imaginary image ************
R_simu_sum = zeros(row,col);
I_simu_sum = zeros(row,col);
% *********first N summation**********
for k1 = 1:row
    for k2 = 1:col
        R_simu_sum(k1,k2,:) = sum(R_simu(k1,k2,1:N));
        I_simu_sum(k1,k2,:) = sum(I_simu(k1,k2,1:N));
    end
end