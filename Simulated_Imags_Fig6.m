% Gaussian kernel
w = 5;
h = fspecial('gaussian',w,1.5);
%% £¨a£©Rayleigh distribution
clc
clear
close all
% parameters
sigma = 1;
N = 2;
row = 100;
col = 100;
num = row*col;
A = 0:0.1:10;
% pdf
b = sqrt(sigma./2);
P_A = raylpdf(A,b);
s = sqrt(sigma./2);
R_simu = normrnd(0,s./sqrt(N),row,col,N);
I_simu = normrnd(0,s./sqrt(N),row,col,N);
w = 5;
h = fspecial('gaussian',w,1.5);
for k = 1:max(max(N))
    R_simu(:,:,k) = imfilter(R_simu(:,:,k),h,'replicate'); 
    [mu,sigma] = normfit(reshape(R_simu(:,:,k),1,[]));
    R_simu(:,:,k) = (R_simu(:,:,k)-mu)./sigma;
    R_simu(:,:,k) = R_simu(:,:,k).*s./sqrt(N);
    I_simu(:,:,k) = imfilter(I_simu(:,:,k),h,'replicate'); 
    [mu,sigma] = normfit(reshape(I_simu(:,:,k),1,[]));
    I_simu(:,:,k) = (I_simu(:,:,k)-mu)./sigma;
    I_simu(:,:,k) = I_simu(:,:,k).*s./sqrt(N);
end
R_simu_sum = zeros(row,col);
I_simu_sum = zeros(row,col);
for k1 = 1:row
    for k2 = 1:col
        R_simu_sum(k1,k2,:) = sum(R_simu(k1,k2,1:N));
        I_simu_sum(k1,k2,:) = sum(I_simu(k1,k2,1:N));
    end
end
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
figure; imagesc(10.*log10(A_simu)); axis equal tight off; colormap('gray'); colorbar; caxis([-10,15])
tightfig;
print('-dtiff','-r300',['Results\','Rayleigh_h'])
%% £¨b£©SaSGR distribution
clc
clear
close all
% parameters
gamma = 10;
alpha = 1.8;
N = 150;
row = 100;
col = 100;
num = row*col;
A = 0:1:2000;
% SaSGR pdf
P_A = SaSGR(gamma,alpha,A);
pd = makedist('Stable','alpha',alpha,'beta',0,'gam',gamma.*(N.^(-1./alpha)),'delta',0);
P_RI = pdf(pd,A);
R_simu = random(pd,row,col,N);
I_simu = random(pd,row,col,N);
w = 5;
h = fspecial('gaussian',w,1.5);
for k = 1:max(max(N))
    R_simu(:,:,k) = imfilter(R_simu(:,:,k),h,'replicate'); 
    gam_old = gamma.*(N.^(-1./alpha));
    gam_new = (sum(sum(abs(h*gam_old).^(alpha)))).^(1/alpha);
    R_simu(:,:,k) = R_simu(:,:,k).*gam_old./gam_new;
    I_simu(:,:,k) = imfilter(I_simu(:,:,k),h,'replicate'); 
    I_simu(:,:,k) = I_simu(:,:,k).*gam_old./gam_new;
end
R_simu_sum = zeros(row,col);
I_simu_sum = zeros(row,col);
for k1 = 1:row
    for k2 = 1:col
        R_simu_sum(k1,k2,:) = sum(R_simu(k1,k2,1:N));
        I_simu_sum(k1,k2,:) = sum(I_simu(k1,k2,1:N));
    end
end
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
figure; imagesc(10.*log10(A_simu)); axis equal tight off; colormap('gray'); colorbar; caxis([5,30])
tightfig;
print('-dtiff','-r300',['Results\','SaSGR_h'])
%% £¨c£©K distribution
clc
clear
close all
% parameters
alpha = 10;
beta = 3;
n = 1;
N_mean = 100;
row = 100;
col = 100;
num = row*col;
A = 0:0.1:3;
% K  pdf
P_A = K_dis(alpha,beta,n,A);
sigma = alpha./beta;
s = sqrt(sigma./2);
r = alpha;
p = alpha./(alpha+N_mean);
N = nbinrnd(r,p,row,col);
R_simu = normrnd(0,s./sqrt(N_mean),row,col,max(max(N)));
I_simu = normrnd(0,s./sqrt(N_mean),row,col,max(max(N)));
w = 5;
h = fspecial('gaussian',w,1.5);
for k = 1:max(max(N))
    R_simu(:,:,k) = imfilter(R_simu(:,:,k),h,'replicate'); 
    [mu,sigma] = normfit(reshape(R_simu(:,:,k),1,[]));
    R_simu(:,:,k) = (R_simu(:,:,k)-mu)./sigma;
    R_simu(:,:,k) = R_simu(:,:,k).*s./sqrt(N_mean);
    I_simu(:,:,k) = imfilter(I_simu(:,:,k),h,'replicate'); 
    [mu,sigma] = normfit(reshape(I_simu(:,:,k),1,[]));
    I_simu(:,:,k) = (I_simu(:,:,k)-mu)./sigma;
    I_simu(:,:,k) = I_simu(:,:,k).*s./sqrt(N_mean);
end
R_simu_sum = zeros(row,col);
I_simu_sum = zeros(row,col);
for k1 = 1:row
    for k2 = 1:col
        R_simu_sum(k1,k2,:) = sum(R_simu(k1,k2,1:N(k1,k2)));
        I_simu_sum(k1,k2,:) = sum(I_simu(k1,k2,1:N(k1,k2)));
    end
end
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
figure; imagesc(10.*log10(A_simu)); axis equal tight off; colormap('gray'); colorbar; caxis([-10,10])
tightfig;
print('-dtiff','-r300',['Results\','K_h'])
%% £¨d£©G0 distribution
clc
clear
close all
% parameters
alpha = 20; % alpha>1
gam = 2;
n = 1;
N_mean = 100;
row = 100;
col = 100;
num = row*col;
A = 0:0.1:10 ;
% G0 pdf
P_A = G0_dis(gam,alpha,n,A);
sigma = gam./(alpha-1);
s = sqrt(sigma./2);
aa = alpha;
bb = 1./N_mean./(alpha-1);
lambda = gamrnd(aa,bb,row,col);
lambda = 1./lambda;
N = poissrnd(lambda,row,col);
R_simu = normrnd(0,s./sqrt(N_mean),row,col,max(max(N)));
I_simu = normrnd(0,s./sqrt(N_mean),row,col,max(max(N)));
w = 5;
h = fspecial('gaussian',w,1.5);
for k = 1:max(max(N))
    R_simu(:,:,k) = imfilter(R_simu(:,:,k),h,'replicate'); 
    [mu,sigma] = normfit(reshape(R_simu(:,:,k),1,[]));
    R_simu(:,:,k) = (R_simu(:,:,k)-mu)./sigma;
    R_simu(:,:,k) = R_simu(:,:,k).*s./sqrt(N_mean);
    I_simu(:,:,k) = imfilter(I_simu(:,:,k),h,'replicate'); 
    [mu,sigma] = normfit(reshape(I_simu(:,:,k),1,[]));
    I_simu(:,:,k) = (I_simu(:,:,k)-mu)./sigma;
    I_simu(:,:,k) = I_simu(:,:,k).*s./sqrt(N_mean);
end
R_simu_sum = zeros(row,col);
I_simu_sum = zeros(row,col);
for k1 = 1:row
    for k2 = 1:col
        R_simu_sum(k1,k2,:) = sum(R_simu(k1,k2,1:N(k1,k2)));
        I_simu_sum(k1,k2,:) = sum(I_simu(k1,k2,1:N(k1,k2)));
    end
end
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
figure; imagesc(10.*log10(A_simu)); axis equal tight off; colormap('gray'); colorbar; caxis([-15,5])
tightfig;
print('-dtiff','-r300',['Results\','G0_h'])
%% £¨e£©U distribution
clc
clear
close all
% parameters
p = 5;
q = 8;
beta = 50;
n = 1;
N_mean = 100;
row = 100;
col = 100;
num = row*col;
A = 0:1:50 ;
% U pdf
P_A = U_dis(beta,p,q,n,A);
sigma = p./(q-1).*beta;
s = sqrt(sigma./2);
lambda = betarnd(p,q,row,col);
lambda = lambda./(1-lambda);
betaN = N_mean.*(q-1)./p;
lambda = betaN.*lambda;
N = poissrnd(lambda,row,col);
R_simu = normrnd(0,s./sqrt(N_mean),row,col,max(max(N)));
I_simu = normrnd(0,s./sqrt(N_mean),row,col,max(max(N)));
w = 5;
h = fspecial('gaussian',w,1.5);
for k = 1:max(max(N))
    R_simu(:,:,k) = imfilter(R_simu(:,:,k),h,'replicate'); 
    [mu,sigma] = normfit(reshape(R_simu(:,:,k),1,[]));
    R_simu(:,:,k) = (R_simu(:,:,k)-mu)./sigma;
    R_simu(:,:,k) = R_simu(:,:,k).*s./sqrt(N_mean);
    I_simu(:,:,k) = imfilter(I_simu(:,:,k),h,'replicate'); 
    [mu,sigma] = normfit(reshape(I_simu(:,:,k),1,[]));
    I_simu(:,:,k) = (I_simu(:,:,k)-mu)./sigma;
    I_simu(:,:,k) = I_simu(:,:,k).*s./sqrt(N_mean);
end
R_simu_sum = zeros(row,col);
I_simu_sum = zeros(row,col);
for k1 = 1:row
    for k2 = 1:col
        R_simu_sum(k1,k2,:) = sum(R_simu(k1,k2,1:N(k1,k2)));
        I_simu_sum(k1,k2,:) = sum(I_simu(k1,k2,1:N(k1,k2)));
    end
end
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
figure; imagesc(10.*log10(A_simu)); axis equal tight off; colormap('gray'); colorbar; caxis([0,20])
tightfig;
print('-dtiff','-r300',['Results\','U_h'])
%% £¨f£©W distribution
clc
clear
close all
% parameters
p = 5;
q = 10;
beta = 50;
n = 1;
N_mean = 100;
row = 100;
col = 100;
num = row*col;
A = 0:1:50 ;
% W pdf
P_A = W_dis(beta,p,q,n,A);
sigma = p./(p+q).*beta;
s = sqrt(sigma./2);
lambda = betarnd(p,q,row,col);
betaN = N_mean.*(p+q)./(p);
lambda = betaN.*lambda;
N = poissrnd(lambda,row,col);
R_simu = normrnd(0,s./sqrt(N_mean),row,col,max(max(N)));
I_simu = normrnd(0,s./sqrt(N_mean),row,col,max(max(N)));
w = 5;
h = fspecial('gaussian',w,1.5);
for k = 1:max(max(N))
    R_simu(:,:,k) = imfilter(R_simu(:,:,k),h,'replicate'); 
    [mu,sigma] = normfit(reshape(R_simu(:,:,k),1,[]));
    R_simu(:,:,k) = (R_simu(:,:,k)-mu)./sigma;
    R_simu(:,:,k) = R_simu(:,:,k).*s./sqrt(N_mean);
    I_simu(:,:,k) = imfilter(I_simu(:,:,k),h,'replicate'); 
    [mu,sigma] = normfit(reshape(I_simu(:,:,k),1,[]));
    I_simu(:,:,k) = (I_simu(:,:,k)-mu)./sigma;
    I_simu(:,:,k) = I_simu(:,:,k).*s./sqrt(N_mean);
end
R_simu_sum = zeros(row,col);
I_simu_sum = zeros(row,col);
for k1 = 1:row
    for k2 = 1:col
        R_simu_sum(k1,k2,:) = sum(R_simu(k1,k2,1:N(k1,k2)));
        I_simu_sum(k1,k2,:) = sum(I_simu(k1,k2,1:N(k1,k2)));
    end
end
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
figure; imagesc(10.*log10(A_simu)); axis equal tight off; colormap('gray'); colorbar; caxis([0,20])
tightfig;
print('-dtiff','-r300',['Results\','W_h'])
