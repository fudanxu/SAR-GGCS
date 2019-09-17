%% £¨a£©Rayleigh distribution
clc
clear
close all
% parameters
sigma = 1;
N = 2;
num = 10000;
A = 0:0.1:3.5;
% Rayleigh pdf
b = sqrt(sigma./2);
P_A = raylpdf(A,b);
% model
s = sqrt(sigma./2);
R_simu = normrnd(0,s./sqrt(N),N,num); 
R_simu = sum(R_simu);
I_simu = normrnd(0,s./sqrt(N),N,num); 
I_simu = sum(I_simu);
A_simu = sqrt(R_simu.^2+I_simu.^2);
theta =  atan2(I_simu,R_simu);
counts = hist(A_simu,A);
P_A_simu = counts./numel(A_simu)./(A(end)-A(end-1));
[counts,centers] = hist(theta,50);
P_Theta_simu = counts./numel(theta)./(centers(end)-centers(end-1));
P_Theta = 1/2/pi.*ones(1,numel(centers));
fig('Rayleigh',A,P_A,P_A_simu,centers,P_Theta,P_Theta_simu,N,0,sigma,0,0,0,0,0,0,0);
%% £¨b£©SaSGR distribution
clc
clear 
close all
% parameters
gamma = 10;
alpha = 1.8;
N = 150;
num = 10000;
A = 0:2:80;
% SaSGR pdf
P_A = SaSGR(gamma,alpha,A);
% model
pd = makedist('Stable','alpha',alpha,'beta',0,'gam',gamma.*(N.^(-1./alpha)),'delta',0);
P_RI = pdf(pd,A);
R_simu = random(pd,N,num);
I_simu = random(pd,N,num);
R_simu = sum(R_simu);
I_simu = sum(I_simu);
A_simu = sqrt(R_simu.^2+I_simu.^2);
theta =  atan2(I_simu,R_simu);
counts = hist(A_simu,A);
P_A_simu = counts./numel(A_simu)./(A(end)-A(end-1));
[counts,centers] = hist(theta,50);
P_Theta_simu = counts./numel(theta)./(centers(end)-centers(end-1));
P_Theta = 1/2/pi.*ones(1,numel(centers));
fig('SaSGR',A,P_A,P_A_simu,centers,P_Theta,P_Theta_simu,N,0,0,alpha,gamma,0,0,0,0,0);
%% £¨c£©K distribution-1
clear
close all
% ²ÎÊý
alpha = 1.5;
beta = 3;
n = 1;
N_mean = 100;
num = 10000;
A = 0:0.06:2.5;
% K  pdf
P_A = K_dis(alpha,beta,n,A);
% model
sigma = alpha./beta;
s = sqrt(sigma./2);
r = alpha;
p = alpha./(alpha+N_mean);
N = nbinrnd(r,p,1,num);
R_simu_sum = zeros(1,num);
for k = 1:num
    R_simu = normrnd(0,s./sqrt(N_mean),1,N(k));
    R_simu_sum(k) = sum(R_simu);
end
I_simu_sum = zeros(1,num);
for k = 1:num
    I_simu = normrnd(0,s./sqrt(N_mean),1,N(k));
    I_simu_sum(k) = sum(I_simu);
end
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
theta =  atan2(I_simu_sum,R_simu_sum);
counts = hist(A_simu,A);
P_A_simu = counts./numel(A_simu)./(A(end)-A(end-1));
[counts,centers] = hist(theta,50);
P_Theta_simu = counts./numel(theta)./(centers(end)-centers(end-1));
P_Theta = 1/2/pi.*ones(1,numel(centers));
fig('K1',A,P_A,P_A_simu,centers,P_Theta,P_Theta_simu,0,N_mean,0,alpha,0,beta,0,0,0,0)
%% £¨d£©K distribution-2
clc
clear
close all
% parameters
alpha = 1.5;
beta = 3;
n = 1;
N_mean = 100;
num = 10000;
A = 0:0.06:2.5;
% K pdf
P_A = K_dis(alpha,beta,n,A);
% model
sigma = alpha./beta;
s = sqrt(sigma./2);
aa = alpha;
bb = N_mean/alpha;
lambda = gamrnd(aa,bb,1,num);
N = poissrnd(lambda,1,num);
N2 = nbinrnd(alpha,alpha./(alpha+N_mean),1,num);
[countsN,centers] = hist(N,100);
countsN2 = hist(N2,centers);
R_simu_sum = zeros(1,num);
for k = 1:num
    R_simu = normrnd(0,s./sqrt(N_mean),1,N(k));
    R_simu_sum(k) = sum(R_simu);
end
I_simu_sum = zeros(1,num);
for k = 1:num
    I_simu = normrnd(0,s./sqrt(N_mean),1,N(k));
    I_simu_sum(k) = sum(I_simu);
end
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
theta = atan2(I_simu_sum,R_simu_sum);
counts = hist(A_simu,A);
P_A_simu = counts./numel(A_simu)./(A(end)-A(end-1));
[counts,centers] = hist(theta,50);
P_Theta_simu = counts./numel(theta)./(centers(end)-centers(end-1));
P_Theta = 1/2/pi.*ones(1,numel(centers));
fig('K2',A,P_A,P_A_simu,centers,P_Theta,P_Theta_simu,0,N_mean,0,alpha,0,beta,0,0,0,0)

%% £¨e£©K distribution-3
clc
clear
close all
% parameters
alpha = 1.5;
beta = 3;
n = 1;
N = 100;
num = 10000;
A = 0:0.06:2.5;
% K pdf
P_A = K_dis(alpha,beta,n,A);
% model
aa = alpha;
b = 2*sqrt(beta);
bb = 4./b.^2;
z = gamrnd(aa,bb,1,num);
z = ones(N,1)*z;
R_simu = normrnd(0,sqrt(z./2./N),N,num);
R_simu = sum(R_simu);
I_simu = normrnd(0,sqrt(z./2./N),N,num);
I_simu = sum(I_simu);
A_simu = sqrt(R_simu.^2+I_simu.^2);
theta =  atan2(I_simu,R_simu);
counts = hist(A_simu,A);
P_A_simu = counts./numel(A_simu)./(A(end)-A(end-1));
[counts,centers] = hist(theta,50);
P_Theta_simu = counts./numel(theta)./(centers(end)-centers(end-1));
P_Theta = 1/2/pi.*ones(1,numel(centers));
fig('K3',A,P_A,P_A_simu,centers,P_Theta,P_Theta_simu,N,0,0,alpha,0,beta,0,0,0,0);
%% £¨f£© G^0 distribution
clc
clear
close all
% parameters
alpha = 3; 
gamma = 2;
n = 1;
N_mean = 100;
num = 10000;
A = 0:0.1:4 ;
% G0 pdf
P_A = G0_dis(gamma,alpha,n,A);
% model
sigma = gamma./(alpha-1);
s = sqrt(sigma./2);
aa = alpha;
bb = 1./N_mean./(alpha-1);
lambda = gamrnd(aa,bb,1,num);
lambda = 1./lambda;
N = poissrnd(lambda,1,num);
R_simu_sum = zeros(1,num);
for k = 1:num
    R_simu = normrnd(0,s./sqrt(N_mean),1,N(k));
    R_simu_sum(k) = sum(R_simu);
end
I_simu_sum = zeros(1,num);
for k = 1:num
    I_simu = normrnd(0,s./sqrt(N_mean),1,N(k));
    I_simu_sum(k) = sum(I_simu);
end
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
theta = atan2(I_simu_sum,R_simu_sum);
counts = hist(A_simu,A);
P_A_simu = counts./numel(A_simu)./(A(end)-A(end-1));
[counts,centers] = hist(theta,50);
P_Theta_simu = counts./numel(theta)./(centers(end)-centers(end-1));
P_Theta = 1/2/pi.*ones(1,numel(centers));
fig('G0',A,P_A,P_A_simu,centers,P_Theta,P_Theta_simu,0,N_mean,0,alpha,gamma,0,0,0,0,0)
%% £¨g£©W distribution
clc
clear
close all
% parameters
p = 5;
q = 9;
beta = 40;
n = 1;
N_mean = 100;
num = 10000;
A = 0:0.3:14 ;
% W pdf
P_A = W_dis(beta,p,q,n,A);
% model
sigma = p./(p+q).*beta;
s = sqrt(sigma./2);
lambda = betarnd(p,q,1,num);
betaN = N_mean.*(p+q)./(p);
lambda = betaN.*lambda;
N = poissrnd(lambda,1,num);
R_simu_sum = zeros(1,num);
for k = 1:num
    R_simu = normrnd(0,s./sqrt(N_mean),1,N(k));
    R_simu_sum(k) = sum(R_simu);
end
I_simu_sum = zeros(1,num);
for k = 1:num
    I_simu = normrnd(0,s./sqrt(N_mean),1,N(k));
    I_simu_sum(k) = sum(I_simu);
end
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
theta = atan2(I_simu_sum,R_simu_sum);
counts = hist(A_simu,A);
P_A_simu = counts./numel(A_simu)./(A(end)-A(end-1));
[counts,centers] = hist(theta,50);
P_Theta_simu = counts./numel(theta)./(centers(end)-centers(end-1));
P_Theta = 1/2/pi.*ones(1,numel(centers));
fig('W',A,P_A,P_A_simu,centers,P_Theta,P_Theta_simu,0,N_mean,0,0,0,beta,p,q,0,0)
%% £¨h£©U distribution
clc
clear
close all
% parameters
p = 5;
q = 9;
beta = 40;
n = 1;
N_mean = 100;
num = 10000;
A = 0:0.5:20 ;
% U pdf
P_A = U_dis(beta,p,q,n,A);
% model
sigma = p./(q-1).*beta;
s = sqrt(sigma./2);
lambda = betarnd(p,q,1,num);
lambda = lambda./(1-lambda);
betaN = N_mean.*(q-1)./p;
lambda = betaN.*lambda;
N = poissrnd(lambda,1,num);
R_simu_sum = zeros(1,num);
for k = 1:num
    R_simu = normrnd(0,s./sqrt(N_mean),1,N(k));
    R_simu_sum(k) = sum(R_simu);
end
I_simu_sum = zeros(1,num);
for k = 1:num
    I_simu = normrnd(0,s./sqrt(N_mean),1,N(k));
    I_simu_sum(k) = sum(I_simu);
end
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
theta = atan2(I_simu_sum,R_simu_sum);
counts = hist(A_simu,A);
P_A_simu = counts./numel(A_simu)./(A(end)-A(end-1));
[counts,centers] = hist(theta,50);
P_Theta_simu = counts./numel(theta)./(centers(end)-centers(end-1));
P_Theta = 1/2/pi.*ones(1,numel(centers));
fig('U',A,P_A,P_A_simu,centers,P_Theta,P_Theta_simu,0,N_mean,0,0,0,beta,p,q,0,0)
%% £¨i£©RiIG distribution-1
clc
clear
close all
% parameters
alpha = 15;
beta = 1;
betax = 0.8;
delta = 20;
N = 100;
num = 10000;
A = 0:0.2:7.5;
% RiIG pdf
P_A = RiIG_dis(alpha,beta,delta,A);
% model
betay = sqrt(beta^2-betax^2);
gamma = sqrt(alpha.^2-beta.^2);
lambda = delta^2;
mu = delta/gamma;
pd = makedist('InverseGaussian','mu',mu,'lambda',lambda);
z = random(pd,1,num);
z = ones(N,1)*z;
R_std = normrnd(0,1,N,num);
I_std = normrnd(0,1,N,num);
R_simu = normrnd(betax.*z./N,sqrt(z./N),N,num);
R_simu = sum(R_simu);
I_simu = normrnd(betay.*z./N,sqrt(z./N),N,num);
I_simu = sum(I_simu);
A_simu = sqrt(R_simu.^2+I_simu.^2);
theta =  atan2(I_simu,R_simu);
counts = hist(A_simu,A);
P_A_simu = counts./numel(A_simu)./(A(end)-A(end-1));
[counts,centers] = hist(theta,50);
P_Theta_simu = counts./numel(theta)./(centers(end)-centers(end-1));
P_Theta = [];
fig('RiIG1',A,P_A,P_A_simu,centers,P_Theta,P_Theta_simu,N,0,0,alpha,0,beta,0,0,0,delta)
%% £¨j£©RiIG distribution-2
clc
clear
close all
% parameters
alpha = 15;
beta = 1;
betax = 0.8;
delta = 20;
num = 10000;
N_mean = 100;
A = 0:0.2:7.5;
% RiIG pdf
P_A = RiIG_dis(alpha,beta,delta,A);
% model
betay = sqrt(beta^2-betax^2);
gamma = sqrt(alpha.^2-beta.^2);
mu = delta/gamma;
sigma = mu*2;
s = sqrt(sigma./2);
lambda = (N_mean.*gamma).^2;
pd = makedist('InverseGaussian','mu',N_mean,'lambda',lambda);
Lam = random(pd,1,num);
N = poissrnd(Lam,1,num);
for k = 1:num
    R_simu = normrnd(betax.*mu./N_mean,s./sqrt(N_mean),1,N(k));
    R_simu_sum(k) = sum(R_simu);
end
I_simu_sum = zeros(1,num);
for k = 1:num
    I_simu = normrnd(betay.*mu./N_mean,s./sqrt(N_mean),1,N(k));
    I_simu_sum(k) = sum(I_simu);
end
A_simu = sqrt(R_simu_sum.^2+I_simu_sum.^2);
theta =  atan2(I_simu_sum,R_simu_sum);
counts = hist(A_simu,A);
P_A_simu = counts./numel(A_simu)./(A(end)-A(end-1));
[counts,centers] = hist(theta,50);
P_Theta_simu = counts./numel(theta)./(centers(end)-centers(end-1));
P_Theta = [];
fig('RiIG2',A,P_A,P_A_simu,centers,P_Theta,P_Theta_simu,0,N_mean,0,alpha,0,beta,0,0,0,delta)

%% £¨k£©Rician distribution
clc
clear
close all
% parameters
sigma = 2;
r = 2;
N = 10;
num = 10000;
A = 0:0.2:8;
centers = -pi:0.2:pi;
% Rician pdf
s = sqrt(sigma./2);
A0 = r*s;
pd = makedist('Rician','s',A0,'sigma',s);
P_A = pdf(pd,A);
P_Theta = Rician_phase_dis(A0,s,centers);
% model
s = sqrt(sigma./2);
A0 = r*s;
R_simu = normrnd(0,s./sqrt(N),N,num); 
R_simu = sum(R_simu)+A0;
I_simu = normrnd(0,s./sqrt(N),N,num); 
I_simu = sum(I_simu);
A_simu = sqrt(R_simu.^2+I_simu.^2);
theta_simu =  atan2(I_simu,R_simu);
counts = hist(A_simu,A);
P_A_simu = counts./numel(A_simu)./(A(end)-A(end-1));
counts = hist(theta_simu,centers);
P_Theta_simu = counts./numel(theta_simu)./(centers(end)-centers(end-1));
fig('Rician',A,P_A,P_A_simu,centers,P_Theta,P_Theta_simu,N,0,sigma,0,0,0,0,0,r,0)

