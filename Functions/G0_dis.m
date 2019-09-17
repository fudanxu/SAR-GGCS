function P = G0_dis(gam,alpha,n,A)
% alpha = 1;
% gam = 1;
% n = 1;
% A = 0:0.01:5;
temp1 = 2.*(n.^n).*gamma(n+alpha).*(gam.^alpha).*(A.^(2.*n-1));
temp2 = gamma(n).*gamma(alpha).*((gam+n.*(A.^2)).^(n+alpha));
P = temp1./temp2;
% figure; plot(A,P)