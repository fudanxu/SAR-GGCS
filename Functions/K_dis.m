function P = K_dis(alpha,beta,n,A)
% alpha = 1;
% beta = 1;
% n = 1;
% A = 0:0.01:5;
temp1 = 4.*beta.*n.*A./(gamma(n).*gamma(alpha));
temp2 = (beta.*n.*(A.^2)).^((alpha+n)./2-1);
temp3 = besselk(alpha-n,2.*A.*sqrt(beta));
P = temp1.*temp2.*temp3;
% figure; plot(A,P)
