function P = W_dis(beta,p,q,n,A)
% beta = 1;
% p = 1;
% q = 1;
% n = 1;
% A = 0:0.1:5;
temp1 = 2.*A.*gamma(p+q)./beta./gamma(n)./gamma(p);
temp2 = (A.^2./beta).^((p+n-3)./2).*exp(-1.*(A.^2)./2./beta);
temp3 = whittakerW((-p-2*q+n+1)./2,(p-n)./2,A.^2./beta);
P = temp1.*temp2.*temp3;
% figure; plot(A,P)