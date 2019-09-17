function P = U_dis(betac,p,q,n,A)
% betac = 1;
% p = 1;
% q = 1;
% n = 1;
% A = 0:0.1:5;
temp1 = 2.*A.*gamma(n+q)./betac./gamma(n)./beta(p,q);
temp2 = (A.^2./betac).^(n-1);
temp3 = kummerU(q+n,1+n-p,A.^2./betac);
P = temp1.*temp2.*temp3;
% figure; plot(A,P)