function P_A = RiIG_dis(alpha,beta,delta,A)
% alpha = 15;
% beta = 9;
% delta = 20;
% A = 0:0.1:30;
if (alpha^2-beta^2) < 0
    disp('Error: alpha must be bigger than beta');
end
gamma = sqrt(alpha.^2-beta.^2);
temp1 = sqrt(2/pi).*(alpha.^(3/2)).*delta.*exp(delta.*gamma);
temp2 = A./((delta.^2+A.^2).^(3/4));
temp3 = besselk(3/2,alpha.*sqrt(delta.^2+A.^2)).*besseli(0,beta.*A);
P_A = temp1.*temp2.*temp3;
% figure; plot(A,P_A)

