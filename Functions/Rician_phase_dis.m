function P_theta = Rician_phase_dis(A0,s,theta)
% r = 1;
% s = 1;
% A0 = r*s;
% theta = -pi:0.01:pi;
temp1 = exp(-1.*(A0.^2)./(2.*s.^2))./2./pi;
temp2 = sqrt(1./(2.*pi)).*A0./s.*exp(-1.*(A0.^2).*(sin(theta).^2)./(2.*s.^2));
temp3 = (1+erf(A0.*cos(theta)./(sqrt(2).*s))).*cos(theta)./2;
P_theta = temp1 + temp2.*temp3;
% figure; plot(theta,P_theta)
