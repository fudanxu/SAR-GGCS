%% SASGR分布参数估计
function [alpha,Gam] = Para_SaSGR(hh,r)
% r = -0.8;
if r>=-1/2 || r<=-2
    fprintf('Wrong r value, r should be -2<r<-1/2')
end

Amp = sqrt(hh.*conj(hh));
Amp_r = mean(mean(Amp.^(r)));
Amp_2r = mean(mean(Amp.^(2*r)));
temp1 = gamma(r+1)*(gamma(-r/2))^2;
temp2 = 2*gamma(-r)*(gamma(r/2+1))^2;
temp = Amp_2r/(Amp_r)^2*temp2/temp1;
funSaSGR = @(alpha)alpha*gamma(-2*r/alpha)/gamma(-r/alpha)-temp;
% options = optimset('Display','iter');
% a =  fzero(funSaSGR,1.6,options);
alpha = 0;
alpha0 = 0.1;
Gam = 0;
while alpha<=0 || alpha>2 || Gam<=0 || Gam==inf
% while alpha<=0 || alpha>2 || Gam<=0
    alpha =  fzero(funSaSGR,alpha0);
    alpha0 = alpha0 + 0.1;
    temp0 = gamma(-r/2)/(2^(r+1))/gamma(r/2+1)*alpha/gamma(-r/alpha);
    Gam = (Amp_r*temp0)^(alpha/r);
end

