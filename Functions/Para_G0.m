function [alpha,Gam] = Para_G0(hh)
Amp = sqrt(hh.*conj(hh));
m_half = mean(mean(Amp.^(1/2)));
m_1 = mean(mean(Amp));
n = 1;
temp1 = (gamma(n+1/4)^2)/gamma(n+1/2)/gamma(n);
funG = @(alpha)temp1*(gamma(-alpha-1/4)^2)/gamma(-alpha-1/2)/gamma(-alpha)/-((m_half^2)/m_1);
alpha = 0;
alpha0 = 0.5;
Gam = 0;
while alpha<=0 || Gam<=0 || Gam==inf
    alpha =  fzero(funG,alpha0);
    alpha0 = alpha0 + 0.1;
    Gam = (m_1*gamma(n)*gamma(-alpha)/gamma(n+1/2)/gamma(-alpha-1/2))^2*n;
end
