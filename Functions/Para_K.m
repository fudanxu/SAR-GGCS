function [alpha,Beta] = Para_K(hh)
Amp = sqrt(hh.*conj(hh));
m_half = mean(mean(Amp.^(1/2)));
m_1 = mean(mean(Amp));
n = 1;
temp1 = gamma(n+1/2)*gamma(n)/(gamma(n+1/4)^2);
funK = @(alpha)temp1*gamma(alpha+1/2)*gamma(alpha)/(gamma(alpha+1/4)^2)-(m_1/(m_half^2));
alpha = 0;
alpha0 = 1;
Beta = 0;
while alpha<=eps*10 || Beta<=eps*10 || Beta==inf
    alpha =  fzero(funK,alpha0);
%     alpha0 = alpha0 + 0.1;
    alpha0 = alpha0 + 1;
    Beta = (gamma(n+1/2)/gamma(n)/m_1*gamma(alpha+1/2)/gamma(alpha))^2/n;
end
