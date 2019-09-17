%% 该函数存在问题尚待解决
function [alpha,Beta] = Para_K_Log(hh)
Amp = sqrt(hh.*conj(hh));
k1 = mean(mean(log(Amp)));
k2 = mean(mean((log(Amp)-k1).^2));
% k3 = mean(mean((log(Amp)-k1).^3));
funK1 = @(alpha)psi(1,1)+psi(1,alpha)-4*k2;
[alpha0,~] = Para_K(hh);
% alpha0 = 0.1;
alpha =  fzero(funK1,alpha0);
Beta = exp(2*k1-psi(alpha)-psi(1));