function c = K_correct_coef(alpha,N_mean)
% alpha = 1.5;
% N_mean = 80;
r = alpha;
p = alpha./(alpha+N_mean);
N1 = nbinrnd(r,p,1,10000);
N2 = nbinrnd(r,p,1,10000);
Nmin = min([N1;N2]);
c = mean(Nmin)/N_mean;
% fprintf('To correct : y = y/a and the corrected coefficient a = 1/c  is %.4f \n',1/c);
