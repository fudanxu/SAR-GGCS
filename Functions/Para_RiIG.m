function [alpha,beta,delta] = Para_RiIG(hh)
Amp = sqrt(hh.*conj(hh));
Amp_2 = mean(mean(Amp.^2));
% alpha = 8;
% beta = 1;
% delta = 20;
alpha = 0.1;
beta = 0.1;
delta = 1;
alphak = [alpha];
betak = [beta];
deltak = [delta];
Maxiter = 1000;
for k = 1:Maxiter
    gammak = sqrt(alpha^2-beta^2);
    vk = sqrt(delta^2+Amp.^2);
    muz = mean(mean(vk./alpha.*(vk.*alpha)./(1+vk.*alpha)));
    temp = alpha./vk.*((vk.*alpha).^2+3.*vk.*alpha+3)./(vk.*alpha.*(1+vk.*alpha));
    muz1 = mean(mean(temp));
    delta = sqrt(abs(1/muz1-1/muz));
    gammak1 = delta/muz;
    alpha = sqrt(gammak1^2+beta^2);
    beta = (Amp_2-2*delta/gammak1)/((gammak1^2/delta^2)+(delta/(gammak1^3)));
    fprintf('Iter: %.0f ; alpha: %.4f ; beta: %.4f ; delta: %.4f \n',k,alpha,beta,delta);
    alphak = [alphak;alpha];
    betak = [betak;beta];
    deltak = [deltak;delta];
end
figure; plot(alphak,'r','Linewidth',2) 
hold on; plot(betak,'b','Linewidth',2)
hold on; plot(deltak,'g','Linewidth',2)
legend({'\alpha','\beta','\delta'},'FontSize',15)
title(['Parameter estimation for RiIG',' \alpha=',num2str(alpha),' \beta=',num2str(beta),' \delta=',num2str(delta)],'FontSize',12)

