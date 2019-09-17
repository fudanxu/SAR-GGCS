
function fig(dis,A,P_A,P_A_simu,centers,P_Theta,P_Theta_simu,N,N_mean,sigma,alpha,gamma,beta,p,q,r,delta)
figure; 
subplot(121)
% plot(A,P_A,'-.','Color',[0 0.447 0.741],'Linewidth',1)
% hold on; plot(A,P_A_simu,'Color',[0.85 0.33 0.1],'Linewidth',1)
bar(A,P_A_simu,'Linewidth',0.5); axis square
hold on; plot(A,P_A,'Linewidth',2);
legend({'Simulation','Theory'},'FontSize',10,'Location','northeast');
titleAmp(dis,N,N_mean,sigma,alpha,gamma,beta,p,q,r,delta);
ylim([0,max([P_A(:);P_A_simu(:)])]); xlim([min(A(:)),max(A(:))]);
set(gca,'FontSize',10)

if strcmp(dis,'RiIG1')||strcmp(dis,'RiIG2')
    subplot(122);
%     plot(centers,P_Theta_simu,'Color',[0.85 0.33 0.1],'Linewidth',1)
    bar(centers,P_Theta_simu,'Linewidth',0.5); axis square
    legend({'Simulation'},'FontSize',10,'Location','northeast')
    titlePhase(dis,N,N_mean,sigma,alpha,gamma,beta,p,q,r,delta);
    ylim([0,0.6]);
    xlim([-pi,pi]);
    set(gca,'FontSize',10)
else
    subplot(122);
%     plot(centers,P_Theta,'-.','Color',[0 0.447 0.741],'Linewidth',1)
%     hold on; plot(centers,P_Theta_simu,'Color',[0.85 0.33 0.1],'Linewidth',1)
    bar(centers,P_Theta_simu,'Linewidth',0.5); axis square
    hold on; plot(centers,P_Theta,'Linewidth',2)
    legend({'Simulation','Theory'},'FontSize',10,'Location','northeast');
    titlePhase(dis,N,N_mean,sigma,alpha,gamma,beta,p,q,r,delta);
%     ylim([0,1]);
    ylim([0,0.5]);
    xlim([-pi,pi]);
    set(gca,'FontSize',10)
end
tightfig;
print('-dtiff','-r300',['Results\',dis])

function titleAmp(dis,N,N_mean,sigma,alpha,gamma,beta,p,q,r,delta)
if strcmp(dis,'Rayleigh')
    title(['Amplitude N =',num2str(N),' \sigma=',num2str(sigma)],'FontSize',10);
elseif strcmp(dis,'SaSGR')
    title(['Amplitude N =',num2str(N),' \alpha=',num2str(alpha),' \gamma=',num2str(gamma)],'FontSize',10)
elseif strcmp(dis,'K1')
    title(['Amplitude N_m =',num2str(N_mean),' \alpha=',num2str(alpha),' \beta=',num2str(beta)],'FontSize',10)
elseif strcmp(dis,'K2')
    title(['Amplitude N_m^\prime = ',num2str(N_mean),' \alpha=',num2str(alpha),' \beta=',num2str(beta)],'FontSize',10)
elseif strcmp(dis,'K3')
    title(['Amplitude N =',num2str(N),' \alpha=',num2str(alpha),' \beta=',num2str(beta)],'FontSize',10)
elseif strcmp(dis,'G0')
    title(['Amplitude N_m^\prime =',num2str(N_mean),' \alpha=',num2str(alpha),' \gamma=',num2str(gamma)],'FontSize',10)
elseif strcmp(dis,'W')
    title(['Amplitude N_m^\prime =',num2str(N_mean),' p=',num2str(p),' q=',num2str(q),' \beta=',num2str(beta)],'FontSize',10)
elseif strcmp(dis,'U')
    title(['Amplitude N_m^\prime =',num2str(N_mean),' p=',num2str(p),' q=',num2str(q),' \beta=',num2str(beta)],'FontSize',10)
elseif strcmp(dis,'Rician')
    title(['Amplitude N =',num2str(N),' \sigma=',num2str(sigma),' r=',num2str(r)],'FontSize',10)
elseif strcmp(dis,'RiIG1')
    title(['Amplitude N =',num2str(N),' \alpha=',num2str(alpha),' \beta=',num2str(beta),' \delta=',num2str(delta)],'FontSize',10)
elseif strcmp(dis,'RiIG2')
    title(['Amplitude N_m^\prime =',num2str(N_mean),' \alpha=',num2str(alpha),' \beta=',num2str(beta),' \delta=',num2str(delta)],'FontSize',10)
end
    
function titlePhase(dis,N,N_mean,sigma,alpha,gamma,beta,p,q,r,delta)
if strcmp(dis,'Rayleigh')
    title(['Phase N =',num2str(N),' \sigma=',num2str(sigma)],'FontSize',10);
elseif strcmp(dis,'SaSGR')
    title(['Phase N =',num2str(N),' \alpha=',num2str(alpha),' \gamma=',num2str(gamma)],'FontSize',10)
elseif strcmp(dis,'K1')
    title(['Phase N_m =',num2str(N_mean),' \alpha=',num2str(alpha),' \beta=',num2str(beta)],'FontSize',10)
elseif strcmp(dis,'K2')
    title(['Phase N_m^\prime =',num2str(N_mean),' \alpha=',num2str(alpha),' \beta=',num2str(beta)],'FontSize',10)
elseif strcmp(dis,'K3')
    title(['Phase N =',num2str(N),' \alpha=',num2str(alpha),' \beta=',num2str(beta)],'FontSize',15)
elseif strcmp(dis,'G0')
    title(['Phase N_m^\prime =',num2str(N_mean),' \alpha=',num2str(alpha),' \gamma=',num2str(gamma)],'FontSize',10)
elseif strcmp(dis,'W')
    title(['Phase N_m^\prime =',num2str(N_mean),' p=',num2str(p),' q=',num2str(q),' \beta=',num2str(beta)],'FontSize',10)
elseif strcmp(dis,'U')
    title(['Phase N_m^\prime =',num2str(N_mean),' p=',num2str(p),' q=',num2str(q),' \beta=',num2str(beta)],'FontSize',10)
elseif strcmp(dis,'Rician')
    title(['Phase N =',num2str(N),' \sigma=',num2str(sigma),' r=',num2str(r)],'FontSize',15)
elseif strcmp(dis,'RiIG1')
    title(['Phase N =',num2str(N),' \alpha=',num2str(alpha),' \beta=',num2str(beta),' \delta=',num2str(delta)],'FontSize',10)
elseif strcmp(dis,'RiIG2')
    title(['Phase N_m^\prime =',num2str(N_mean),' \alpha=',num2str(alpha),' \beta=',num2str(beta),' \delta=',num2str(delta)],'FontSize',10)
end