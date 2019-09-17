function [R_simu_sum,I_simu_sum] = Clip_RIN(row,col,betax,betay,mu,s,hr,N_mean,N)
R_simu = normrnd(betax.*mu./N_mean,s./sqrt(N_mean),row,col,max(max(N)));
I_simu = normrnd(betay.*mu./N_mean,s./sqrt(N_mean),row,col,max(max(N)));
% ********************************generate the correlated gaussian scattering field********************************
h0 = waitbar(0,'Please wait...');
for k = 1:max(max(N))
    R_simu(:,:,k) = imfilter(R_simu(:,:,k),hr,'replicate'); 
    [m,sigma] = normfit(reshape(R_simu(:,:,k),1,[]));
    R_simu(:,:,k) = (R_simu(:,:,k)-m)./sigma;
    R_simu(:,:,k) = R_simu(:,:,k).*s./sqrt(N_mean)+betax.*mu./N_mean;
    I_simu(:,:,k) = imfilter(I_simu(:,:,k),hr,'replicate'); 
    [m,sigma] = normfit(reshape(I_simu(:,:,k),1,[]));
    I_simu(:,:,k) = (I_simu(:,:,k)-m)./sigma;
    I_simu(:,:,k) = I_simu(:,:,k).*s./sqrt(N_mean)+betay.*mu./N_mean;
    waitbar(k/max(max(N)),h0);
end
close(h0)
%% ********obtain the real image and imaginary image ************
R_simu_sum = zeros(row,col);
I_simu_sum = zeros(row,col);
% *********first N summation**********
for k1 = 1:row
    for k2 = 1:col
        R_simu_sum(k1,k2,:) = sum(R_simu(k1,k2,1:N(k1,k2)));
        I_simu_sum(k1,k2,:) = sum(I_simu(k1,k2,1:N(k1,k2)));
    end
end