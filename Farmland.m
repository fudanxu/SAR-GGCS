clc
clear
close all
load('.\Data\farmland.mat')
load('.\Data\farmland_label.mat')
% *************************
addpath(genpath(pwd));

c1 = 25; c2 = 35;
classnum = max(max(farmland_label));
[row,col] = size(farmland);
L = zeros(row,col,classnum);
for k = 1:classnum
   ind = find(farmland_label==k);
   temp = zeros(size(farmland_label));
   temp(ind) = 1;
   L(:,:,k) = temp;
end

[row,col] = size(farmland);
R_simu_sum = zeros(row,col,classnum);
I_simu_sum = zeros(row,col,classnum);
N = 100;
for k = 1:classnum
    datak = L(:,:,k).*farmland;
    figure; imagesc(10.*log10(datak.*conj(datak)));axis tight equal; caxis([20,35])
    fprintf('Please select the homogenous samples for this kind of farmland'); 
    h = imrect;
    position = wait(h);
    rmin = ceil(position(2));
    rmax = rmin + fix(position(4));
    cmin = ceil(position(1));
    cmax = cmin + fix(position(3));
    hh = farmland(rmin:rmax,cmin:cmax);
    [R_simu_sum(:,:,k),I_simu_sum(:,:,k)]  = Simulation_correlated_Rayleigh_clutter_farmland(hh,N,row,col);
end
close all;
R_simu = sum(R_simu_sum.*L,3);
I_simu = sum(I_simu_sum.*L,3);
A_simu = sqrt(R_simu.^2+I_simu.^2);
za = sqrt(farmland.*conj(farmland));
theta = atan2(I_simu_sum,R_simu_sum);
zr = real(farmland);
zi = imag(farmland);
A = linspace(min(za(:)), max(za(:)), 50);
counts = hist(za(:),A);
P_A_real = counts./numel(za)./(A(end)-A(end-1));
[counts,centers] = hist(atan2(zi(:),zr(:)),100);
P_Theta_real = counts./numel(zi)./(centers(end)-centers(end-1));
counts = hist(A_simu(:),A);
P_A_simu = counts./numel(A_simu)./(A(end)-A(end-1));
[counts,centers] = hist(theta(:),100);
P_Theta_simu = counts./numel(theta)./(centers(end)-centers(end-1));
%% Plot figures
%% £¨a£©real data
figure; imagesc(20.*log10(za(1:row,1:col))); axis equal tight off; colorbar; colormap('gray'); caxis([c1,c2]);
tightfig;
print('-dtiff','-r300',['Results\','Farmland_real'])
%% £¨b£©simulated data
figure; imagesc(20.*log10(A_simu)); axis equal tight off; colorbar; colormap('gray'); caxis([c1,c2]);
tightfig;
print('-dtiff','-r300',['Results\','Farmland_simulated'])

