% 生成单个散射体的相关瑞利幅度分布
function ai = corrected_Amp(row,col,bi,ha)
% sigma_ai = sigma./N;
% b = sqrt(sigma_ai./2);
G = wgn(row,col,0);
F = ifftshift(fft2(fftshift(G)));
Ha = ifftshift(fft2(fftshift(ha)));
GF = Ha.*F;
G2 = ifftshift(ifft2(fftshift(GF)));
G2(imag(G2)<1e-10) = real(G2(imag(G2)<1e-10));

ai0 = raylrnd(bi,row,col);
Lai0 = 20.*log10(ai0);
[counts,centers_ai] = hist(Lai0(:),100);
Pai = counts./numel(Lai0)./(centers_ai(end)-centers_ai(end-1));
Fai = cdf_pr(Pai,centers_ai);
[Lai,~] = DisTrans_to_Fz(Fai,centers_ai,G2);
ai = 10.^(Lai./20);

