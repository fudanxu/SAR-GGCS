function [rho,R] = Cof_estimate_freq_1(z)
Z = ifftshift(fft2(fftshift(z-mean(z(:)))));
R = Z.*conj(Z);
rho = fftshift(ifft2(ifftshift(R)));
rho = rho./max(max(rho));
