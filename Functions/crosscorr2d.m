function R=crosscorr2d(f,g)

% Compute the 2D spatial autocorrelation of a matrix or image I using the 
% Wiener - Khintchine Theorem. The output is the normalized correlation
% coefficient -1 < C < 1.
%

%
% ref: http://mathworld.wolfram.com/Wiener-KhinchinTheorem.html
%
% %EX 1: 2D zero mean Gaussian process with variance = 1 
% %The result is a delta-function at zero, with spurious correlations
% around zero.
%
% I=randn(501,501);
% A=autocorr2d(I);
%

f=double(f); %convert to double
f=f-mean(f(:)); %subtract mean
f=f/sqrt(sum(f(:).*conj(f(:)))); %normalize magnitude

g=double(g); %convert to double
g=g-mean(g(:)); %subtract mean
g=g/sqrt(sum(g(:).*conj(g(:)))); %normalize magnitude

F = fft2(f);
G = fft2(g);

R = fftshift(ifft2(F.*conj(G))); %compute crosscorrelation