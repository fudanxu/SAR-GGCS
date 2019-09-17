function [sigma,A0] = Para_Rice(hh)
I = hh.*conj(hh);
I_mean = mean(mean(I));
I2_mean = mean(mean(I.^2));
A0 = (abs(2*I_mean^2-I2_mean))^(1/4);
sigma = I_mean-A0^2;