function P = SaSGR(gamma,alpha,R)
% gamma = 1;
% alpha = 1;
% R = 0:0.5:5;
f = @(s,r) r.*s.*exp(-(gamma.*s).^alpha).*besselj(0,s.*r);
m = length(R);
P = zeros(1,m);
for k = 1:m
    r = R(k);
    P(k) = quadgk(@(s)f(s,r),0,Inf,'MaxIntervalCount',1500);
end
