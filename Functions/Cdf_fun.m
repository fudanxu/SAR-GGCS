function y = Cdf_fun(centers_x,Fx,x)
centers_gap = centers_x(2)-centers_x(1);
Ind = ceil((x-min(centers_x))./centers_gap);
y = zeros(size(x));
for k = 1:length(y)
    if Ind(k)==0
%         kl = (Fx(Ind(k)+1))/centers_gap;
%         y(k) = kl*x(k);
        y(k) = Fx(1);
    elseif Ind(k)==length(Fx)
        y(k) = Fx(end);
    elseif Ind(k)>=1 && Ind(k)<length(Fx)
        kl = (Fx(Ind(k)+1)-Fx(Ind(k)))/centers_gap;
        y(k) = Fx(Ind(k))+kl*(x(k)-centers_x(Ind(k)));
    end
end
% centers_y = zeros(size(centers_x));
% for k = 1:length(centers_y)
%     centers_y(k) = find()
% end