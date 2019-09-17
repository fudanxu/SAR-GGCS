function z = InvCdf_fun(Fzm,centers_z,u)
Ind = Location(Fzm,u);
z = zeros(size(u));
for k = 1:length(z)
    if Ind(k)==0
        z(k) = centers_z(1);
    elseif Ind(k)==length(centers_z)
        z(k) = centers_z(end);
    elseif Ind(k)>=1 && Ind(k)<length(centers_z)
        Fzm_gap = Fzm(Ind(k)+1)-Fzm(Ind(k));
        centers_gap = centers_z(Ind(k)+1)-centers_z(Ind(k));
        kl = centers_gap/Fzm_gap;
        z(k) = centers_z(Ind(k))+kl*(u(k)-Fzm(Ind(k)));
    end
end

