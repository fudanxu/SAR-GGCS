function Ind = Location(Fzm,u)
sz = size(u);
u = reshape(u,1,[]);
Ind = zeros(size(u));
for k = 1:length(u)
    Fz_a = [Fzm,u(k)];
    [~,I] = sort(Fz_a);
    Ind(k) = find(I==length(Fz_a))-1;
end
Ind = reshape(Ind,sz);