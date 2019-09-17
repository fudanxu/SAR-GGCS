%% 已知概率求分布函数
function cdf = cdf_pr(P,centers)
    n = max(size(centers));
    centers_gap = centers(2)-centers(1);
    cdf = zeros(1,n);
    cdf(1) = P(1)*centers_gap;
    for k = 2:n
        cdf(k) = P(k)*centers_gap+cdf(k-1);
    end
        
