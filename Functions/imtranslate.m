function outimage = imtranslate(rho_h,rs,cs)
[r,c] = size(rho_h);
outimage = zeros(r,c);
if cs >= 0 && rs >= 0
    for m = 1:r-rs
        for n = 1:c-cs
            outimage(m+rs,n+cs) = rho_h(m,n);
        end
    end
elseif cs < 0 && rs < 0
    for m = 1-rs:r
        for n = 1-cs:c
            outimage(m+rs,n+cs) = rho_h(m,n);
        end
    end
elseif cs >= 0 && rs < 0
    for m = 1-rs:r
        for n = 1:c-cs
            outimage(m+rs,n+cs) = rho_h(m,n);
        end
    end
 elseif cs < 0 && rs >= 0
    for m = 1:r-rs
        for n = 1-cs:c
            outimage(m+rs,n+cs) = rho_h(m,n);
        end
    end   
end