function h0 = Select_center_win(h,w)
[r,c] = find(h ==max(max(h)));
if mod(w,2)==0
    h0 = h(r-fix(w/2)+1:r+fix(w/2),c-fix(w/2)+1:c+fix(w/2));
elseif mod(w,2)==1
    h0 = h(r-fix(w/2):r+fix(w/2),c-fix(w/2):c+fix(w/2));
end