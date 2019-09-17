% 生成相关均匀分布的相位
function [X,Y] = Phase_uniform(row,col,hp)
X = normrnd(0,1,row,col);
Y = normrnd(0,1,row,col);
C = X+1i.*Y;
C = conj(imfilter(C,hp,'replicate')); 
X = real(C); Y = imag(C);
[mu,sig] = normfit(X(:));
X = (X-mu)./sig;
[mu,sig] = normfit(Y(:));
Y = (Y-mu)./sig;
% C = X+1i.*Y;