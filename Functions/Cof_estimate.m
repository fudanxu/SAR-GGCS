% M = 200;
% N = 300;
% g = wgn(M,N,0,1,1);
% z = abs(g);
% n_v = 5;
function rho = Cof_estimate(z,n_v)
M = size(z,1);
N = size(z,2);
n_c = fix(N/(2*n_v));
n_f = fix(M/(2*n_v));
x = zeros(n_v,n_v,n_f*n_c);
for p = 0:(n_f-1)
    for q = 0:(n_c-1)
        x(:,:,p*n_c+q+1) = z(2*n_v*p+1:n_v*(2*p+1),2*n_v*q+1:n_v*(2*q+1));
    end
end
m = mean(x,3);
mm = m(:,:,ones(1,n_f*n_c));
s = sqrt(sum((x - mm).*conj(x-mm),3));
xx = x - mm;
cov_z = zeros(n_v^2,n_v^2);
for n = 1:n_f*n_c
    temp = reshape(xx(:,:,n),[],1);
    cov_z = cov_z + temp*temp';
end
s_z = reshape(s,[],1);
s_z = s_z*s_z';
Rho = cov_z./s_z;
rho = reshape(Rho(1,:),n_v,n_v);
% figure; imagesc(real(rho));axis equal tight off; 

% c = cell(n_f,n_c);
% for p = 0:(n_f-1)
%     for q = 0:(n_c-1)
%         c{p+1,q+1} = z(2*n_v*p+1:2*n_v*(p+1),2*n_v*q+1:2*n_v*(q+1));
%     end
% end
% z_v = cell(n_f,n_c);
% for p = 0:(n_f-1)
%     for q = 0:(n_c-1)
%         z_v{p+1,q+1} = z(2*n_v*p+1:n_v*(2*p+1),2*n_v*q+1:n_v*(2*q+1));
%     end
% end