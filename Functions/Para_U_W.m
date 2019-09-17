function [beta,p,q] = Para_U_W(hh)
Iten = hh.*conj(hh);
mu1 = mean(mean(Iten));
mu2 = mean(mean(Iten.^2));
mu3 = mean(mean(Iten.^3));
mu4 = mean(mean(Iten.^4));
a = mu2/2/(mu1^2);
b = mu3/3/mu1/mu2;
m1 = mu1/gamma(1+1);
m2 = mu2/gamma(2+1);
m3 = mu3/gamma(3+1);
m4 = mu4/gamma(4+1);
Beta1 = m3^2/m2^3;
Beta2 = m4/m2^2;
A = Beta1*((Beta2+3)^2)/4/(4*Beta2-3*Beta1)/(2*Beta2-3*Beta1-6);
fprintf('The coefficient A is %.4f \n',A);
if 2*Beta2-3*Beta1-6==0
    fprintf('It is K distribution \n');
    alpha = 1/(a-1);
    beta = 1/(mu1*(a-1));
end 
if A==1
    fprintf('It is G0 distribution \n');
    alpha = (a+1)/a;
    gam = mu1;
end
if A<0
    fprintf('It is W distribution \n');
    beta = mu1*(b*(a-2)+a)/(2*a-b-1);
    p = 2*(a-b)/(b*(2-a)-a);
    q = 2*(a-b)*(a-1)*(b-1)/(2*a-b-1)/(2*b-a*b-a);
end
if A>1
    fprintf('It is U distribution \n');
    beta = mu1*(b*(a-2)+a)/(2*a-b-1);
    p = 2*(a-b)/(b*(2-a)-a);
    q = (4*a-3*b-1)/(2*a-b-1);
end
