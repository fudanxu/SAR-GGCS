function p = Probability_N(A,N)
a = 1;
m = length(A);
p = zeros(1,m);
for n = 1:m
    if N > 3
      p(n) = abs(quadgk(@(rho)probabilityfun(rho,A(n),N,a),0,inf));
    else
      p(n) = abs(quadgk(@(rho)probabilityfun(rho,A(n),N,a),0,200));
%     p(n) = quadgk(@(rho)probabilityfun(rho,A(n),N,a),0,inf);
%     p(n) = integral(@(rho)probabilityfun(rho,A(n),N,a),0,inf);
%     p(n) = integral(@(rho)probabilityfun(rho,A(n),N,a),0,100);
      if N==1
          p(n) = p(n)/200;
      end
    end
end
end

function f = probabilityfun(rho,A,N,a)
J0 = besselj(0,2*pi*A*rho);
JN = besselj(0,2*pi*a*rho/sqrt(N)).^(N);
f = rho.*JN.*J0.*4.*(pi^2).*A;
end

