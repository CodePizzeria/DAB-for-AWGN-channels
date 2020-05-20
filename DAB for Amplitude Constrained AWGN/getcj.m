function[c_vector] = getcj(pX, xSupport, N,n)
% TO DO:
%   set more dynamic limits for the integral 
c_vector = zeros(n, 1);
figaussian = @(x,mu,var) 1/sqrt(2*pi*var)*exp(-(x-mu).^2/2/var);
Q = @(y,x) gaussian(y, x, N);
sumpQ = @(y) pX.'*Q(y,xSupport);

for i = 1:n
    xi = xSupport(i);
    arg = @(y) Q(y,xi).*log2(Q(y,xi)./sumpQ(y)); 
    c_vector(i) = integral(arg, -1-10*sqrt(N),1+10*sqrt(N));
end
c_vector = 2.^(c_vector);
end
        