function [D_PYgivenX_PY] = find_D_PYgivenX_PY(PY, PYgivenX, x, xSupport, N)
% options
%   - pass the PY function handle, and p(Y|x) function handle 

%function [D_PYgivenX_PY] = find_D_PYgivenX_PY(n, PY, x)
% Computes the relative entropy between PYgivenX,
% which is the AWGN and the output distribution PY 
% induced by the selected input distribution on x.
PYgivenXeqx = @(y) PYgivenX(y,x);
pad = 5*sqrt(N);
xmin = min(xSupport);
xmax = max(xSupport);
arg = @(y) PYgivenXeqx(y).*log2(PYgivenXeqx(y))-PYgivenXeqx(y).*log2(PY(y));
D_PYgivenX_PY = integral(arg, xmin - pad, xmax + pad);
end
