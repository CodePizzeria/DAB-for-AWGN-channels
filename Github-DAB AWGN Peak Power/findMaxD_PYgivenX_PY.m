function [Maximizing_x_values, D_values, Maximizing_x_value, MaximumD] = findMaxD_PYgivenX_PY(PY, PYgivenX, xSupport, N)
% function [Maxima] = FindMaximaD_PYgivenX_PY(PY, Maximum_x)

% 1: find sign change indices of derivative of D(p(y|x)||p(y)) 

xmax = max(xSupport);
xmin = min(xSupport);
x_values=xmin:.01:xmax;
d_dx_D_PYgivenX_PY = zeros(1, length(x_values));
for i = 1:length(x_values)
    d_dx_D_PYgivenX_PY(i)=find_d_dx_D_PYgivenX_PY(PY, PYgivenX, x_values(i), xSupport, N);
end
indices = signChange(d_dx_D_PYgivenX_PY);
% 2. find the exact values of x where d/dx changes using fzero
Maximizing_x_values = [0];
fun = @(xx) find_d_dx_D_PYgivenX_PY(PY, PYgivenX, xx, xSupport, N);
for j = 1:length(indices)
    a_maximizing_x=fzero(fun, [x_values(indices(j)),x_values(indices(j)+1)]);
    Maximizing_x_values = [Maximizing_x_values a_maximizing_x];
end
Maximizing_x_values = [Maximizing_x_values xmax];
% 3. find the maximum D(p(y|x)||p(y)) value from the prev list
MaximumD=0;
D_values = zeros(1, length( Maximizing_x_values));
for i=1:length( Maximizing_x_values)
    D_values(i) = find_D_PYgivenX_PY(PY, PYgivenX, Maximizing_x_values(i), xSupport, N);
    if (D_values(i) > MaximumD)
        MaximumD = D_values(i);
        Maximizing_x_value = Maximizing_x_values(i);
    end
end