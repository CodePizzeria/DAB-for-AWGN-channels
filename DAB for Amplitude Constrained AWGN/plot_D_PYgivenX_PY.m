function [x, D_PYgivenX_PY, The_Max] = plot_D_PYgivenX_PY(figure_number, title_text, N, PY,PYgivenX, MI, xSupport, Maximizing_x_values)
x=0:max(xSupport)/1000:max(xSupport);
D_PYgivenX_PY= zeros(1, length(x));
for i=1:length(x)
D_PYgivenX_PY(i)= find_D_PYgivenX_PY(PY, PYgivenX, x(i), xSupport, N);
%n, PY, x(i), 10^(-10));
end
% [EmpiricalMax Empirical_Index] = max(D_PYgivenX_PY)
% EmpiricalDifference = EmpiricalMax - MI
figure(figure_number)
clf
plot(x, D_PYgivenX_PY, '.')
hold on 
plot( [0 max(xSupport)], [MI MI], 'g--')
plot(xSupport,  MI * ones(1,length(xSupport)),'o', 'MarkerFaceColor', 'g')
The_Max=0;
Dmax = zeros(1, length( Maximizing_x_values));
for i=1:length( Maximizing_x_values)
    Dmax(i) = find_D_PYgivenX_PY(PY, PYgivenX, Maximizing_x_values(i), xSupport, N);
    The_Max = max(The_Max, Dmax(i));
end
plot(Maximizing_x_values, Dmax, 'o', 'MarkerFaceColor', 'r')
plot_positions = mat2str(xSupport);
plot_positions = regexprep(plot_positions, ' ', '\n');
plot_positions = erase(plot_positions, '[');
plot_positions = erase(plot_positions, ']');
text(0.3,2.05,plot_positions)
grid on
title (title_text,  'FontSize', 16)
xlabel ('x', 'FontSize', 16)
ylabel ('D(Py|x||Py)', 'FontSize', 16)