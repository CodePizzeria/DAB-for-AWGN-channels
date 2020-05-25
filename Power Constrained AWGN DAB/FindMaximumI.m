
function [Maximizing_x_star, Maximum_I] = FindMaximumI(xSupport, InputPMF, index, N, E)

m = length(xSupport);
xmax = max(xSupport);
xmin = min(xSupport);
if index ~= 1 
    e = InputPMF.*xSupport.^2;
    if m == 4 
        P_in = 0;
        E_in = 0;
    else
        P_in = sum(InputPMF([index+1:m-index]));
        E_in = sum(e([index+1:m-index]));
    end

    xSearchLim = max(xSupport(index-1),-sqrt((E-E_in)/(1-P_in)));

    fun = @(Qx_star) -Find_symmetric_I_x_star_out(xSupport, InputPMF, index, Qx_star, N, E);
    [Maximizing_x_star, Fval, ExitFlag] = fminbnd(fun, xSearchLim, xSupport(index)+0.1);
else 
    xSearchLim = -40; %-1 for peak limit, really large peak limit for only avg power constraint
    if m == 3
        InputPMF([1,3]) = InputPMF([1,3])-1e-3;
        InputPMF(2) = InputPMF(2)+2e-3;
    end
    fun = @(Qx_star) -Find_symmetric_I_x_star_in(xSupport, InputPMF, index, Qx_star, N, E);
    [Maximizing_x_star, Fval, ExitFlag] = fminbnd(fun, xSearchLim, xSupport(index)+0.1);
end
% over here, xSearchLim was greater than xSupport(index)
Maximum_I =  -Fval;
end
