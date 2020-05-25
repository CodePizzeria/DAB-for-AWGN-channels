function [Maximizing_x_star, Maximum_I] = FindMaximumI_via_fminbnd(XSupport, InputPMF, index,N)
% function [Maximizing_x_star, Maximum_I] = FindMaximumI(n,XSupport, InputPMF, index, min_x_star, max_x_star)
% This function finds maximum mutual information that can be obtained by symmetrically
% moving a binomial distribution mass point, and the maximizing mass point
% location
myfun = @(QXSupport, QInputPMF, Qindex, Qx_star)  -Find_symmetric_I_x_star(QXSupport, QInputPMF, Qindex, Qx_star,N);
fun = @(Qx_star) myfun(XSupport, InputPMF, index, Qx_star);
[Maximizing_x_star, Fval, ExitFlag] = fminbnd(fun, XSupport(index-1), XSupport(index));
Maximum_I =  -Fval;