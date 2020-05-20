function[sMI] = Find_symmetric_I_x_star(XSupport, InputPMF, index, x_star, N)
m = length(XSupport);
% Compute PYgivenX and especially PYgivenX_x_star and PYgivenX_1_minus_x_star
xmax = max(XSupport);
xmin = min(XSupport);

XSupport(index) = x_star;

XSupport(m-index+1) = xmax-x_star+xmin;
% Compute Mutual Information
sMI = MutualInformation (InputPMF, XSupport, N);



