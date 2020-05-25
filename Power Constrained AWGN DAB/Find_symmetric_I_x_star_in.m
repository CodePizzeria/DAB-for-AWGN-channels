function[sMI] = Find_symmetric_I_x_star_in(xSupport, InputPMF, index, x_star, N, E)
% USE THIS TO DO THE SAME THING AS THE ORIGINAL FUNCTION, EXCEPT YOU MOVE
% AN OUTER POINT OUTWARDS, WHILE FLOWING PROBABILITY TO/FROM INTERNAL
% POINTS. 
m = length(xSupport);
xmax = max(xSupport);
xmin = min(xSupport);
% new code where we deal with rescaling pmf. 


xSupport(index) = x_star;
xSupport(m-index+1) = xmax-x_star+xmin;
e = InputPMF.*xSupport.^2;
P_in = sum(InputPMF([index+1:m-index]));
P_star = sum(InputPMF([index,m-index+1]));
P_out = sum(InputPMF([1:index-1,m-index+2:m]));

E_in = sum(e([index+1:m-index]));
E_star = sum(e([index,m-index+1]));
E_out = sum(e([1:index-1,m-index+2:m]));

if m == 3
    scale = E/E_star;
    scale_in = (1-scale*P_star)/P_in;
else
    scale = 1/P_in - E/E_in;
    scale = scale/(P_star/P_in-E_star/E_in);
    scale_in = 1/P_in - P_star/P_in*scale;
end

InputPMF([index,m-index+1]) = scale*InputPMF(index);
InputPMF(2:end-1) = scale_in*InputPMF(2:end-1);
% Compute Mutual Information
sMI = MutualInformation(InputPMF, xSupport, N);
end