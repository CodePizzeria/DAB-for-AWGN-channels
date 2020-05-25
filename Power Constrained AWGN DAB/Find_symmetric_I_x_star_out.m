function[sMI] = Find_symmetric_I_x_star_out(xSupport, InputPMF, index, x_star, N, E)
% checked: this calculates a new distribution where the average energy is
%   maintained by "stealing" probability from the points further out from the
%   points being moved 

m = length(xSupport);
xmax = max(xSupport);
xmin = min(xSupport);
% new code where we deal with rescaling pmf. 


xSupport(index) = x_star;
xSupport(m-index+1) = xmax-x_star+xmin;
e = InputPMF.*xSupport.^2;
if m == 4 
    P_in = 0;
else
    P_in = sum(InputPMF([index+1:m-index]));
end
P_star = sum(InputPMF([index,m-index+1]));
P_out = sum(InputPMF([1:index-1,m-index+2:m]));

if m == 4 
    E_in = 0;
else
    E_in = sum(e([index+1:m-index]));
end
E_star = sum(e([index,m-index+1]));
E_out = sum(e([1:index-1,m-index+2:m]));

scale = (1-P_in)/P_out - (E-E_in)/E_out;
scale = scale/(P_star/P_out-E_star/E_out);

scale_out = (1-P_in)/P_out - P_star/P_out*scale;

InputPMF([index,m-index+1]) = scale*InputPMF(index);
InputPMF([1:index-1,m-index+2:m]) = InputPMF([1:index-1,m-index+2:m])*scale_out;
% Compute Mutual Information
sMI = MutualInformation(InputPMF, xSupport, N);
end