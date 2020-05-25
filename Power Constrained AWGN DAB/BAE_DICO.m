function[p_star, MI, E,iter] = BAE_DICO(pX, xSupport, N,s, BAETolerance)
% WARNING:
%   WILL NOT WORK WHEN NOISE IS VERY LOW RELATIVE TO SPACING: IE 4 points
%   with MI = 1.73 bits works well, but anything higher leads to bad
%   numerical things. 
% TO DO:
%   try to decouple generation of the noise with the actual BA algorithm,
%   so that it can be run with any discrete input, continuous output
%   channel. 

% BA for Discrete Input Continuous Output
% in BlahutArimotoCMM, you can generate the matrix P(Y|X), and simply pass
% P(Y|X) and the probabilities P(X) to the BA algorithm. However, we have
% to change this a little, since there is no matrix in this case to pass

% We pass xSupport (the means), N (the variance), and pX, then
% calculate everything internally. We will want to generalize this algorithm
% to any discrete input, continuous output channel. but this has not been
% done yet. 

% xSupport is a column vector including the and pX are column vectors

% ---------START CODE----------
% constants 

n = length(xSupport);
e = xSupport.^2;
% outputs 
C = 0;
p_star = pX;

error = BAETolerance*2;
iter = 0;
C_data = [];
I_U_data = [];
while error > BAETolerance
    iter = iter + 1;
    if iter ~= 1
        p_star = p_star.*c/2^(I_L);
    end
    % note I_L = C when convergence happens. So we just use C for both
    c = getcj(p_star, xSupport, N, n);
    c = 2.^(log2(c)-s*e);
    % calculate I_L
    I_L = log2(p_star.'*c);
    % calculate I_U
    I_U = log2(max(c));
    error = I_U - I_L;
    if iter == 1000
        break
    end
end
E = p_star.'*e;
MI = s*E + I_L;
end