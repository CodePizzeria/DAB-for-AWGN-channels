function[p_star, MI,iter] = BA_DICO(pX, xSupport, N)
% BA for Discrete Input Continuous Output
% This is the well-known Blahut Arimoto Algorithm. it is outlined in
% Blahut's paper titled "Computation of Capacity and Rate Distortion Functi
% xSupport and pX are column vectors

% START CODE 
% inputs 
% constants 
target = 1e-10; %tolerance of upper and lower bound on capacity for BA (DAB)
n = length(xSupport);
% outputs 
C = 0;
p_star = pX;

error = target*2;
iter = 0;
while error > target
    iter = iter + 1;
    if iter ~= 1
        p_star = p_star.*c/2^(C);
    end
    % note C = I_L at convergence. So we call I_L C in this case.
    c = getcj(p_star, xSupport, N, n);
    % calculate I_L
    C = log2(p_star.'*c);
    % calculate I_U
    I_U = log2(max(c));
    error = I_U - C;
    if iter == 1000
        break
    end
end
MI = C;
end