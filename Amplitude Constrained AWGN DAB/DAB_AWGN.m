function [xSupport, m, InputPMF, MI, Iter, idiffs, MIs, MaxDs, x_maxs, focus_indices, active_xs] = DAB_AWGN(N, xSupport, Maximum_x,Minimum_x, BATolerance, ITolerance, plot_flag, Max_Iter)
% ==================== TO DO ===================
% for call to 
%   - pass BA tolerance top down 
%   - pass integral tolerance top down (+/- gaussian)
%   - fix integral problems at low noise 
%   - major speed up 
%   - save plotting variables in a vector and plot afterwards.
%   - fix how it looks like two points are being added at the same time
% removed 
%   - output Time
% ================== functions ===================
% [D_PYgivenX_PY] = find_D_PYgivenX_PY(PY, PYgivenX, x, xSupport, N)

% initialize storage vectors. 
special_pause = 2;
idiffs = zeros(1,Max_Iter);
MaxDs = zeros(1,Max_Iter);
MIs = zeros(1,Max_Iter);
focus_indices = zeros(1,Max_Iter);
x_maxs = zeros(1,Max_Iter);
active_xs = zeros(1,Max_Iter);

TSTART = tic;
m = length(xSupport);
xSupport_odd = mod(m,2);
Iter = 1;
% 1. finds p* and I using BA
[InputPMF,MI] = BA_DICO((1/m)*ones(m,1), xSupport, N);
% PY = find_q_from_p(xSupport,InputPMF',n);

% ===== 2. finds the location and value of maximum D(p(y|x)||p(y)), bookkeeping,
% ===== plot D(p(y|x)||p(y)) if plot_flag is on 
PY = @(y) InputPMF.'*gaussian(y,xSupport,N);
PYgivenX = @(y,x) gaussian(y,x,N);
[Maximizing_x_values, D_values, Maximizing_x_value, MaximumD] = findMaxD_PYgivenX_PY(PY, PYgivenX, xSupport, N);

Idifference= MaximumD - MI;
DBound = MaximumD;
x_maxs(Iter)= Maximizing_x_value;
MaxDs(Iter)= MaximumD;
MIs(Iter)= MI;
idiffs(Iter) = Idifference;
if(plot_flag)
    figure(75)
    title_text = sprintf('D(Py|x||Py)from iteration %d, N=%d, I difference = %1.7f', Iter ,10*log10(N), Idifference);
    [return_x, D_PYgivenX_PY, The_Max] = plot_D_PYgivenX_PY(75, title_text, N, PY, PYgivenX, MI, xSupport, Maximizing_x_value);
    pause(0.01);
end

split_flag = 0;

if(xSupport_odd)
    % m is length of support
    last_left_index = (m-1)/2;
else
    last_left_index = m/2;
end
% 3: 
while (Idifference>ITolerance) && (Iter < Max_Iter) && (split_flag==0)
    % Identify Focus Index
    % This is a value between 2 and middle index
    focus_index = 0;
    for check_index = 2:last_left_index
        % check_index
        if (Maximizing_x_value < xSupport(check_index)) && (Maximizing_x_value > xSupport(check_index -1))
            focus_index = check_index;
        elseif (Maximizing_x_value > xSupport(m - check_index + 1)) && (Maximizing_x_value < xSupport(m - check_index+2))
            focus_index = check_index;
        end
    end
    % focus_index
    focus_indices(Iter) = focus_index;
    if (focus_index == 0)
        if(plot_flag)
        figure(74)
        title_text = sprintf('BEFORE D(Py|x||Py) iter=%d, N=%d, I_diff=%1.7f', Iter ,10*log10(N), Idifference);
        [return_x, D_PYgivenX_PY, The_Max] = plot_D_PYgivenX_PY(74, title_text, N, PY, PYgivenX, MI, xSupport, Maximizing_x_value);
        pause(special_pause);
        end
        split_flag = 1;
        break;
    end
    
    % Now we use the function FindMaximumI to find the x_star that maximizes mutual information
    [Maximizing_x_star, Maximum_I] = FindMaximumI_via_fminbnd(xSupport, InputPMF, focus_index,N);
    
    % Now we reset Blahut Arimoto.
    xSupport(focus_index) = Maximizing_x_star;
    xSupport(m - focus_index+1) = Maximum_x - Maximizing_x_star +Minimum_x;
    active_xs(Iter) = xSupport(focus_index);
    
    Iter = Iter + 1;
    % we now compute the optimal probability allocation for these mass points
    [InputPMF, MI] = BA_DICO(InputPMF, xSupport, N);  %xSupport is a row vector
    MIs(Iter) = MI;
    PY = @(y) InputPMF.'*gaussian(y,xSupport,N);
    [Maximizing_x_values, D_values, Maximizing_x_value, MaximumD] = findMaxD_PYgivenX_PY(PY, PYgivenX, xSupport, N);
    %Lets find out how close this solution is to the relative entropy bound
    DBound = min(DBound,MaximumD);
    Idifference = DBound - MI;
    idiffs(Iter) = Idifference;
    MaxDs(Iter) = MaximumD;
    x_maxs(Iter)= Maximizing_x_value;
    if(plot_flag)
        figure(75)
        title_text = sprintf('D(Py|x||Py)from iteration %d, N=%d, I difference = %1.7f', Iter ,10*log10(N), Idifference);
        [return_x, D_PYgivenX_PY, The_Max] = plot_D_PYgivenX_PY(75, title_text, N, PY, PYgivenX, MI, xSupport, Maximizing_x_value);
        pause(0.01);
        
    end
end

% ===================== STOP HERE ===========================
% 6: Now we split the middle point or add a point at 1/2  if needed
if(split_flag==1)
    Maximizing_x_value = round(Maximizing_x_value, 10);
    % Code below assumes m is odd
    if(length(xSupport)~=3 || Maximizing_x_value ~= (Maximum_x+Minimum_x)/2)
        if(xSupport_odd)
            middle_index = (m+1)/2;
            middle_probability = InputPMF(middle_index);
            for i = m:-1:(middle_index)
                xSupport(i+1) = xSupport(i);
                InputPMF(i+1) = InputPMF(i);
            end
            xSupport(middle_index) = (Maximum_x+Minimum_x)/2;
            xSupport(middle_index+1) = (Maximum_x+Minimum_x)/2;
            InputPMF(middle_index) = middle_probability/2;
            InputPMF(middle_index+1) = middle_probability/2;
            m=m+1;
            xSupport_odd=0;
        else
            middle_index = m/2;
            for i = m:-1:(middle_index+1)
                xSupport(i+1) = xSupport(i);
                InputPMF(i+1) = InputPMF(i);
            end
            xSupport(middle_index+1) = (Maximum_x+Minimum_x)/2;
            m=m+1;
            xSupport_odd=1;
            for i=1:m
                InputPMF(i) = 1/m;
            end
        end
        % 6c
        if(xSupport_odd)
            last_left_index = (m-1)/2;
        else
            last_left_index = m/2;
        end
    end

%     xSupport
%     InputPMF
    
    Iter = Iter + 1;
    % we now compute the optimal probability allocation for these mass points
    [InputPMF, MI] = BA_DICO(InputPMF, xSupport, N);  %xSupport is a row vector
    MIs(Iter) = MI;
    PY = @(y) InputPMF.'*gaussian(y,xSupport,N);
    [Maximizing_x_values, D_values, Maximizing_x_value, MaximumD] = findMaxD_PYgivenX_PY(PY, PYgivenX, xSupport, N);
    %Lets find out how close this solution is to the relative entropy bound
    DBound = min(DBound,MaximumD);
    Idifference = DBound - MI;
    idiffs(Iter) = Idifference;
    MaxDs(Iter) = MaximumD;
    
if(plot_flag)
    figure(76)
    title_text = sprintf('SPLIT AFTER D(Py|x||Py) iter=%d, N=%d, I_diff = %1.7f', Iter ,10*log10(N), Idifference);
    [return_x, D_PYgivenX_PY, The_Max] = plot_D_PYgivenX_PY(76, title_text, N, PY, PYgivenX, MI, xSupport, Maximizing_x_value);
    pause(special_pause);
end 
end
% Now we continue regular iterations after the split.
while (Idifference>ITolerance) && (Iter < Max_Iter) && (split_flag==1)
    % Identify Focus Index
    % This is a value between 2 and middle index
    focus_index = 0;
    for check_index = 2:last_left_index
        %check_index
        if (Maximizing_x_value < xSupport(check_index)) && (Maximizing_x_value > xSupport(check_index -1))
            focus_index = check_index;
        elseif (Maximizing_x_value > xSupport(m - check_index + 1)) && (Maximizing_x_value < xSupport(m - check_index+2))
            focus_index = check_index;
        end
    end
    %focus_index
    focus_indices(Iter) = focus_index;
    if (focus_index == 0)
        split_flag = 2;
        break;
    end
    
    % Now we use the function FindMaximumI to finnd the x_star that maximizes mutual information
    [Maximizing_x_star, Maximum_I] = FindMaximumI_via_fminbnd(xSupport, InputPMF, focus_index,N);
    % problem is maxxstar is -1, and 1-it is equal to 2
    
    % Now we reset the xSupport vector to use the best positions for 2 and 4
    % and re-run Blahut Arimoto.
    
    xSupport(focus_index) = Maximizing_x_star;
    xSupport(m - focus_index+1) = Maximum_x - Maximizing_x_star + Minimum_x;
    active_xs(Iter) = xSupport(focus_index);
    Iter = Iter + 1;
    % we now compute the optimal probability allocation for these mass points
    [InputPMF, MI] = BA_DICO(InputPMF, xSupport, N);  %xSupport is a row vector
    MIs(Iter) = MI;
    PY = @(y) InputPMF.'*gaussian(y,xSupport,N);
    [Maximizing_x_values, D_values, Maximizing_x_value, MaximumD] = findMaxD_PYgivenX_PY(PY, PYgivenX, xSupport, N);
    DBound = min(DBound,MaximumD);
    Idifference = DBound - MI;
    idiffs(Iter) = Idifference;
    MaxDs(Iter) = MaximumD;
    x_maxs(Iter)= Maximizing_x_value;
    if(plot_flag)
        figure(75)
        title_text = sprintf('D(Py|x||Py)from iteration %d, N=%d, I difference = %1.7f', Iter ,10*log10(N), Idifference);
        [return_x, D_PYgivenX_PY, The_Max] = plot_D_PYgivenX_PY(75, title_text, N, PY, PYgivenX, MI, xSupport, Maximizing_x_value);
        pause(0.01);
    end
%     Maximizing_x_value
%     xSupport
%    Time = toc(TSTART);
end
