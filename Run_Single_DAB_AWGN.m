clear
close all 
ITolerance = 10^(-3);
BATolerance = 10^(-5);% Blahut Arimoto runs until MI achieves is within 10^-7 of capacity
plot_flag=0;% plot_flag=0 turns off plotting
Max_Iterations = 500; % This is the maximum number of iterations.

Maximum_x = 1;
Minimum_x = -1;
XSupport = [Minimum_x Maximum_x].';
NdB = -4:-0.05:-25
imax = length(NdB);
% bookkeeping 
peakSNRplot = [];
XSupportplot = [];
InputPMFplot = [];
mplot = [];
SNRplot = [];
MIplot = [];

for i = 1:imax
    N = 10^(NdB(i)/10);
    [XSupport, m, InputPMF, MI, Iterations, idiffs, MIs, MaxDs, x_maxs, focus_indices, active_xs] = DAB_AWGN(N, XSupport, Maximum_x, Minimum_x, BATolerance, ITolerance, plot_flag, Max_Iterations);
    % live printing of input distribution for debugging
    displaytext = [XSupport XSupport InputPMF];
    displaytext(:,1)= NdB(i);
    VarNames = {'N', 'XSupport','PMF'};
    table(displaytext(:,1), displaytext(:,2),displaytext(:,3),'VariableNames',VarNames)
    % bookkeeping
    peakSNRplot = [peakSNRplot -NdB(i)];
    MIplot = [MIplot MI];
    E = sum(InputPMF.*(XSupport.^2));
    SNR = 10*log10(E/N);
    SNRplot = [SNRplot SNR];
    m = length(XSupport);
    mplot = [mplot m];
    % let's plot them plots boissss
    for j = 1:m
        XSupportplot = [XSupportplot XSupport(j)];
        InputPMFplot = [InputPMFplot InputPMF(j)];
        figure(1)
        hold on 
        plot(-NdB(i), XSupport(j), 'ko', 'MarkerSize', 25*sqrt(InputPMF(j))/2, 'MarkerFaceColor', 'r')
        hold off
    end
    figure(2)
    hold on 
    subplot(2,2,1)
    plot(-NdB(i), MI,'o')
    title("rate")
    hold off
    
    hold on
    subplot(2,2,2)
    plot(-NdB(i), m, 'o')
    title("number of bits")
    hold off 
    
    hold on
    subplot(2,2,[3,4])
    plot(-NdB(i), SNR, 'o')
    title("SNR")
    hold off
    
    drawnow
end
