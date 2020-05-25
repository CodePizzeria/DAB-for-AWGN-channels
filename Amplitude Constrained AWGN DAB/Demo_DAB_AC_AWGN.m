clear
close all 
ITolerance = 10^(-4); % DAB converges if UB-LB is les than this. 
BATolerance = 10^(-7);% Blahut Arimoto runs until MI achieves is within 10^-7 of capacity
plot_flag=0;% plot_flag=0 turns off plotting
Max_Iterations = 500; % This is the maximum number of iterations.

Maximum_x = 1;
Minimum_x = -1;
XSupport = [Minimum_x Maximum_x].';
NdB = -3:-0.1:-17

% bookkeeping 
XSupportplot = [];
InputPMFplot = [];
mplot = [];
SNRplot = [];
MIplot = [];
Hplot = [];
H = @(pX) sum(pX.*log2(1./pX))
    
% run DAB for a whole range of SNRs. i is 

set(gcf, 'Position', get(0, 'Screensize'));

for i = 1:length(NdB)
    % run DAB for a specified SNR = 1/N
    N = 10^(NdB(i)/10);
    [XSupport, m, InputPMF, MI, Iterations, idiffs, MIs, MaxDs, x_maxs, focus_indices, active_xs] = DAB_AWGN(N, XSupport, Maximum_x, Minimum_x, BATolerance, ITolerance, plot_flag, Max_Iterations);

    % creating values to be plotted. 
    MIplot = [MIplot MI];
    E = sum(InputPMF.*(XSupport.^2));
    Hplot = [Hplot H(InputPMF)];
    SNR = 10*log10(E/N);
    SNRplot = [SNRplot SNR];
    m = length(XSupport);
    mplot = [mplot m];
    
    
    % let's plot them plots boissss
    for j = 1:m
        XSupportplot = [XSupportplot XSupport(j)];
        InputPMFplot = [InputPMFplot InputPMF(j)];
        figure(1)
        subplot(2,4,[1,2,5,6])
        hold on 
        plot(-NdB(i), XSupport(j), 'ko', 'MarkerSize', 25*sqrt(InputPMF(j))/2, 'MarkerFaceColor', 'r')
        hold off
    end
    if i ==1, title('Capacity Achieving Distributions for AC-AWGN Channel')
    legend('Area indicates probability','Location','northwest','AutoUpdate','off')
    xlabel('1/N (dB)')
    ylabel('Mass Point Locations')
    xlim(-[NdB(1),NdB(end)])
    end
       
    subplot(2,4,[3,4])
    hold on
    plot(-NdB(1:i), log2(mplot),'m')
    plot(-NdB(1:i), Hplot,'k')
    plot(-NdB(1:i), MIplot,'b')
    hold off 
    if i == 1, title("log2|X|, input entropy, and C_{peak}",'interpreter','tex')
    ylim([0,inf])
    xlim(-[NdB(1),NdB(end)])
    ylim([0,3])
    xlabel('peak SNR (dB)')
    ylabel('bits')
    legend('log_2(|X|)','input entropy','AC-AWGN capacity','Location','northwest','AutoUpdate','off')
    grid on, end
    
    subplot(2,4,[7,8])
    hold on
    plot(-NdB(1:i), SNRplot,'k')
    plot(-NdB, -NdB,'k--')
    hold off
    
    if i == 1, title("SNR vs peak SNR")
    legend('SNR', 'peak SNR','Location','northwest','AutoUpdate','off')
    xlabel('peak SNR, i.e. 1/N (dB)')
    ylabel('SNR')
    xlim(-[NdB(1),NdB(end)])
    ylim(-[NdB(1),NdB(end)]), end
    
    
    drawnow
end
