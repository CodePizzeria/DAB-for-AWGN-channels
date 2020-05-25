close all 
dBs = 0:1:30;
ms = [2,4,8,16];
startgap = .15;
endgap = 1e-2;%1e-2 usually
tolerance = 1e-6;
ETolerance = 1e-8; % 1e-6 orignally 
BAETolerance = 1e-8;
plotting = true;
maxIters = 200;


N = 0.2;

figure(1)
plot(dBs, 1/2*log2(1+10.^(dBs/10)),'k')
grid on
position = get(0, 'Screensize');
position(3) = position(3)/2;
set(gcf, 'Position', position);
title('Rate Achieved by DAB optimated input PMFs')
xlabel('SNR (dB)')
ylabel('Rate (bits)')



endpt = 3;
for m = ms
    
    xSupport = linspace(-endpt,endpt,m).';
    xmax = max(xSupport);
    xmin = min(xSupport);
    
    InputPMF = 1/m*ones(m,1);

    inputPMFs = [];
    xSupports = [];
    Unger_MI_plot = [];
    i_plot = [];
    
    for dB = dBs
        E = 10.^(dB/10)*N;
        MI_plot = [];
        for i = 1:maxIters
            if m == 2
                xSupport(1) = -sqrt(E);
                xSupport(2) = sqrt(E);
                InputPMF = 1/m*ones(m,1);
                MI = MutualInformation(InputPMF, xSupport, N);
                MI_plot = [MI_plot MI];
                break
            end
            % run BAE_DICO_search to find the optimal probability assignment on the current support set 
            [InputPMF, MI, EE, s] = BAE_DICO_search(InputPMF, xSupport, N,E, ETolerance,BAETolerance, false);
            MI_plot = [MI_plot MI]
            i_plot = [i_plot i];
            if plotting 
                figure(2)
                hold on
                for j = 1:m
                    plot(i, xSupport(j),'ko', 'MarkerSize', 25*sqrt(InputPMF(j))/2+1e-10, 'MarkerFaceColor', 'r')
                end
                hold off
                title('Finding Cardinality ' +string(m) +' PC-AWGN Channel Capacity approaching input PMF for SNR = ' + string(dB) + ' dB')
                xlabel('Iteration')
                ylabel('Mass Point Location')
                position = get(0, 'Screensize');
                position([1,2]) = position([3,4])/2;
                position([3,4]) = position([3,4])/2;
                set(gcf, 'Position', position);
                drawnow, end
            
            % move each point and its symmetric point. 
            for index = 1:floor(m/2)
                [Maximizing_x_star, Maximum_I] = FindMaximumI(xSupport, InputPMF, index, N, E);
                xSupport(index) = Maximizing_x_star;
                xSupport(m-index+1) = xmax - (xSupport(index)- xmin);
            end

            % I am not updating the p_star here, because the next BAE iteration
            % will scrap it and find the new p_star anyways.
            if i ~= 1,if abs(MI_plot(end) - MI_plot(end-1)) < tolerance, break, end,end
            if i == maxIters, disp("reached max iterations"), end
        end
        
        inputPMFs = [inputPMFs InputPMF];
        xSupports = [xSupports xSupport];
        Unger_MI_plot = [Unger_MI_plot MI_plot(end)];
        
        dBplotend = find(dBs == dB)
        figure(1)
        hold on 
        plot(dBs(1:dBplotend), Unger_MI_plot,'b')
        hold off
        legend('PC-AWGN Channel Capacity','Rate Achieved by DAB optimized input PMFs','AutoUpdate','off','Location','northwest')

        figure(2), clf
        position = get(0, 'Screensize');
        position([1,2]) = position([3,4])/2;
        position([3,4]) = position([3,4])/2;
        set(gcf, 'Position', position);

        figure(3)
        subplot(2,2,log2(m))
        xlim([0,max(dBs)])
        ylim([-10,10])
        hold on
        for j = 1:m
            plot(dB, xSupport(j),'ko', 'MarkerSize', 25*sqrt(InputPMF(j))/2+1e-10, 'MarkerFaceColor', 'r')
        end
        hold off
        position = get(0, 'Screensize');
        position(1) = position(3)/2;
        position([3,4]) = position([3,4])/2;
        set(gcf, 'Position', position);
        title('DAB Optimized input PMFs with cardinality ' + string(m))
        xlabel('SNR (dB)')
        ylabel('Mass Point Locations')
        drawnow

        if log2(m)-Unger_MI_plot(end) < endgap, dBs_end = find(dBs == dB);break; end
    end 
    C = 1/2*log2(1+10.^(dBs(1:dBs_end)/10));
    new_dBs_start = find(C - Unger_MI_plot > startgap);
    new_dBs_start = new_dBs_start(1);
    
    figure(2), clf
    dBs = dBs(new_dBs_start:end);
    
    endpt = -xSupports(1,new_dBs_start)*1.5;
end

