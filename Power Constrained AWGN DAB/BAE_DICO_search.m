function[p_star, MI, EE, s] = BAE_DICO_search(pX, xSupport, N,E, ETolerance,BAETolerance,plot_flag)
    EE_plot = [];
    s_plot = [];
    s = 0;
    f = @(s) BAE_DICO(pX, xSupport, N,s,BAETolerance);
    [p_star, MI, EE, iter] = f(s);
    EE_plot = [EE_plot EE];
    s_plot = [s_plot s];
    if plot_flag, figure(2), hold on, end
    if EE < E
        fprintf('BAE finds p(x) with E = %f, constraint irrelevant', EE)
    else
        i = 1;
        s = 0.1;
        [p_star, MI, EE, iter] = f(s);
        EE_plot = [EE_plot EE];
        s_plot = [s_plot s];
        if plot_flag, plot(s_plot(end), EE_plot(end),'o'), end 
        while abs(EE-E) > ETolerance
            step_correction = 1/(1+20*0.1^i);
            i = i+1;
            ds = ( (EE_plot(end)-EE_plot(end-1)) / (s_plot(end)-s_plot(end-1)) );
            
            step = max(-1,min(1,step_correction*(EE-E)/ds));
            s = max(0, s - step);
            [p_star, MI, EE, iter] =f(s);
            EE_plot = [EE_plot EE];
            s_plot = [s_plot s];
            if s >4
                display("s broke")
            end
            if i > 100, break, end
            if plot_flag, plot(s_plot(end), EE_plot(end),'o'), drawnow, pause(1), end
        end
    end
    if plot_flag, hold off, pause(5), close(gcf), end
end 
 