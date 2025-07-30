% This script is used to generate the plot of example non-stationary
% fluctuation analysis shown in Figure 3H.

% peak-scaled non-stationary fluctuation analysis
    load('Fig3_NSFA_Data.mat')
    aligned_trace = Data.aligned_trace;
    i_peak = Data.i_peak;
    t_peak = Data.t_peak;

% configuration
    fs = 10;                % data aquisition frequency (kHz)
    nbin = 80;              % number of bins for mean trace amplitude
    Vm = -68;               % holding voltage
    Ecl = -15;              % Ecl 
    t_win = [-5 50];        % time window relative to max rise time for trace alignment    
    i_axis = [-800, 250];   % amplitude window of plot 
    scalex = 2;             % time (ms) for a scale bar
    scaley = 20;            % amplitdue (pA) for a scale bar

% an average trace and its peak amplitude and time
    mean_trace = struct;
    mean_trace.trace = mean(aligned_trace,2); % mean trace
    mean_trace.i_peak = unique(min(mean_trace.trace(round(abs(t_win(1))*0.95*fs):round(abs(t_win(1))*1.2*fs)))); % peak amplitude of the mean trace 
    mean_trace.t_peak = min(find(mean_trace.trace == mean_trace.i_peak))/fs; % time at peak amplitude of the mean trace
        
% generate mean trace scaled to peak of for each trace
    mean_trace.ps_trace = mean_trace.trace.*(i_peak./mean_trace.i_peak);

% generate variance around the mean
    var_trace = var(aligned_trace-mean_trace.ps_trace,0,2);

% get mean and variance in a set range of mean amplitude (10-90% of max)
    % set range of amplitude to fit 
    aoi.range = [0.1 0.9]; 
    
    % find the corresponding x coordinates (aoi= amplitude of interest; xoi= x point of interest)       
    aoi.max = mean_trace.trace(mean_trace.t_peak*fs)*aoi.range(1);
    aoi.min = mean_trace.trace(mean_trace.t_peak*fs)*aoi.range(2);
    
    xoi.lin = find(mean_trace.trace(mean_trace.t_peak*fs:end) < aoi.max & mean_trace.trace(mean_trace.t_peak*fs:end) > aoi.min);
    xoi.lin = xoi.lin+mean_trace.t_peak*fs-1;

    % restrict xoi to be continous from max
    x_broken = min(find(diff(xoi.lin) > 1)); 
    xoi.lin(x_broken+1:end) = [];

    % bin xoi by amplitude 
    aoi.bin = linspace(aoi.min, aoi.max, nbin);
    xoi.bin = {};
    
    for xx = 1:length(aoi.bin)-1
        xoi.bin{xx} = find(mean_trace.trace(xoi.lin) >= aoi.bin(xx) &...
                          mean_trace.trace(xoi.lin) <  aoi.bin(xx+1));
        xoi.bin{xx} = xoi.bin{xx} + min(xoi.lin) - 1;
    end
    xoi.bin = xoi.bin(~cellfun('isempty', xoi.bin)); % remove empty bins
  
    % get mean amp and variance for each bin
    nsfa = struct;
    for xx = 1:length(xoi.bin)
        nsfa.mean_bin(xx) = mean(mean_trace.trace(xoi.bin{xx}));
        nsfa.var_bin(xx)  = mean(var_trace(xoi.bin{xx}));
    end
                
% fit a quadratic function 
    nsfa.fitcoef = polyfit(-nsfa.mean_bin,nsfa.var_bin,2);
    nsfa.fitline = polyval(nsfa.fitcoef,[aoi.min:aoi.max]*-1);

% get N, i, gamma, P(open,peak) from coefficients
    nsfa.N  = -1/nsfa.fitcoef(1);
    nsfa.i  = nsfa.fitcoef(2);               % (pA)
    nsfa.g  = abs(nsfa.i/(Vm - Ecl))*1000;   % (pS)
    nsfa.Po = -mean_trace.i_peak/(nsfa.N * nsfa.i); % P(open, peak)
        
           
% figure: aligned traces and the mean trace 
    figure;    
    for xx = 1:size(aligned_trace,2)
        h = plot([1:size(aligned_trace,1)]/fs,aligned_trace(:,xx),'Color',[0.2 0.2 0.2]); hold on;
        h.Color(4) = 0.3;
    end   
    plot([1:size(mean_trace.trace,1)]/fs,mean_trace.trace,'Color','g','linewidth',2) % mean trace
    plot(xoi.lin/fs,mean_trace.trace(xoi.lin),'Color','r','linewidth',2) % mark time window of analysis 
    xlabel('Time (ms)')
    ylabel('Amplitude (pA)')
    box off
    axis off
    xlim([0  25])   

    % Scale bar 
    scalebar_locx = xlim;
    scalebar_locy = ylim;
    scalex_loc(1) = (scalebar_locx(2)-scalebar_locx(1))*0.6;  % time scale (left)
    scalex_loc(2) = scalex_loc(1)+scalex;                     % time scale (right)
    scaley_loc(1) = scalebar_locy(1)+abs(scalebar_locy(2)-scalebar_locy(1))*0.6;  % amplitude scale (low)
    scaley_loc(2) = scaley_loc(1)+scaley;                                         % amplitude scale (high)
    plot([scalex_loc(1); scalex_loc(2)], [scaley_loc(1); scaley_loc(1)], '-k', ...             % time
         [scalex_loc(2); scalex_loc(2)], [scaley_loc(1); scaley_loc(2)], '-k', 'LineWidth', 2) % amplitude

    % text
    tx_scalex = text(scalex_loc(1)+scalex*0.5, scaley_loc(1)-abs(i_axis(2)-i_axis(1))*0.018,...    % time text location
                [num2str(scalex) ' ms'], 'HorizontalAlignment','center','FontSize',18);            % time text

    tx_scaley = text(scalex_loc(2)+(scalebar_locx(2)-scalebar_locx(1))*0.035,scaley_loc(1)+scaley*0.5,...    % amplitude text location
                [num2str(scaley) ' pA'], 'HorizontalAlignment','left','FontSize',18);                        % amplitude text
                         
    set(gcf,'PaperUnits','Inch','PaperPosition',[0 0 3 2.8])
    saveas(gcf,'Fig3H_NSFA_traces.svg')

% figure: variance-amplitude plot 
    figure;
    plot(-nsfa.mean_bin,nsfa.var_bin,'o','MarkerSize',8,'Color','k'); hold on;
    plot([aoi.min:aoi.max]*-1,nsfa.fitline,'Color','r','LineWidth',2);
    yl = ylim;
    xlim([0 50])
    ylim([0 yl(2)])
    xlabel('Amplitude (pA)')
    ylabel('Variance (\sigma^2)')
    box off
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 2;
    ax.XTick = [0:10:50];
    ax.YTick = [0:10:50];
   
    set(gcf,'PaperUnit','Inch','PaperPosition',[0 0 3.5 2.8]) 
    saveas(gcf,'Fig3H_VarAmp_plot.svg')

