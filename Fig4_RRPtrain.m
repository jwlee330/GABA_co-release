% This script is used to generate the plots shown in Figure 4A-C. 
    genotype = 'dKO';   % options: 'WT','65KO','dKO'

    nDataPts = 50;      % number of first data points to use for analysis
    nRegression = 10;   % number of last data points to use for regression
    fs = 10;            % sampling frequency (kHz)
    
    load('Fig4_RRPtrain_Data.mat')
    mean_trace = Data.mean_trace{strcmp(Data.genotype,genotype)};
    stim_xloc = Data.stim_xloc{strcmp(Data.genotype,genotype)};
    peak_amp = Data.peak_amp{strcmp(Data.genotype,genotype)};
    mIPSCAmp = Data.mIPSCAmp(strcmp(Data.genotype,genotype));

% amplitude from the averaged trace
    nStim      = [0:nDataPts-1]';
    absAmp     = peak_amp/1000;  %(nA)
    absVesicle = peak_amp/mIPSCAmp;
    cumAmp     = cumsum(absAmp);
    cumVesicle = cumsum(absVesicle);
    
% fit a linear model
    [ampLine.fit, ampLine.stat] = fit(nStim(end-10+1:end),cumAmp([nStim(end-10+1:end)]+1),'poly1');
    [vesLine.fit, vesLine.stat] = fit(nStim(end-10+1:end),cumVesicle([nStim(end-10+1:end)]+1),'poly1');   
    fitLineAbs  = ampLine.fit.p1 * nStim + ampLine.fit.p2;
    fitLineVesicle = vesLine.fit.p1 * nStim + vesLine.fit.p2;
  
% set figure color
    switch genotype
        case 'WT';   linecolor = 'k'; 
        case '65KO'; linecolor = [0.6 0.6 0.6];  
        case 'dKO';  linecolor = [0.0000, 0.4470, 0.7410]; 
    end

% plot trace
    for xx = 1:nDataPts
        mean_trace(stim_xloc(xx):stim_xloc(xx)+10) = nan;
    end

    figure;
    subplot('Position',[0.05 0.05 0.9 0.8])
    f = plot(mean_trace,'-'); hold on;
    f.Color = linecolor;
    f.LineWidth = 1.5;          
    xlim([0 length(mean_trace)])
    ylim([-4600 50])
    axis off

    % Scale bar 
    scalex = 50;   % Time (ms)
    scaley = 500;  % Amplitdue (pA)
    xl = xlim; yl = ylim;
    scalex_loc(1) = (xl(2)-xl(1))*0.65;           % time scale (left)
    scalex_loc(2) = scalex_loc(1)+scalex*fs;      % time scale (right)
    scaley_loc(1) = yl(1)+abs(yl(2)-yl(1))*0.15;  % amplitude scale (low)
    scaley_loc(2) = scaley_loc(1)+scaley;         % amplitude scale (high)

    plot([scalex_loc(1); scalex_loc(2)], [scaley_loc(1); scaley_loc(1)], '-k', ...               % time
         [scalex_loc(2); scalex_loc(2)], [scaley_loc(1); scaley_loc(2)], '-k', 'LineWidth', 1.5) % amplitude

    % text
    tx_scalex = text(scalex_loc(1)+scalex*5, scaley_loc(1)-abs(yl(2)-yl(1))*0.04,...    % time text location
                [num2str(scalex) ' ms'], 'HorizontalAlignment','center','FontSize',20);         % time text

    tx_scaley = text(scalex_loc(2)+(xl(2)-xl(1))*0.01,scaley_loc(1)+scaley*0.5,...      % amplitude text location
                [num2str(scaley/1000) ' nA'], 'HorizontalAlignment','left','FontSize',20);           % amplitude text
    
    % save figure
    set(gcf, 'PaperUnits', 'centimeter','PaperPosition',[0 0 12 8]);
    saveas(gcf,['Fig4A_RRPtrain_trace_' genotype '.svg'])

% plot amplitude
    figure;
    subplot('position',[0.3 0.3 0.65 0.65])
    plot(nStim,cumAmp(1:nDataPts),'ok'); hold on;
    plot(nStim,fitLineAbs,'-','LineWidth',2.5, 'Color',linecolor);
    xlim([0 50])
    ylim([0 100])
    xlabel('Stimulation Number')
    ylabel('Cumulative IPSC (nA)')
    set(gca,'XTick',0:10:50,'LineWidth',2,'FontSize',18')
    box off
    set(gcf, 'PaperUnits', 'centimeter','PaperPosition',[0 0 12 12]);
    saveas(gcf,['Fig4B_RRPtrain_IPSC_' genotype '.svg'])
    
% plot vesicles
    figure;
    subplot('position',[0.3 0.3 0.65 0.65])
    plot(nStim,cumVesicle(1:nDataPts),'ok'); hold on;
    plot(nStim,fitLineVesicle,'-','LineWidth',2.5, 'Color',linecolor);
    xlim([0 50])
    ylim([0 3000])
    xlabel('Stimulation Number')
    ylabel({'Cumulative Number of'; 'Vesicle Release'})
    set(gca,'XTick',0:10:50,'YTick',0:500:5000,'LineWidth',2,'FontSize',18')
    box off
    set(gcf, 'PaperUnits', 'centimeter','PaperPosition',[0 0 12 12]);
    saveas(gcf,['Fig4B_RRPtrain_vesicle_' genotype '.svg'])
 
% result
    RRPtrain = struct;
    RRPtrain.RRP = vesLine.fit.p2;
    RRPtrain.Pr = cumAmp(1)/ampLine.fit.p2;
    RRPtrain.replenishment_percent = (ampLine.fit.p1*1000/mIPSCAmp)/RRPtrain.RRP*100;
    disp(RRPtrain)