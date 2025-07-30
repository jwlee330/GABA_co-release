%% Part 1
% This script is used to generate the traces and plots of MNTB-evoked responses shown in Figure 2B.

load('Fig2_MNTBInputs_Data.mat','Data_PSCs')

for xx = 1:size(Data_PSCs,1)
% get data for the select genotype
    genotype = Data_PSCs.genotype{xx};
    rec = Data_PSCs.rec{xx};
    stim = Data_PSCs.stim{xx};
    peakamp = Data_PSCs.peakamp{xx};
    GMM = Data_PSCs.GMM(xx);

% plot MNTB-evoked response traces 
    fs = 10;                   % data aquisition frequency (kHz)
    t_axis = [96.8, 136.8];    % time window to plot (ms)
    i_axis = [-7000, 250];     % amplitude window to plot (pA)
    switch genotype
        case 'WT'; fcolor = [0 0 0];
        case 'dKO'; fcolor = [0.0000, 0.4470, 0.7410];
        case '65KO'; fcolor = [0.6 0.6 0.6];
    end

    % plot
    figure;
    plot(rec,'Color',[fcolor, 0.4],'LineWidth',1.5); hold on;
    plot(1000, i_axis(2), 'v','Color',[0.6350, 0.0780, 0.1840],'MarkerFaceColor',[0.6350, 0.0780, 0.1840],'MarkerSize',10) % add stimulus time
    axis([t_axis*fs, i_axis]);
    set(gca, 'Visible', 'off','color', 'none')

    % add scale bar 
    scalex = 5;     % time (ms)
    scaley = 500;   % amplitdue (pA)

    scalex_loc(1) = (t_axis(2)-t_axis(1))*0.6*fs;             % time scale (left)
    scalex_loc(2) = scalex_loc(1)+scalex*fs;                  % time scale (right)
    scaley_loc(1) = i_axis(1)+abs(i_axis(2)-i_axis(1))*0.7;   % amplitude scale (low)
    scaley_loc(2) = scaley_loc(1)+scaley;                     % amplitude scale (high)
    plot([scalex_loc(1), scalex_loc(2)], [scaley_loc(1), scaley_loc(1)], '-k', ...             % time
         [scalex_loc(2), scalex_loc(2)], [scaley_loc(1), scaley_loc(2)], '-k', 'LineWidth', 2) % amplitude

    % text on scale bar
    tx_scalex = text(scalex_loc(1)+scalex*5, scaley_loc(1)-abs(i_axis(2)-i_axis(1))*0.06,...    
                [num2str(scalex) ' ms'], 'HorizontalAlignment','center','FontSize',18);                     
    tx_scaley = text(scalex_loc(2)+(t_axis(2)-t_axis(1))*0.35,scaley_loc(1)+scaley*0.7,...      
                [num2str(scaley/1000) ' nA'], 'HorizontalAlignment','left','FontSize',18);                  

    % save      
    set(gcf, 'PaperUnits', 'centimeter','PaperPosition',[0 0 8.5 8.5]);
    saveas(gcf,['Fig2B_traces_' genotype '.svg'])

% fit GMM to peak amplitudes (commented out to use already existing data)
    % GMM = fitGMM(peakamp);
    
% generate a GMM plot
    figure;
    a_axis = [0 7]; % common amplitude axis

    % scatter plot on left
    subplot('position',[0.22 0.20 0.53 0.75])
    f1 = scatter(stim,peakamp/1000,30,fcolor);
    f1.MarkerFaceAlpha = 0.5;
    ax1 = gca;
    set(ax1,'fontsize',18,'LineWidth',2,'XTickLabelRotationMode','manual','YLim',a_axis,'YTick',0:1:7)  
    switch genotype
        case 'WT'; ax1.XTick = 0:200:600; xlim([0 600]); 
        case 'dKO'; ax1.XTick = 0:200:600; xlim([0 600]);
        case '65KO'; ax1.XTick = 0:100:200; xlim([0 200]);
    end
    xlabel('Stim. Intensity (\muA)');
    ylabel('Peak Amplitude (nA)');
    box off
    
    % histogram on right
    subplot('position',[0.75 0.20 0.19 0.75])
    f2 = histogram(peakamp/1000,0:0.15:a_axis(2),...
            'normalization','probability','orientation','horizontal'); hold on;
    f2.EdgeColor = [1 1 1];
    f2.FaceColor = fcolor;

    pdfx   = [0:max(peakamp)*1000*1.1]';
    pdfy   = pdf(GMM.GMMbest,pdfx);
    pdfy   = pdfy/max(pdfy)*max(f2.Values);

    f3 = plot(pdfy,pdfx/1000,'Color',fcolor,'LineWidth',2);
    ax3 = gca;
    ax3.Visible = 'off';
    ax3.YLim = a_axis;
    ax3.XTickLabel = [];
    ax3.YTickLabel = [];

    % mark the center of each Gaussian
    subplot('position',[0.94 0.20 0.02 0.75])
    
    ydata = GMM.GMMbest_mu/1000;
    xdata = ones(1,length(ydata));       
    
    f4 = plot(xdata,ydata,'*','Color',fcolor,'MarkerSize',8);
    ax4 = gca;
    ax4.YLim = a_axis;
    ax4.Visible = 'off';

    % save
    set(gcf, 'PaperUnits', 'centimeter','PaperPosition',[0 0 8.8 9.5]);
    saveas(gcf,['Fig2B_GMMplot_' genotype '.svg'])
end

%% Part 2
% This script is used to generate the distribution of fiber strengths shown in Fig 2F and 2G.

% load data
    load('Fig2_MNTBInputs_Data.mat','Data_FiberStrengths')

% process data
    group = {'WT','65KO','dKO'}; 
    FS = Data_FiberStrengths;
    for xx = 1:size(FS,1)
        FS.strongest{xx} = max(FS.fiber_strengths{xx});
        FS.nonstrongest{xx} = FS.fiber_strengths{xx}(FS.fiber_strengths{xx}~=FS.strongest{xx});
        FS.norm_strengths{xx} = FS.fiber_strengths{xx}/FS.strongest{xx};
        FS.norm_nonstrongest{xx} = FS.norm_strengths{xx}(FS.norm_strengths{xx} ~= 1);
    end
    
    % sort by strongest fiber strength for each genotype
    for genotype = group
        idx = cellfun(@(x) strcmp(x,genotype),FS.genotype);
        FS(idx,:) = sortrows(FS(idx,:),'strongest');
    end

    % group by genotype
    fiber_strengths   = cell(1,numel(group));
    strongest         = cell(1,numel(group));
    norm_strengths    = cell(1,numel(group));
    norm_nonstrongest = cell(1,numel(group));
    
    for gp = 1:numel(group)
        idx = find(cellfun(@(x) strcmp(x,group(gp)),FS.genotype));
        for xx = 1:numel(idx)
            add_strengths         = [ones(numel(FS.fiber_strengths{idx(xx)}),1)*(xx-1)+1, FS.fiber_strengths{idx(xx)}'];
            add_strongest         = [ones(numel(FS.strongest{idx(xx)}),1)*(xx-1)+1, FS.strongest{idx(xx)}'];
            add_norm_strengths    = [ones(numel(FS.norm_strengths{idx(xx)}),1)*(xx-1)+1, FS.norm_strengths{idx(xx)}'];
            add_norm_nonstrongest = [ones(numel(FS.norm_nonstrongest{idx(xx)}),1)*(xx-1)+1, FS.norm_nonstrongest{idx(xx)}'];
            
            fiber_strengths{gp}   = [fiber_strengths{gp}; add_strengths]; 
            strongest{gp}         = [strongest{gp}; add_strongest];
            norm_strengths{gp}    = [norm_strengths{gp}; add_norm_strengths];
            norm_nonstrongest{gp} = [norm_nonstrongest{gp}; add_norm_nonstrongest];
        end
    end

% figure setting
    % color
    pColor = {'k',[0.6 0.6 0.6],[0.0000, 0.4470, 0.7410]}; 

    % fontsize
    fontsz.tick   = 16;
    fontsz.label  = 18;
    fontsz.legend = 14;

% absolute fiber strength distribution
    figure;
    subplot_pos{1} = [0.17 0.48 0.80 0.50];
    subplot_pos{2} = [0.17 0.16 0.80 0.27];
    
    % raster plot of fiber strengths on top
    subplot('Position',subplot_pos{1})
    plot([0 5],[-25 -25],'k:','LineWidth',1); hold on;
    plot([0 5],[-41 -41],'k:','LineWidth',1);
    
    spacing = 1;
    for gp = 1:numel(group)
        plot(fiber_strengths{gp}(:,2)/1000,fiber_strengths{gp}(:,1)*-1-spacing,'o','MarkerSize',2,'Color',pColor{gp}); hold on;
        plot(strongest{gp}(:,2)/1000,strongest{gp}(:,1)*-1-spacing,'o','MarkerSize',2,'Color',pColor{gp},'MarkerFaceColor',pColor{gp}); % mark strongest fiber
        height = [min(fiber_strengths{gp}(:,1)*-1-spacing),max(fiber_strengths{gp}(:,1)*-1-spacing)];
        text(4,mean(height),group{gp},'FontSize',fontsz.legend,'HorizontalAlignment','left','VerticalAlignment','middle','Color',pColor{gp})
        spacing = spacing + size(strongest{gp},1)+3;
    end
    box on;
    ax = gca;
    ax.YTick = [-64, -39, -23,-2]; 
    ax.YTickLabel = {'57', '36', '23','1'}; 
    ax.XTickLabel = [];
    ax.YAxis.TickLength = [0 0];
    
    xlim([0 5]); ylim([-66,0])
    set(gca,'FontSize',14,'LineWidth',1.5)
    ylabel('Neuron Number','FontSize',fontsz.label)
    
    % cumulative plot on bottom
    subplot('Position',subplot_pos{2})
    cdf_amp = cell(1,numel(group));
    cdf_p   = cell(1,numel(group));
    for gp = 1:numel(group)
        [cdf_p{gp}, cdf_amp{gp}] = ecdf(fiber_strengths{gp}(:,2));
        plot(cdf_amp{gp}/1000,cdf_p{gp},'LineWidth',2,'Color',pColor{gp}); hold on;
    end
    
    xlim([0 5]); ylim([0 1])
    set(gca,'FontSize',fontsz.tick,'LineWidth',1.5); 
    xlabel('Fiber Strength (nA)','FontSize',fontsz.label); 
    ylabel('Probability','FontSize',fontsz.label); 
    box off; hold off

    set(gcf,'PaperUnits','Centimeter','PaperPosition',[0 0 11 11])
    saveas(gcf,'Fig2FG_FiberStrengths.svg')

% normalized fiber strength distribution
    figure;
    subplot_pos{1} = [0.17 0.48 0.80 0.50];
    subplot_pos{2} = [0.17 0.16 0.80 0.27];

    % raster plot of all normalized strengths on top
    subplot('Position',subplot_pos{1})
    plot([0 1],[-25 -25],'k:','LineWidth',1); hold on;
    plot([0 1],[-41 -41],'k:','LineWidth',1);
    
    spacing = 1;  
    for gp = 1:numel(group)
        plot(norm_nonstrongest{gp}(:,2),norm_nonstrongest{gp}(:,1)*-1-spacing,'o','MarkerSize',2,'Color',pColor{gp}); hold on;
        spacing = spacing + size(strongest{gp},1)+3;
    end
    
    ax = gca;
    ax.YTick = [-64, -39, -23,-2]; 
    ax.YTickLabel = {'57', '36', '23','1'}; 
    box on;
    ax.XTickLabel = [];
    ax.YAxis.TickLength = [0 0];
    
    xlim([0 1]); ylim([-66,0])
    set(gca,'FontSize',14,'LineWidth',1.5)
    ylabel('Neuron Number','FontSize',fontsz.label)

    % cumulative plot on bottom
    subplot('Position',subplot_pos{2})
    cdf_amp = cell(1,numel(group));
    cdf_p   = cell(1,numel(group));
    for gp = 1:numel(group)
        [cdf_p{gp}, cdf_amp{gp}] = ecdf(norm_nonstrongest{gp}(:,2));
        plot(cdf_amp{gp},cdf_p{gp},'LineWidth',2,'Color',pColor{gp}); hold on;
    end

    xlim([0 1]); ylim([0 1])
    set(gca,'FontSize',fontsz.tick,'LineWidth',1.5,'XTick',[0.2 0.4 0.6 0.8 1]);
    xlabel('Normalized Fiber Strength','FontSize',fontsz.label); 
    ylabel('Probability','FontSize',fontsz.label); 
    box off; hold off

    set(gcf,'PaperUnits','Centimeter','PaperPosition',[0 0 11 11])
    saveas(gcf,'Fig2FG_NormalizedFiberStrengths.svg')


%% fit GMM to peak amplitudes 
function GMM = fitGMM(peakamp)
    niter = 100;        % number of iterations
    numc = 20;          % max. number of components
   
    GMMtry = cell(niter,numc);
    BICtry = zeros(niter,numc);
    
    fitopts = statset('MaxIter',1000);
    h = waitbar(0,'Iteration...');
    
    figure; 
    for n = 1:1:niter
        rng('shuffle'); 
        for k = 1:1:numc 
            try
                GMMtry{n,k} = fitgmdist(peakamp,k,'Options',fitopts); % fit a k-component GMM
                BICtry(n,k) = GMMtry{n,k}.BIC;
            catch
                GMMtry{n,k} = [];
                BICtry(n,k) = 1e10;
            end
        end  

        plot(BICtry(n,1:k),'o'); hold on;
        xlabel('Number of GMM components'); ylabel('BIC')
        xlim([1,20])
        ylim([min(min(BICtry(1:n,1:k)))*0.99, min(min(BICtry(1:n,1:k)))*1.1])        
        waitbar(n/niter,h);
        pause(0.1)
    end
    close(h);
    
    % find the best model (most complex model within minimum BIC + lenience)
    GMM.BIC_lenience = 20;
    GMM.BICmin = min(min(BICtry));    
    GMM.BICdelta = BICtry - GMM.BICmin;
    GMM.BICdelta(GMM.BICdelta < GMM.BIC_lenience) = 0;              
    GMM.bestncomp_candidates = find(~all(GMM.BICdelta,1));
    GMM.bestncomp = max(GMM.bestncomp_candidates);   
    [GMM.idx(1), ~] = find(BICtry == min(BICtry(:,GMM.bestncomp)),1);
    GMM.idx(2) = GMM.bestncomp;
    GMM.GMMbest = GMMtry{GMM.idx(1), GMM.idx(2)};
    
    GMM.GMMbest_mu = sort(GMM.GMMbest.mu)';
    GMM.fiber_strengths = [GMM.GMMbest_mu(1), diff(GMM.GMMbest_mu)];
    GMM.mean_fiber_strength = mean(GMM.fiber_strengths);
    GMM.std_fiber_strength = std(GMM.fiber_strengths);
    GMM.cv_fiber_strength = GMM.std_fiber_strength/GMM.mean_fiber_strength;
end