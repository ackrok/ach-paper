%% create structure with infusion photometry data
s(1).s = acsf; s(2).s = d1d2; s(3).s = glu;
s(1).inf = 'aCSF'; s(2).inf = 'D1R/D2R'; s(3).inf = 'NMDA/AMPA';
s(1).win = [20 40; 20 40; 20 40; 20 40]; s(2).win = s(1).win; 
s(3).win = s(1).win + 10; s(3).win(4,:) = [20 40];
 
%% mu, std of rest ACh distribution during infusions
mu = []; sigma = []; % initialize/clear matrix
fp_cell = {}; fp_base_cell = {};
% fp_1 = {}; fp_2 = {};
z = 1; % CHANGE ACH or DA
 
for y = 1:length(s) % iterate over infusion
    beh = s(y).s;
    for x = 1:length(beh) % iterate over animal
        %%
        ii = find(strcmp({cannula.rec},beh(x).rec)); ii = ii(1); % identify matching baseline recording from cannula larger structure
        fp_base = cannula(ii).FP{z}; % extract full baseline photometry signal
        idx_base_rest = extractEventST([1:length(fp_base)]', cannula(ii).onRest, cannula(ii).offRest, 1);
        fp_base_rest_mu = nanmean(fp_base(idx_base_rest)); % mean of baseline photometry recording
        fp_base_rest_std = nanstd(fp_base(idx_base_rest));
        
        fp_base_z = (fp_base(idx_base_rest) - fp_base_rest_mu)/fp_base_rest_std;
        fp_base = fp_base(idx_base_rest) - fp_base_rest_mu; % subtract mean of baseline recording to center around 0%
        
        %%
        fp = beh(x).FP{z}; % extract full photometry signal
        % fp = fp - fp_base_rest_mu; % subtract mean of baseline recording to center around 0%
        fp = (fp - fp_base_rest_mu)/fp_base_rest_std; % Z-SCORE
        idx_inf = [s(y).win(x,1).*(50*60) : s(y).win(x,2).*(50*60)]'; % infusion window
        idx_inf_rest = extractEventST(idx_inf, beh(x).onRest, beh(x).offRest, 1); % index: infusion + rest
        fp_inf_rest = fp(idx_inf_rest); % extract photometry during infusion window
        mu_inf_rest = nanmean(fp_inf_rest); % MU of fp rest/inf values
        sigma_inf_rest = nanstd(fp_inf_rest); % STD of fp rest/inf values
        %%
        mu(x,y) = mu_inf_rest; % STORE
        sigma(x,y) = sigma_inf_rest; 
        fp_cell{x,y} = fp_inf_rest; fp_base_cell{x,y} = fp_base_z;
        
        % fp_1{x,y} = beh(x).FP{z}(extractEventST([1:length(fp)]', beh(x).onRest, beh(x).offRest, 1)); % rest baseline
        % fp_2{x,y} = cannula(ii).FP{z}(extractEventST([1:length(fp_base)]', cannula(ii).onRest, cannula(ii).offRest, 1)); % rest baseline
        % fp_1{x,y} = fp_1{x,y} - fp_base_rest_mu; fp_2{x,y} = fp_2{x,y} - fp_base_rest_mu; % subtract baseline
    end
end

%% PLOT violinplot STD
figure;
violinplot(sigma); xticklabels({s.inf});
title('STD of ACh rest distribution (%dF/F)'); 
ylabel('standard deviation (%dF/F)'); ylim([0 2.5]);
 
%% PLOT histograms
figure;
bin = 0.2;
for x = 1:size(fp_cell,1); for y = 1:size(fp_cell,2)
    z = x + (4*(y-1)); sp(z) = subplot(3,4,z); hold on
    % histogram(fp_cell{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',[0.05 0.75 0.45],'FaceAlpha',0.2,'EdgeAlpha',0.1);
    % histogram(fp_2{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeAlpha',0.1);
    % histogram(fp_1{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',[0.05 0.75 0.45],'FaceAlpha',0.2,'EdgeAlpha',0.1);
    title(sprintf('%s - %s',strtok(s(y).s(x).rec,'_'),s(y).inf));
    end; end
linkaxes(sp,'x');

%% EXAMPLES: AK190
x = 2; % AK190
clr = [0.05 0.75 0.45; 0.5 0.2 0.55; 0.85 0.35 0.1]; % teal, purple, orange
fig = figure; fig.Position(3) = 1375;
for y = 1:3
    sp(y) = subplot(1,3,y); hold on; bin = 0.2;
    histogram(fp_base_cell{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeAlpha',0.1);
    yyaxis right; histogram(fp_cell{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',clr(y,:),'FaceAlpha',0.2,'EdgeAlpha',0.1); 
    % histogram(fp_2{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeAlpha',0.1);
    % yyaxis right; histogram(fp_1{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',clr(y,:),'FaceAlpha',0.2,'EdgeAlpha',0.1); 
    legend({'pre-infusion',s(y).inf});
    xlabel('rDA1m fluorescence (dF/F %)'); xlim([-10 15]);
    axis('square');
    title(sprintf('%s - %s',strtok(s(y).s(x).rec,'_'),s(y).inf));
end; 
linkaxes(sp,'x');
movegui(gcf,'center');