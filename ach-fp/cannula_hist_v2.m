raw = struct;
for y = 1:3
    %%
    fPath = 'C:\Users\Anya\Desktop\FP_LOCAL\4_glutamate\rawdata\'; 
    [fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
    if ~iscell(fName); fName = {fName}; end
    %%
    tmp = struct;
    for x = 1:length(fName) 
        load(fullfile(fPath,fName{x})); 
        [an,b] = strtok(fName{x},'_'); day = strtok(b,'_');
        tmp(x).rec = [an,'_',day]; 
        tmp(x).site = 'DLS';
        tmp(x).task = 'reward';
        tmp(x).Fs = 50;
        tmp(x).demod = data.final.nbFP{1}; %Store the demodulated signal in the non-baseline data
        tmp(x).baseline = (data.final.FP{1}./100).*data.final.FPbaseline{1}; % Baselined, but NOT dF/F
        tmp(x).fp = data.final.FP{1};
        tmp(x).reward = data.final.rew.onset;
        tmp(x).on = data.final.mov.onsets;
        tmp(x).off = data.final.mov.offsets;
        tmp(x).onRest = data.final.mov.onsetsRest;
        tmp(x).offRest = data.final.mov.offsetsRest;
        fprintf('Extracted from %s\n',fName{x});
    end
    raw(y).raw = tmp;
    %%
    demod = []; baseline = []; fp = [];
    for x = 1:4
        demod(:,x) = raw(y).raw(x).demod;
        baseline(:,x) = raw(y).raw(x).baseline; 
        fp(:,x) = raw(y).raw(x).fp; 
    end
    raw(y).time = data.final.time;
    raw(y).demod = demod; raw(y).baseline = baseline; raw(y).fp = fp;
end

raw(1).win = [20 40; 20 40; 20 40; 20 40]; raw(2).win = raw(1).win; 
raw(3).win = raw(1).win + 10; raw(3).win(4,:) = [20 40];

%%
fig = figure; fig.Position(3) = 1375;
for y = 1:3
    sp(y) = subplot(1,3,y); hold on;
    plot(raw(y).time/60, raw(y).demod);
    xlabel('Time (min)'); ylabel('demod V');
    axis('square');
    title(sprintf('%s demod',raw(y).inf));
end
linkaxes(sp,'y');
movegui(gcf,'center');

%% photometry signal V for histogram
fp_hist_V = cell(4,3);
for x = 1:4
    for y = 1:3
        fp = raw(y).raw(x).baseline;
        fp = fp - prctile(fp, 5); % subtract bottom 1% of points
        idx = [raw(y).win(x,1).*(50*60) : raw(y).win(x,2).*(50*60)]'; % infusion window
        idx_rew = extractEventST(idx,floor(raw(y).raw(x).reward),floor(raw(y).raw(x).reward)+100,1); % reward periods during infusion window
        idx = idx(~ismember(idx, idx_rew)); % remaining index are during infusion window and NOT during reward
        idx = extractEventST(idx,raw(y).raw(x).onRest,raw(y).raw(x).offRest,1); % rest periods during infusion window
        fp = fp(idx).*100;
        fp_hist_V{x,y} = fp;
    end  
end

%% PLOT: histogram
fig = figure; % fig.Position(3) = 1375;
clr = [0.05 0.75 0.45; 0.5 0.2 0.55; 0.85 0.35 0.1]; % teal, purple, orange
bin = 0.02;
for x = 1:4
    sp(x) = subplot(2,2,x); hold on
    for y = [1 3]
        histogram(fp_hist_V{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',clr(y,:),'FaceAlpha',0.2,'EdgeAlpha',0.1); 
    end  
    xlabel('photometry signal (V)'); xlim([-0.5 1.5]); legend({raw([1 3]).inf});
    % axis('square');
    title(sprintf('%s',strtok(raw(y).raw(x).rec,'_')));
end
linkaxes(sp,'x');
movegui(gcf,'center');

%% PLOT: histogram stats
group = [1*ones(4,1);2*ones(4,1);3*ones(4,1)];
xx = cellfun(@nanmean, fp_hist_V); [~,~,stats] = anova1(xx(:),group,'off'); [c_mean] = multcompare(stats,'Display','off');
xx = cellfun(@nanstd, fp_hist_V); [~,~,stats] = anova1(xx(:),group,'off'); [c_std] = multcompare(stats,'Display','off');

fig = figure; fig.Position(3) = 1375;
subplot(1,3,1);
violinplot(cellfun(@nanmean, fp_hist_V)); xticklabels({'acsf','d1d2-ant','nbqx'});
ylabel('mean'); axis('square');
title(sprintf('MEAN p-value: a/d - %1.3f | a/n - %1.3f | d/n - %1.3f',c_mean(1,6),c_mean(2,6),c_mean(3,6))); 
subplot(1,3,2);
violinplot(cellfun(@nanstd, fp_hist_V)); xticklabels({'acsf','d1d2-ant','nbqx'});
ylabel('stdev'); axis('square');
title(sprintf('STDEV p-value: a/d - %1.3f | a/n - %1.3f | d/n - %1.3f',c_std(1,6),c_std(2,6),c_std(3,6))); 

movegui(gcf,'center');

%%
x = 2; % AK 197
fig = figure; fig.Position(3) = 1375;
clr = [0.05 0.75 0.45; 0.5 0.2 0.55; 0.85 0.35 0.1]; % teal, purple, orange
subplot(1,3,1); hold on
histogram(fp_hist_V{x,1},'BinWidth',bin,'Normalization','probability','FaceColor',clr(1,:),'FaceAlpha',0.2,'EdgeAlpha',0.1); 
yyaxis right; histogram(fp_hist_V{x,2},'BinWidth',bin,'Normalization','probability','FaceColor',clr(2,:),'FaceAlpha',0.2,'EdgeAlpha',0.1); 
xlabel('photometry signal (V)'); xlim([-0.5 1.5]); legend({raw([1 2]).inf});
axis('square');
title(sprintf('%s',strtok(raw(y).raw(x).rec,'_')));

x = 4; % AK 197
subplot(1,3,2); hold on
histogram(fp_hist_V{x,1},'BinWidth',bin,'Normalization','probability','FaceColor',clr(1,:),'FaceAlpha',0.2,'EdgeAlpha',0.1); 
yyaxis right; histogram(fp_hist_V{x,3},'BinWidth',bin,'Normalization','probability','FaceColor',clr(3,:),'FaceAlpha',0.2,'EdgeAlpha',0.1); 
xlabel('photometry signal (V)'); xlim([-0.5 1.5]); legend({raw([1 3]).inf});
axis('square');
title(sprintf('%s',strtok(raw(y).raw(x).rec,'_')));

movegui(gcf,'center');


%% photometry signal dF/F for histogram
fp_hist_df = cell(4,3);
for x = 1:4
    for y = 1:3
        fp = raw(y).raw(x).fp;
        fp = fp - s(y).base_mu_rest(x);
        fp = fp - prctile(fp, 1); % subtract bottom 1% of points
        idx = [raw(y).win(x,1).*(50*60) : raw(y).win(x,2).*(50*60)]'; % infusion window
        idx_rew = extractEventST(idx,floor(raw(y).raw(x).reward),floor(raw(y).raw(x).reward)+100,1); % reward periods during infusion window
        idx = idx(~ismember(idx, idx_rew)); % remaining index are during infusion window and NOT during reward
        idx = extractEventST(idx,raw(y).raw(x).onRest,raw(y).raw(x).offRest,1); % rest periods during infusion window
        fp = fp(idx);
        fp_hist_df{x,y} = fp;
    end  
end
%% PLOT: histogram
fig = figure; % fig.Position(3) = 1375;
clr = [0.05 0.75 0.45; 0.5 0.2 0.55; 0.85 0.35 0.1]; % teal, purple, orange
bin = 0.2;
for x = 1:4
    sp(x) = subplot(2,2,x); hold on
    for y = [1 2]
        histogram(fp_hist_df{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',clr(y,:),'FaceAlpha',0.2,'EdgeAlpha',0.1); 
    end  
    xlabel('dF/F'); legend({raw([1 2]).inf});
    % axis('square');
    title(sprintf('%s',strtok(raw(y).raw(x).rec,'_')));
end
linkaxes(sp,'x');
movegui(gcf,'center');
%% PLOT: histogram stats
group = [1*ones(4,1);2*ones(4,1);3*ones(4,1)];
xx = cellfun(@nanmean, fp_hist_df); [~,~,stats] = anova1(xx(:),group,'off'); [c_mean] = multcompare(stats,'Display','off');
xx = cellfun(@nanstd, fp_hist_df); [~,~,stats] = anova1(xx(:),group,'off'); [c_std] = multcompare(stats,'Display','off');

fig = figure; fig.Position(3) = 1375;
subplot(1,3,1);
violinplot(cellfun(@nanmean, fp_hist_df)); xticklabels({'acsf','d1d2-ant','nbqx'});
ylabel('mean'); axis('square');
title(sprintf('MEAN p-value: a/d - %1.3f | a/n - %1.3f',c_mean(1,6),c_mean(2,6))); 
subplot(1,3,2);
violinplot(cellfun(@nanstd, fp_hist_df)); xticklabels({'acsf','d1d2-ant','nbqx'});
ylabel('stdev'); axis('square');
title(sprintf('STDEV p-value: a/d - %1.3f | a/n - %1.3f',c_std(1,6),c_std(2,6))); 

movegui(gcf,'center');
%%
x = 2; % AK 190
fig = figure; fig.Position(3) = 1375;
clr = [0.05 0.75 0.45; 0.5 0.2 0.55; 0.85 0.35 0.1]; % teal, purple, orange
subplot(1,3,1); hold on
histogram(fp_hist_df{x,1},'BinWidth',bin,'Normalization','probability','FaceColor',clr(1,:),'FaceAlpha',0.2,'EdgeAlpha',0.1); 
yyaxis right; histogram(fp_hist_df{x,2},'BinWidth',bin,'Normalization','probability','FaceColor',clr(2,:),'FaceAlpha',0.2,'EdgeAlpha',0.1); 
xlabel('ACh (%dF/F)'); xlim([-5 15]); legend({raw([1 2]).inf});
axis('square');
title(sprintf('%s',strtok(raw(y).raw(x).rec,'_')));

x = 2; % AK 190
subplot(1,3,2); hold on
histogram(fp_hist_df{x,1},'BinWidth',bin,'Normalization','probability','FaceColor',clr(1,:),'FaceAlpha',0.2,'EdgeAlpha',0.1); 
yyaxis right; histogram(fp_hist_df{x,3},'BinWidth',bin,'Normalization','probability','FaceColor',clr(3,:),'FaceAlpha',0.2,'EdgeAlpha',0.1); 
xlabel('ACh (%dF/F)'); xlim([-5 15]); legend({raw([1 3]).inf});
axis('square');
title(sprintf('%s',strtok(raw(y).raw(x).rec,'_')));

movegui(gcf,'center');