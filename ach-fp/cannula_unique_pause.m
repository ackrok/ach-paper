%%
thres = -2; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
width = 8; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)

h = waitbar(0, 'pauses');
for y = 1:length(s) % iterate over infusion
    %%
    beh = s(y).s;
    % [align_full, time] = plot_fp2event(beh,[-6 2],0);
    % mag_mu = []; mag_sigma = []; num_pause = []; 
    mag_pause = cell(length(beh),2); 
    dur_pause = cell(length(beh),2); freq_pause = [];
    ach2achpause = cell(length(beh),2);

    for x = 1:length(beh) % iterate over animal
        fp = beh(x).FP{1}; Fs = 50; % extract full photometry signal
        idx_inf = [s(y).win(x,1).*(50*60) : s(y).win(x,2).*(50*60)]'; % infusion window
        rewWindow = 50; % how many samples after reward delivery is the reward window
        idx_inf_rew = extractEventST(idx_inf, floor(beh(x).reward), floor(beh(x).reward)+rewWindow, 1); % identify sample during reward
        idx_inf_rest = extractEventST(idx_inf, beh(x).onRest, beh(x).offRest, 1); % index: infusion + rest
        idx_inf_rest = idx_inf_rest(~ismember(idx_inf_rest, idx_inf_rew)); % exclude reward, include rest
        % fp = fp - s(y).base_mu_rest(x); % subtract mean of baseline recording to center around 0%
        fp = fp - nanmean(fp);
        
    %% pause identification
        [idxpause_imm, imm_good, imm_valmag, ~, imm_crossstart] = findIdxPause(fp(idx_inf_rest),'pause',thres,width);
        [idxpause_rew, rew_good, rew_valmag, ~, rew_crossstart] = findIdxPause(fp(idx_inf_rew),'pause',thres,width);
        
    %% ensure reward pause is minimum within reward window
        rew_startGood = rew_crossstart(rew_good,:); % Extract start/end times for good pauses
        rew_valGood = rew_valmag(rew_good);
        tmp = [1:rewWindow+1:length(idx_inf_rew)]'; % Reward start index, in index values of fp_rew
        tmp = ceil(rew_startGood./(rewWindow+1)); % Divide by length of window to find which reward this pause is in
        tmp((diff(tmp,1,2) ~= 0),:) = nan; % NaN if pause start/end not within same reward window
        uni = unique(tmp); uni = uni(~isnan(uni)); % Unique reward windows that have crossings
        rew_valGoodUnique = nan(length(uni),1); 
        rew_startGoodUnique = nan(length(uni),2);
        rew_idxMaxUnique = nan(length(uni),1);
        for c = 1:length(uni)
            cc = find(tmp(:,1) == uni(c)); % idxGood that matches unique reward windows
            [rew_valGoodUnique(c),ii] = min(rew_valGood(cc)); % Keep minimum for multiple crosses in single reward period
            rew_startGoodUnique(c,:) = rew_startGood(cc(ii),:);
            rew_idxMaxUnique(c) = idxpause_rew(cc(ii));
        end
        
    %% ensure imm pauses are also separated by window
%         tmp = find(diff(idxpause_imm) < 50); % Find immobility pauses that are within 1s
%         imm_good(tmp) = nan;
        
    %% ACh to ACh pause
        [ach2achpause{x,1},sta_time] = getSTA(fp(idx_inf_rest), idxpause_imm/Fs, Fs, [-2 2]);
        ach2achpause{x,2} = getSTA(fp(idx_inf_rew), rew_idxMaxUnique/Fs, Fs, [-2 2]);
        
    %%
        % mag_mu(x,:) = [nanmean(imm_valmag), nanmean(rew_valmag)];  % MEAN of pause max amplitude
        % mag_sigma(x,:) = [nanstd(imm_valmag), nanstd(rew_valmag)];  % STD of pause max amplitude
        % num_pause(x,:) = [length(imm_good), length(rew_good)]; % Number of pauses satisfying conditions
        mag_pause{x,1} = imm_valmag(imm_good); % Pause magnitude
        mag_pause{x,2} = rew_valGoodUnique; % Reward magnitude
        dur_pause{x,1} = diff(imm_crossstart(imm_good,:),1,2); % Duration of pause that exceeds threshold
        dur_pause{x,2} = diff(rew_startGoodUnique,1,2); % Duration of pause that exceeds threshold
        freq_pause(x) = (1/mean(diff(idxpause_imm)))*Fs; % Frequency of maximum deflections

    waitbar(x/length(beh),h);
%%
    end
    s(y).mag_pause = mag_pause; s(y).dur_pause = dur_pause; s(y).ach2achpause = ach2achpause;
    s(y).freq_pause = freq_pause;
end
close(h);

%% PLOT: Pause Amplitude errorbars
% figure; hold on; 
% errorbar(mag_mu(:,1),mag_sigma(:,1),'c'); errorbar(mag_mu(:,2),mag_sigma(:,2),'k')
% ylabel('Pause Amplitude (%dF/F)'); legend({'rew','non-rew'});
% title(sprintf('Pause Amplitude (thres = %1.1f) (p = %1.2f)',thres,signrank(mag_mu(:,1), mag_mu(:,2))))
%  
% %% PLOT: HISTOGRAM
% figure;
% for x = 1:length(beh)
%     sp(x) = subplot(4,6,x); hold on
%     %histogram(pause_rew{x},'BinWidth',0.5,'Normalization','probability','FaceColor','b','FaceAlpha',0.2,'EdgeAlpha',0.1);
%     %if ~isempty(pause_nonRew{x})
%     %yyaxis right; histogram(pause_nonRew{x},'BinWidth',0.5,'Normalization','probability','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.1);
%     %end
%     violinplot(nonrew_dur{x});
%     title(sprintf('%s',beh(x).rec));
% end
 
%% n = X mice
y = 2; beh = s(y).s;
ach2achpause = s(y).ach2achpause;
mag_pause = s(y).mag_pause;
dur_pause = s(y).dur_pause;

tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
mag_an_rew = cell(nAn,1); mag_an_nonrew = cell(nAn,1); % Initialize cell arrays for by animal grouping
dur_an_rew = cell(nAn,1); dur_an_nonrew = cell(nAn,1);
ach2pause_an = cell(2,1);
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    ach2pause_an{1}(:,x) = nanmean([ach2achpause{ii,1}],2);
    ach2pause_an{2}(:,x) = nanmean([ach2achpause{ii,2}],2);
    for jj = 1:length(ii)
        mag_an_rew{x} = [mag_an_rew{x}; mag_pause{ii(jj),2}];
        mag_an_nonrew{x} = [mag_an_nonrew{x}; mag_pause{ii(jj),1}];
        dur_an_rew{x} = [dur_an_rew{x}; dur_pause{ii(jj),2}];
        dur_an_nonrew{x} = [dur_an_nonrew{x}; dur_pause{ii(jj),1}];
    end
end

mag_an_mu = [cellfun(@nanmean, mag_an_rew), cellfun(@nanmean, mag_an_nonrew)]; 
dur_an_mu = [cellfun(@nanmean, dur_an_rew), cellfun(@nanmean, dur_an_nonrew)]; 
dur_an_mu = dur_an_mu.*(1/Fs); % convert to seconds

%% PLOT: n = X mice -- Pause MAGNTIDUE + DURATION errorbars
fig = figure; fig.Position(3) = 1375;

subplot(1,3,1); hold on; clr = {'k','b'};
plot([0 0],[-6 3],'--k'); plot([-1 1],[-2 -2],'--k');  
for x = 1:2
shadederrbar(sta_time, nanmean(ach2pause_an{x},2), SEM(ach2pause_an{x},2), clr{x}); hold on
end
xlabel('Latency to ACh pause (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('ACh3.0 (%dF/F)'); yticks([-6:2:3]);
axis('square')
title(sprintf('%s: ACh to ACh pause (n = %d)',nAn));

subplot(1,3,2); hold on
plot_me = mag_an_mu; % CHANGE
%
p = signrank(plot_me(:,1),plot_me(:,2));%,'tail','left'); % null hypothesis: mean difference = 0 VS difference < 0 (reward pause MAGNITUDE is more positive than non reward pause)
errorbar(abs(nanmean(plot_me)'),abs(SEM(plot_me,1)'),'.k')
plot([1.15; 1.85].*ones(2,nAn),abs([plot_me']),'.:k');
ylabel('magnitude (%dF/F)'); ylim([2 7]); yticks([3:6]);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'rew','immobility'});
axis('square')
title(sprintf('Pause Magnitude (thres = %d) (p = %1.2f)',thres,p));

subplot(1,3,3); hold on
plot_me = dur_an_mu; % CHANGE
%
p = signrank(plot_me(:,1),plot_me(:,2)); %,'tail','right'); % null hypothesis: mean difference = 0 VS difference < 0 (reward pause DURATION is more positive than non reward pause)
errorbar(abs(nanmean(plot_me)'),abs(SEM(plot_me,1)'),'.k')
plot([1.15; 1.85].*ones(2,nAn),abs([plot_me']),'.:k');
ylabel('duration (s)'); ylim([0 0.4]); yticks([0.1:0.1:0.3]);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'rew','immobility'});
axis('square')
title(sprintf('Pause Duration (width = %d) (p = %1.2f)',width,p));

movegui(gcf,'center');

%% PLOT: n = X mice -- HISTOGRAM
figure;
for x = 1:nAn
    sp(x) = subplot(3,4,x); hold on
    % histogram(pause_an_rew{x},'BinWidth',0.5,'Normalization','probability','FaceColor','b','FaceAlpha',0.2,'EdgeAlpha',0.1);
    % yyaxis right; histogram(pause_an_nonrew{x},'BinWidth',0.5,'Normalization','probability','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.1);
    violinplot(dur_an_nonrew{x});
    title(sprintf('%s',uni{x}));
end
linkaxes(sp,'y');

%% PLOT 1: ACh to ACh pause
fig = figure; fig.Position(3) = 1375;
clr = [0.05 0.75 0.45; 0.5 0.2 0.55; 0.85 0.35 0.1];
for y = 1:length(s)
    sp(y) = subplot(1,3,y); hold on
    mat = cell(2,1);
    for x = 1:4
        mat(:,x) = nanmean(s(y).ach2achpause{x,1},2);
        % plot(sta_time, nanmean(s(y).ach2achpause{x,1},2), 'Color', clr(y,:)); 
        % shadederrbar(sta_time, nanmean(s(y).ach2achpause{x,2},2), SEM(s(y).ach2achpause{x,2},2), clr(y,:)); hold on
    end
    xlabel('Latency to ACh pause (s)'); ylabel('ACh (%dF/F)'); axis('square')
    title(sprintf('%s',s(y).inf));
end
linkaxes(sp,'y'); linkaxes(sp,'x'); xlim([-1 1]);
movegui(gcf,'center');

%% PLOT 2: pause MAGNITUDE
fig = figure; fig.Position(3) = 1375;
for y = 1:length(s)
    sp(y) = subplot(1,3,y); hold on
    plotme = cellfun(@nanmean, s(y).mag_pause);
    violinplot(plotme);
    xticklabels({'imm','rew'}); ylabel('Pause Magnitude (%dF/F)'); axis('square')
    title(sprintf('%s (p = %1.3f)',s(y).inf,signrank(plotme(:,1),plotme(:,2))));
end
linkaxes(sp,'y');
movegui(gcf,'center');

%% PLOT 2: pause DURATION
fig = figure; fig.Position(3) = 1375;
for y = 1:length(s)
    sp(y) = subplot(1,3,y); hold on
    plotme = cellfun(@nanmean, s(y).dur_pause);
    violinplot(plotme.*(1000/Fs));
    xticklabels({'imm','rew'}); ylabel('Pause Duration (ms)'); axis('square')
    title(sprintf('%s (p = %1.3f)',s(y).inf,signrank(plotme(:,1),plotme(:,2))));
end
linkaxes(sp,'y');
movegui(gcf,'center');