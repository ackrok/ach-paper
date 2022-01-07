% a = [10:16,20,22:30,36:37,39:43];
% beh = modACh(a); % Extract recordings with reward
% [align_full, time] = plot_fp2event(beh,[-6 2],0); % Align photometry to events

% a = [7:12,14:16];
% beh = allDMS(a); % Extract recordings with reward
% [align_full, t] = plot_fp2event(beh,[-6 2]); % Align photometry to events

%%
% modACh DLS: -3, 7
% allDMS: -1, 7
thres = -1; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
width = 7; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)

rew_valMax = cell(length(beh),1); rew_dur = cell(length(beh),1);
imm_valMax = cell(length(beh),1); imm_dur = cell(length(beh),1);
freq = [];
% rew_mu = []; rew_sigma = []; rew_num = [];
% nonrew_freq = []; nonrew_dur = {};
ach2achpause = cell(length(beh),2);

h = waitbar(0, 'pauses');
for x = 1:length(beh); y = 1;
 
    fp = beh(x).FP{1}; Fs = beh(x).Fs; % extract photometry signal from structure
    fp_mu = nanmean(fp); % mean of entire photometry signal
    fp = fp - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
    idx_rew = extractEventST([1:length(fp)]', floor(beh(x).reward), floor(beh(x).reward)+100, 1); % identify sample during reward
    idx_mov = extractEventST([1:length(fp)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_imm = extractEventST([1:length(fp)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    fp_imm = fp(idx_imm_nonRew);
 
    %%
%     a_rew = align_full{x,y}([find(time==0):find(time==1)],:) - fp_mu; % SUBTRACT baseline%
%     a_rew = a_rew(:,all(~isnan(a_rew))); % remove NaN (unrewarded) trials
%     [m, i] = min(a_rew); i = i+find(time==0)-1; % pause window: [0 1]s
%     rew_valMax{x} = m(:); % STORE
    
    fp_rew = fp(idx_rew); % photometry during onset window: [-0.2 1]
    [rew_idxMax, rew_idxGood, rew_valMaxAll, rew_crossGroups, rew_crossStart] = findIdxPause(fp_rew, 'pause', thres, width);
    rew_startGood = rew_crossStart(rew_idxGood,:); % Extract start/end times for good pauses
    rew_valGood = rew_valMaxAll(rew_idxGood);
    tmp = [1:101:length(fp_rew)]'; % Reward start index, in index values of fp_rew
    tmp = ceil(rew_startGood./101); % Divide by length of window to find which reward this pause is in
    tmp((diff(tmp,1,2) ~= 0),:) = nan; % NaN if pause start/end not within same reward window
    uni = unique(tmp); uni = uni(~isnan(uni)); % Unique reward windows that have crossings
    rew_valGoodUnique = nan(length(uni),1); 
    rew_startGoodUnique = nan(length(uni),2);
    rew_idxMaxUnique = nan(length(uni),1);
    for c = 1:length(uni)
        cc = find(tmp(:,1) == uni(c)); % idxGood that matches unique reward windows
        [rew_valGoodUnique(c),ii] = min(rew_valGood(cc)); % Keep minimum for multiple crosses in single reward period
        rew_startGoodUnique(c,:) = rew_startGood(cc(ii),:);
        rew_idxMaxUnique(c) = rew_idxMax(cc(ii));
    end
    rew_valMax{x} = rew_valGoodUnique;
    
%     pause_tmp = [];
%     for c = 1:length(rew_crossGroups)
%         pause_tmp(c) = min(fp_rew(rew_crossGroups{c})); % Identify minima of threshold crossing groups
%     end 
%     rew_idxMatch = ismember(pause_tmp(:), m(:)); % Where is max onset
    % pause_group = cross_cell(where); 
    % rew_crossStart = rew_crossStart(rew_idxMatch,:); % Extract grouped crossings that match maximum onset group
    
    %%
    [imm_idxMax, imm_idxGood, imm_valMax{x}, ~, imm_crossStart] = findIdxPause(fp_imm, 'pause', thres, width);

%%
    rew_dur{x} = diff(rew_startGoodUnique,1,2); % Duration of onset peak that exceeds threshold
    imm_dur{x} = diff(imm_crossStart(imm_idxGood,:),1,2); % Duration of pause/peak that exceeds threshold
    freq(x) = (1/mean(diff(imm_idxMax)))*Fs; % Frequency of maximum deflections
    
    % rew_mu(x,:) = [nanmean(rew_valMax{x}), nanmean(imm_valMax{x})]; % MEAN of pause max amplitude
    % rew_sigma(x,:) = [nanstd(rew_valMax{x}), nanstd(imm_valMax{x})]; % STD of pause max amplitude
    % rew_num(x,:) = [length(find(~isnan(rew_valMax{x}))), length(imm_idxGood)]; % Number of pauses satisfying conditions
    % nonrew_freq(x) = (1/mean(diff(imm_idxMax)))*Fs; % Frequency of maximum deflections
    
    %% ach to pause
    % idxpause_imm = findIdxPause(fp(idx_imm_nonRew),'pause',thres,width);
    % idxpause_rew = findIdxPause(fp(idx_rew),'pause',thres,width);
    [ach2achpause{x,1},sta_time] = getSTA(fp(idx_imm_nonRew), imm_idxMax/Fs, Fs, [-2 2]);
    ach2achpause{x,2} = getSTA(fp(idx_rew), rew_idxMaxUnique/Fs, Fs, [-2 2]);
    
    waitbar(x/length(beh),h);
%%
end
close(h);
 
%% PLOT: Pause Amplitude errorbars
% figure; hold on; 
% errorbar(rew_mu(:,1),rew_sigma(:,1),'c'); errorbar(rew_mu(:,2),rew_sigma(:,2),'k')
% ylabel('Pause Amplitude (%dF/F)'); legend({'rew','non-rew'});
% title(sprintf('Pause Amplitude (thres = %1.1f) (p = %1.2f)',thres,signrank(rew_mu(:,1), rew_mu(:,2))))
 
%% PLOT: HISTOGRAM
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
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
mag_an_rew = cell(nAn,1); mag_an_nonrew = cell(nAn,1); % Initialize cell arrays for by animal grouping
dur_an_rew = cell(nAn,1); dur_an_nonrew = cell(nAn,1);
ach2pause_an = cell(2,1);
freq_an = [];
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    ach2pause_an{1}(:,x) = nanmean([ach2achpause{ii,1}],2);
    ach2pause_an{2}(:,x) = nanmean([ach2achpause{ii,2}],2);
    for y = 1:length(ii)
        mag_an_rew{x} = [mag_an_rew{x}; rew_valMax{ii(y)}];
        mag_an_nonrew{x} = [mag_an_nonrew{x}; imm_valMax{ii(y)}];
        dur_an_rew{x} = [dur_an_rew{x}; rew_dur{ii(y)}];
        dur_an_nonrew{x} = [dur_an_nonrew{x}; imm_dur{ii(y)}];
    end
    freq_an(x) = nanmean(freq(ii));
end

mag_an_mu = [cellfun(@nanmean, mag_an_rew), cellfun(@nanmean, mag_an_nonrew)]; 
dur_an_mu = [cellfun(@nanmean, dur_an_rew), cellfun(@nanmean, dur_an_nonrew)]; 
dur_an_mu = dur_an_mu.*(1/Fs); % convert to seconds

%% PLOT: n = X mice -- Pause MAGNTIDUE + DURATION errorbars
fig = figure; fig.Position(3) = 1375;

subplot(1,3,1); hold on; clr = {'k','b'};
plot([0 0],[-6 3],'--k'); plot([-1 1],[-2 -2],'--k');  
for y = 1:length(ach2pause_an)
shadederrbar(sta_time, nanmean(ach2pause_an{y},2), SEM(ach2pause_an{y},2), clr{y}); hold on
end
xlabel('Latency to ACh pause (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('ACh3.0 (%dF/F)'); yticks([-6:2:3]);
axis('square')
title(sprintf('ACh to ACh pause (n = %d mice)',nAn));

subplot(1,3,2); hold on
%
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
% PLOT: n = X mice -- Pause Duration errorbars
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
% figure;
% for x = 1:nAn
%     sp(x) = subplot(3,4,x); hold on
%     % histogram(pause_an_rew{x},'BinWidth',0.5,'Normalization','probability','FaceColor','b','FaceAlpha',0.2,'EdgeAlpha',0.1);
%     % yyaxis right; histogram(pause_an_nonrew{x},'BinWidth',0.5,'Normalization','probability','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.1);
%     violinplot(dur_an_nonrew{x});
%     title(sprintf('%s',uni{x}));
% end
% linkaxes(sp,'y');

%% PLOT: ACh to ACh pause
% figure; hold on; clr = {'k','b'};
% plot([0 0],[-6 3],'--k'); plot([-1 1],[-2 -2],'--k');  
% for y = 1:length(ach2pause_an)
% shadederrbar(sta_time, nanmean(ach2pause_an{y},2), SEM(ach2pause_an{y},2), clr{y}); hold on
% end
% xlabel('Latency to ACh pause (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
% ylabel('ACh3.0 (%dF/F)'); yticks([-6:2:3]);
% axis('square')
% title(sprintf('ACh to ACh pause (n = %d mice)',nAn));
