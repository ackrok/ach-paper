% a = [1:length(modACh)]; rmv = [32 38]; a(rmv) = [];
% beh = modACh(a); % Extract recordings with reward
% [align_full, time] = plot_fp2event(beh,[-4 2],0); % Align photometry to events

%%
%modACh DLS: thres = 2; width = 7
thres = 2; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
width = 7; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)

peak_onset = cell(1,length(beh));
peak_non = cell(1,length(beh));
peak_mu = []; peak_sigma = []; peak_num = []; peak_dur = {}; prc_captured = [];
non_freq = []; non_dur = {};
ach2achpeak = cell(length(beh),2);

h = waitbar(0, 'peaks - onset vs non');
for x = 1:length(beh); y = 1;

    fp = beh(x).FP{1}; Fs = beh(x).Fs; % extract photometry signal from structure
    fp_mu = nanmean(fp); % mean of entire photometry signal
    fp = fp - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
    idx_mov = extractEventST([1:length(fp)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_imm = extractEventST([1:length(fp)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    fp_imm = fp(idx_imm);

    %%
    a_on = align_full{x,y}([find(time == 0):find(time == 0.5)],:) - fp_mu; % SUBTRACT baseline%
    base = nanmean(align_full{x,y}([find(time == -4):find(time == -2)],:),1);
    a_on = a_on - base;
    a_on = a_on(:,all(~isnan(a_on))); % remove NaN trials
    [m, i] = max(a_on); i = i+find(time == 0)-1; % onset window: [0 1]s
    peak_onset{x} = m(:); % STORE
    
    idx_onset = extractEventST([1:length(fp)]', beh(x).on, beh(x).on+25, 1);
    fp_onset = fp(idx_onset); % photometry during onset window: [0 0.5]
    cross = find(fp_onset > thres); 
    [cross_cell, cross_pts] = consecutive_vec2cell(cross); % Extract grouped crossings
    peak_tmp = []; peak_tmp_idx = [];
    for c = 1:length(cross_cell)
        [peak_tmp(c),ii] = max(fp_onset(cross_cell{c})); 
        peak_tmp_idx(c) = cross_cell{c}(ii);
    end % Identify maximum of threshold crossing groups
    where = ismember(peak_tmp(:), m(:)+base); % Where is max onset
    peak_group = cross_cell(where); peak_pts = cross_pts(where,:); % Extract grouped crossings that match maximum onset group
    prc_captured(x) = length(find(where))/size(a_on,2); % Percent of onsets captured by threshold
    
    % peak_onset{x} = peak_tmp(where)'; % STORE only peaks that are above threshold
    
    %%
    % [idxPeakMax, idxPeakStart, peakMag, crossGroups, crossStart] = findIdxPeak(fpInputVec, thres, width)
    [idx_c_max, idx_c, peak_non{x}, ~, cross_pts] = findIdxPause(fp_imm, 'peak', thres, width); 

%%
    peak_mu(x,:) = [nanmean(peak_onset{x}), nanmean(peak_non{x})]; 
    peak_sigma(x,:) = [nanstd(peak_onset{x}), nanstd(peak_non{x})]; 
    peak_num(x,:) = [length(find(~isnan(peak_onset{x}))), length(find(~isnan(peak_non{x})))]; % Number of peaks satisfying conditions
    peak_dur{x} = peak_pts(:,2)-peak_pts(:,1); % Duration of onset peak that exceeds threshold
    non_freq(x) = (1/mean(diff(idx_c_max)))*Fs; % Frequency of maximum deflections
    non_dur{x} = cross_pts(idx_c,2)-cross_pts(idx_c,1); % Duration of pause/peak that exceeds threshold
    
%% ach to peaks
    idxpeak_imm = findIdxPause(fp(idx_imm), 'peak', thres, width);
    % idxpeak_mov = findIdxPause(fp(idx_onset), 'peak', thres, width);
    [ach2achpeak{x,1},sta_time] = getSTA(fp(idx_imm), idxpeak_imm/Fs, Fs, [-2 2]);
    % ach2achpeak{x,2} = getSTA(fp(idx_onset), idxpeak_mov/Fs, Fs, [-2 2]);
    ach2achpeak{x,2} = getSTA(fp, idx_onset(peak_tmp_idx)/Fs, Fs, [-2 2]);
%% 
waitbar(x/length(beh),h);
end
close(h);

%% PLOT: Pause Amplitude errorbars
% figure; hold on; 
% errorbar(peak_mu(:,1),peak_sigma(:,1),'g'); errorbar(peak_mu(:,2),peak_sigma(:,2),'k')
% ylabel('Peak Amplitude (%dF/F)'); legend({'onset','non-movement'});
% title(sprintf('Peak Amplitude (thres = %1.1f) (p = %1.2f)',thres,signrank(peak_mu(:,1), peak_mu(:,2))))

%% PLOT: HISTOGRAM
% figure;
% for x = 1:length(beh)
%     sp(x) = subplot(4,6,x); hold on
%     histogram(peak_onset{x},'BinWidth',0.5,'Normalization','probability','FaceColor','g','FaceAlpha',0.2,'EdgeAlpha',0.1);
%     if ~isempty(peak_non{x})
%     yyaxis right; histogram(peak_non{x},'BinWidth',0.5,'Normalization','probability','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.1);
%     end
%     title(sprintf('%s',beh(x).rec));
% end

%% n = X mice
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
onset_an = [];
mag_an_on = cell(nAn,1); mag_an_non = cell(nAn,1); % Initialize cell arrays for by animal grouping
dur_an_on = cell(nAn,1); dur_an_non = cell(nAn,1);
freq_an = [];
capture_an = [];
ach2peak_an = cell(2,1);

for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    ach2peak_an{1}(:,x) = nanmean([ach2achpeak{ii,1}],2);
    ach2peak_an{2}(:,x) = nanmean([ach2achpeak{ii,2}],2);
    onset_tmp = [];
    for y = 1:length(ii)
        mag_an_on{x} = [mag_an_on{x}; peak_onset{ii(y)}];
        mag_an_non{x} = [mag_an_non{x}; peak_non{ii(y)}];
        dur_an_on{x} = [dur_an_on{x}; peak_dur{ii(y)}];
        dur_an_non{x} = [dur_an_non{x}; non_dur{ii(y)}];
        onset_tmp = [onset_tmp, align_full{ii(y),1}];
    end
    onset_an(:,x) = nanmean(onset_tmp,2);
    freq_an(x) = nanmean(non_freq(ii));
    capture_an(x) = nanmean(prc_captured(ii));
end

mag_an_mu = [cellfun(@nanmean, mag_an_on), cellfun(@nanmean, mag_an_non)]; 
dur_an_mu = [cellfun(@nanmean, dur_an_on), cellfun(@nanmean, dur_an_non)]; 
dur_an_mu = dur_an_mu.*(1/Fs); % convert to seconds

%% PLOT: n = X mice -- MAGNTIDUE + DURATION errorbars
fig = figure; fig.Position(3) = 1375;

subplot(1,3,1); hold on; clr = {'k','g'};
plot([0 0],[-3 13],'--k'); plot([-1 1],[4 4],'--k');  
for y = 1:length(ach2peak_an)
shadederrbar(sta_time, nanmean(ach2peak_an{y},2), SEM(ach2peak_an{y},2), clr{y}); hold on
end
xlabel('Latency to ACh peak (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('ACh3.0 (%dF/F)'); yticks([-4:4:14]);
axis('square')
title(sprintf('ACh to ACh (n = %d) (every %1.2f s)',nAn,1/nanmean(freq_an)));

subplot(1,3,2); hold on
plot_me = mag_an_mu; % CHANGE
%
p = signrank(plot_me(:,1),plot_me(:,2));%,'tail','left'); % null hypothesis: mean difference = 0 VS difference < 0 (onset peak MAGNITUDE is more positive than rest)
errorbar(abs(nanmean(plot_me)'),abs(SEM(plot_me,1)'),'.k')
plot([1.15; 1.85].*ones(2,nAn),abs([plot_me']),'.:k');
ylabel('peak amplitude (%dF/F)'); ylim([0 10]); yticks([0:4:10]);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'onset','immobility'});
axis('square')
title(sprintf('Peak Amplitude (thres = %1.1f) (p = %1.2f)',thres,p));

subplot(1,3,3); hold on
plot_me = dur_an_mu; % CHANGE
%
p = signrank(plot_me(:,1),plot_me(:,2)); %,'tail','right'); % null hypothesis: mean difference = 0 VS difference < 0 (onset peak DURATION is more positive than rest)
errorbar(abs(nanmean(plot_me)'),abs(SEM(plot_me,1)'),'.k')
plot([1.15; 1.85].*ones(2,nAn),abs([plot_me']),'.:k');
ylabel('peak duration (s)'); ylim([0 0.5]); yticks([0:0.1:0.5]);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'onset','immobility'});
axis('square')
title(sprintf('Peak Duration (width = %1.1f) (p = %1.2f)',width,p));

movegui(gcf,'center');

%% PLOT: n = X mice -- HISTOGRAM
% figure;
% for x = 1:nAn
%     sp(x) = subplot(3,4,x); hold on
%     % histogram(peak_an_on{x},'BinWidth',0.5,'Normalization','probability','FaceColor','g','FaceAlpha',0.2,'EdgeAlpha',0.1);
%     % yyaxis right; histogram(peak_an_non{x},'BinWidth',0.5,'Normalization','probability','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.1);
%     errorbar(1, nanmean(peak_an_dur{x}), nanstd(peak_an_dur{x}), 'g')
%     errorbar(2, nanmean(non_an_dur{x}), nanstd(non_an_dur{x}), 'k'); xlim([0.5 2.5]); xticks([1 2]);
%     title(sprintf('%s',uni{x}));
% end
% linkaxes(sp,'y');
% 
%% PLOT: freq, duration
% figure; violinplot(non_dur); ylabel('(samples)'); title('Average Duration of Non-Onset Peaks');
% figure; violinplot(non_freq); ylabel('Freq (Hz)'); title('Frequency of Non-Onset Peaks');
% 