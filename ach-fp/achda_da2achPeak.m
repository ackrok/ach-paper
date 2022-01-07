%% Align photometry data to acceeration
% good = [1:length(modAChDA)]; rmv = [3 5 18 22:24 27:28 30 33:34 36 41:42]; good(rmv) = []; % acceleration
good = [10:16,20,22:29,32:35,40:41,44:47]; % reward
beh = modAChDA(good);

%%
thres = 4; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
width = 5; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)

da2achpeak = cell(length(beh),3);
ach2achpeak = cell(length(beh),3);
lbl = {'immobility','locomotion','reward'};

h = waitbar(0, 'cross corr');
for x = 1:length(beh)
    
    %% acetylcholine
    fp_ach = beh(x).FP{1}; Fs = beh(x).Fs; % extract photometry signal from structure
    fp_mu = nanmean(fp_ach); % mean of entire photometry signal
    fp_ach = fp_ach - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
    idx_rew = extractEventST([1:length(fp_ach)]', floor(beh(x).reward), floor(beh(x).reward)+100, 1); % identify sample during reward
    idx_mov = extractEventST([1:length(fp_ach)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_imm = extractEventST([1:length(fp_ach)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
      
    fp_da = beh(x).FP{2}; % dopamine
    fp_da = fp_da - nanmean(fp_da);
 
    %% identify pauses
    % idxPauseMin = findIdxPause(fpInputVec, thres, width)
    peak_imm = findIdxPeak(fp_ach(idx_imm_nonRew), thres, width);
    peak_mov = findIdxPeak(fp_ach(idx_mov), thres, width);
    peak_rew = findIdxPeak(fp_ach(idx_rew), thres, width);
    
    %% extract ACh, DA signal during pauses
    % we will use idx_pause to center our window
    [da2achpeak{x,1}, sta_time] = getSTA(fp_da(idx_imm_nonRew), peak_imm/Fs, Fs, [-2 2]);
    da2achpeak{x,2} = getSTA(fp_da(idx_mov), peak_mov/Fs, Fs, [-2 2]);
    da2achpeak{x,3} = getSTA(fp_da(idx_rew), peak_rew/Fs, Fs, [-2 2]);
    
    ach2achpeak{x,1} = getSTA(fp_ach(idx_imm_nonRew), peak_imm/Fs, Fs, [-2 2]);
    ach2achpeak{x,2} = getSTA(fp_ach(idx_mov), peak_mov/Fs, Fs, [-2 2]);
    ach2achpeak{x,3} = getSTA(fp_ach(idx_rew), peak_rew/Fs, Fs, [-2 2]);
    
%%
    waitbar(x/length(beh),h);
end
close(h);
 
%% N = X mice: DA to ACh pause
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
da2peak_an = cell(3,1); ach2peak_an = cell(3,1);
max_val = []; max_time = [];
for y = 1:2
    for x = 1:nAn
        ii = find(strcmp(tmp,uni{x}));
        da2peak_an{y}(:,x) = nanmean([da2achpeak{ii,y}],2);
        ach2peak_an{y}(:,x) = nanmean([ach2achpeak{ii,y}],2);
    end
    [max_val(:,y), ii] = max(da2peak_an{y});
    max_time(:,y) = sta_time(ii);
end

%% PLOT: STATS -- DA to ACh pause -- comparing peak magnitude, peak latency
group = [1*ones(nAn,1);2*ones(nAn,1)];%;3*ones(nAn,1)];
[~,~,stats] = anova1(max_time(:),group,'off'); [c_time] = multcompare(stats,'Display','off');
[~,~,stats] = anova1(max_val(:),group,'off'); [c_val] = multcompare(stats,'Display','off');

fig = figure; fig.Position([3 4]) = [1000 900];
subplot(2,2,1);
violinplot(max_time); 
xticklabels({'immobility','locomotion','reward'});
ylabel('Latency to ACh pause (s)'); % ylim([-0.4 0]); yticks([-0.4:0.1:0]);
axis('square');
title(sprintf('DA to ACh pause - peak latency \n p-value: I/L - %1.3f',c_time(6))); 

subplot(2,2,2);
violinplot(max_val); 
xticklabels({'immobility','locomotion','reward'});
ylabel('rDA1m (%dF/F)'); %ylim([-0.4 0]); yticks([-0.4:0.1:0]);
axis('square');
title(sprintf('DA to ACh pause - peak magnitude \n p-value: I/L - %1.3f',c_time(6))); 

movegui(gcf,'center');

% PLOT: AVERAGES for all groups
clr = {'r','g','b'};

% fig = figure; fig.Position([3 4]) = [1000 420];
subplot(2,2,3);
plot([0 0],[-6 6],'--k'); hold on; 
for y = 1:2
shadederrbar(sta_time, nanmean(da2peak_an{y},2), SEM(da2peak_an{y},2), clr{y}); hold on
end
xlabel('Latency to ACh peak (s)'); ylabel('rDA1m (%dF/F)'); xlim([-1 1]);
axis('square')
title(sprintf('DA to ACh peak (n = %d mice)',nAn));

subplot(2,2,4);
plot([0 0],[-4 14],'--k'); hold on; 
for y = 1:2
shadederrbar(sta_time, nanmean(ach2peak_an{y},2), SEM(ach2peak_an{y},2), clr{y}); hold on
end
xlabel('Latency to ACh peak (s)'); ylabel('ACh3.0 (%dF/F)'); xlim([-1 1]);
axis('square')
title(sprintf('ACh to ACh peak (n = %d mice)',nAn));
movegui(gcf,'center');

%%
figure; 
subplot(1,2,1);
plot(max_time(:), max_val(:), '.m', 'MarkerSize',15);
xlim([-0.4 0]); xticks([-0.4:0.1:0])
ylim([0 14]);
axis('square')

subplot(1,2,2); hold on;
for x = 1:3; plot(max_time(:,x), max_val(:,x), clr{x}, 'LineStyle','none','Marker','.','MarkerSize',15); end
xlim([-0.4 0]); xticks([-0.4:0.1:0])
ylim([0 14]);
axis('square')

%% findIdxPeak
function [idxPeakMax, idxPeakStart, peakMag, crossGroups, crossStart] = findIdxPeak(fpInputVec, thres, width)
    % thresholds = [1:20]; % thresholds: defined as deviation from new baseline
    % peak_test = cell(1,length(thresholds));
    % mu = []; sigma = [];
    % for t = 1:length(thresholds) % Iterate over possible thresholds
    %     thres = thresholds(t);
    cross = find(fpInputVec > thres); % Find all samples below threshold
    [crossGroups, crossStart] = consecutive_vec2cell(cross); % Extract grouped crossings
    peakMag = nan(length(crossGroups),1); % Initialize vector for amplitudes
    idxPeakStart = []; idxPeakMax = [];
    for c = 1:length(crossGroups) % Iterate over discrete pauses
        fp_vec = [fpInputVec(crossGroups{c})]; % Vector for photometry values that cross during this pause/peak
        if length(crossGroups{c}) > width % Only find amplitudes for peaks that are X samples long (1 sample = 20ms)
            [peakMag(c),ii] = max(fp_vec); % Peaks amplitude
            idxPeakStart = [idxPeakStart; c]; % Index of crosses that satisfy parameters (width)
            idxPeakMax = [idxPeakMax; crossGroups{c}(ii)]; % Index of max point in cross
        end
    end
    %         mu(t) = nanmean(peak_test{t});
    %         sigma(t) = nanstd(peak_test{t});
    %     end
    %     figure; errorbar(mu, sigma); hold on; errorbar(0,nanmean(m),nanstd(m))
end