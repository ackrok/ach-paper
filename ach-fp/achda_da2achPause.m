%% Align photometry data to acceeration
% good = [1:length(modAChDA)]; rmv = [3 5 18 22:24 27:28 30 33:34 36 41:42]; good(rmv) = []; % acceleration
good = [10:16,20,22:29,32:35,40:41,44:47]; % reward
beh = modAChDA(good);

%%
thres = -2; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
width = 8; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)

da2achpause = cell(length(beh),3);
ach2achpause = cell(length(beh),3);
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
    pause_imm = findIdxPause(fp_ach(idx_imm_nonRew), thres, width);
    pause_mov = findIdxPause(fp_ach(idx_mov), thres, width);
    pause_rew = findIdxPause(fp_ach(idx_rew), thres, width);
    
    %% extract ACh, DA signal during pauses
    % we will use idx_pause to center our window
    [da2achpause{x,1}, sta_time] = getSTA(fp_da(idx_imm_nonRew), pause_imm/Fs, Fs, [-2 2]);
    da2achpause{x,2} = getSTA(fp_da(idx_mov), pause_mov/Fs, Fs, [-2 2]);
    da2achpause{x,3} = getSTA(fp_da(idx_rew), pause_rew/Fs, Fs, [-2 2]);
    
    ach2achpause{x,1} = getSTA(fp_ach(idx_imm_nonRew), pause_imm/Fs, Fs, [-2 2]);
    ach2achpause{x,2} = getSTA(fp_ach(idx_mov), pause_mov/Fs, Fs, [-2 2]);
    ach2achpause{x,3} = getSTA(fp_ach(idx_rew), pause_rew/Fs, Fs, [-2 2]);
    
%%
    waitbar(x/length(beh),h);
end
close(h);
 
% N = X mice: DA to ACh pause
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
da2pause_an = cell(3,1); ach2pause_an = cell(3,1);
max_val = []; max_time = [];
for y = 1:3
    for x = 1:nAn
        ii = find(strcmp(tmp,uni{x}));
        da2pause_an{y}(:,x) = nanmean([da2achpause{ii,y}],2);
        ach2pause_an{y}(:,x) = nanmean([ach2achpause{ii,y}],2);
    end
    [max_val(:,y), ii] = max(da2pause_an{y});
    max_time(:,y) = sta_time(ii);
end

%% PLOT: STATS -- DA to ACh pause -- comparing peak magnitude, peak latency
group = [1*ones(nAn,1);2*ones(nAn,1);3*ones(nAn,1)];
[~,~,stats] = anova1(max_time(:),group,'off'); [c_time] = multcompare(stats,'Display','off');
[~,~,stats] = anova1(max_val(:),group,'off'); [c_val] = multcompare(stats,'Display','off');

fig = figure; fig.Position([3 4]) = [1000 900];
subplot(2,2,1);
violinplot(max_time); 
xticklabels({'immobility','locomotion','reward'});
ylabel('Latency to ACh pause (s)'); ylim([-0.4 0]); yticks([-0.4:0.1:0]);
axis('square');
title(sprintf('DA to ACh pause - peak latency \n p-value: I/L - %1.3f | I/R - %1.3f | L/R - %1.3f',c_time(1,6),c_time(2,6),c_time(3,6))); 

subplot(2,2,2);
violinplot(max_val); 
xticklabels({'immobility','locomotion','reward'});
ylabel('rDA1m (%dF/F)'); %ylim([-0.4 0]); yticks([-0.4:0.1:0]);
axis('square');
title(sprintf('DA to ACh pause - peak magnitude \n p-value: I/L - %1.3f | I/R - %1.3f | L/R - %1.3f',c_val(1,6),c_val(2,6),c_val(3,6))); 

movegui(gcf,'center');

% PLOT: AVERAGES for all groups
clr = {'r','g','b'};

% fig = figure; fig.Position([3 4]) = [1000 420];
subplot(2,2,3);
plot([0 0],[-2 8],'--k'); hold on; 
for y = 1:3
shadederrbar(sta_time, nanmean(da2pause_an{y},2), SEM(da2pause_an{y},2), clr{y}); hold on
end
xlabel('Latency to ACh pause (s)'); ylabel('rDA1m (%dF/F)'); xlim([-1 1]);
axis('square')
title(sprintf('DA to ACh pause (n = %d mice)',nAn));

subplot(2,2,4);
plot([0 0],[-6 3],'--k'); hold on; 
for y = 1:3
shadederrbar(sta_time, nanmean(ach2pause_an{y},2), SEM(ach2pause_an{y},2), clr{y}); hold on
end
xlabel('Latency to ACh pause (s)'); ylabel('ACh3.0 (%dF/F)'); xlim([-1 1]);
axis('square')
title(sprintf('ACh to ACh pause (n = %d mice)',nAn));
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

%% findIdxPause
function idxPauseMin = findIdxPause(fpInputVec, thres, width)
    cross = find(fpInputVec < thres); % Find all samples below threshold
    [cross_cell, ~] = consecutive_vec2cell(cross); % Extract grouped crossings
    pause_nonRew = nan(length(cross_cell),1); % Initialize vector for pause amplitudes
    idx_c = []; idxPauseMin = [];
    for c = 1:length(cross_cell) % Iterate over discrete pauses
        fp_vec = [fpInputVec(cross_cell{c})]; % Vector for photometry values that cross during this pause
        if length(cross_cell{c}) > width % Only find amplitudes for pauses that are X samples long (1 sample = 20ms)
            [pause_nonRew(c),ii] = min(fp_vec); % Pause amplitude
            idx_c = [idx_c; c]; % Index of crosses that satisfy parameters (width)
            idxPauseMin = [idxPauseMin; cross_cell{c}(ii)]; % Index of minimum point in cross
        end
    end
end