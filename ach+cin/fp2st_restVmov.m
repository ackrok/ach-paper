%% (1) Load data
load('fullcin_Feb21_ACh+DA+PF.mat')
load('beh_wt_ACh+DA+PF.mat')
sub = cinwt([1:19, 30:36, 63:73]);  % Extract CINs from recordings with ACh
beh = behwt([1:19]); % Extract ACh recordings

%%
velThres = 0.25;
Fs = 50;
minRestTime = 4*Fs; minRunTime = 1*Fs;
timeThres = 4*Fs; timeShift = 0.5*Fs;

for idx_b = 1:length(beh)
    vel = beh(idx_b).vel; vel = abs(vel);
    [onsetInd,offsetInd] = getOnsetOffset(-vel,-velThres,minRunTime,minRestTime,1);
    [onsetInd,offsetInd] = adjOnsetOffset(onsetInd,offsetInd,timeThres,vel);
    onsetInd = onsetInd+timeShift; offsetInd = offsetInd-timeShift;
    beh(idx_b).onRest = onsetInd; beh(idx_b).offRest = offsetInd;
end

%% (2) STA
Fs = 50; time = [-1:1/Fs:1]; nShuff = 50; % CHANGE: window for STA
mat = struct; % Initialize structure
h = waitbar(0, 'STA: signal to CIN st');
for x = 1:length(sub)
    idx = find(strcmp({beh.rec},sub(x).rec));
    if isempty(idx) || isempty(sub(x).burst); continue; end

    %% Select event times and signal to align to
    %st = sub(x).st; % Aligning to spike times
    st = sub(x).burst.Windows(:,1); % Aligning to burst onset times
    
    %sig = [beh(idx).vel(1); diff(movmean(beh(idx).vel,10))]; % Aligning acceleration signal
    sig = beh(idx).FP{1}; % Aliging photometry signal
    
    %% Extract event times during rest & movement
    st_sub = cell(1,2); st_sub_forShuff = cell(1,2);
    st_sub{1} = extractEventST(st, beh(idx).on/Fs, beh(idx).off/Fs, 0); % Event times during movement
    st_sub{2} = extractEventST(st, beh(idx).onRest/Fs, beh(idx).offRest/Fs, 0); % Event times during rest
    st_sub_forShuff{1} = extractEventST(st, beh(idx).on/Fs, beh(idx).off/Fs, 1); % Event times during movement
    st_sub_forShuff{2} = extractEventST(st, beh(idx).onRest/Fs, beh(idx).offRest/Fs, 1); % Event times during rest
    
    %% 
    sta_raw = cell(length(st_sub),1); sta_z = sta_raw; sta_shuff = sta_raw;
    for y = 1:length(st_sub) % Iterate over rest and movement
        if isempty(st_sub{y}); continue; end
        if length(st_sub_forShuff{1}) < 2; continue; end
        sta_raw{y} = getSTA(sig, st_sub{y}, Fs, [time(1), time(end)]);
        st_new = shuffleST(st_sub_forShuff{1}, nShuff); % Shuffle event times n times
        mat_shuff = [];
        for z = 1:nShuff
            tmp_shuff = getSTA(sig, st_new{z}, Fs, [time(1), time(end)]); % Aligning signal to shuffled event times
            mat_shuff = [mat_shuff, nanmean(tmp_shuff,2)]; % Mean across all events
        end
        mu = nanmean(mat_shuff(:)); sigma = nanmean(nanstd(mat_shuff,[],2)); % Mean, Std across number of shuffles
        sta_z{y} = (sta_raw{y} - nanmean(mu))./nanmean(sigma); % Z-score using mu, sigma from shuffled STA
        sta_shuff{y} = mat_shuff;
    end
    %% Load into output structure
    mat(x).rec = sub(x).rec; mat(x).n = sub(x).n; mat(x).FPnames{1} = beh(idx).FPnames{1};
    mat(x).sta_mvmt = sta_raw{1}; mat(x).sta_rest = sta_raw{2};
    mat(x).staZ_mvmt = sta_z{1}; mat(x).staZ_rest = sta_z{2};
    mat(x).shuff_mvmt = sta_shuff{1}; mat(x).shuff_rest = sta_shuff{2};
    
    waitbar(x/length(sub),h);
end
close(h); fprintf('Done: aligning photometry to CIN burst onsets \n');
if isempty(mat(10).n); mat(10) = []; end

%% (3) Extract from output structure
sta_mvmt = []; sta_rest = [];
staZ_mvmt = []; staZ_rest = []; shuff_mvmt = []; shuff_rest = [];
for x = 1:length(mat)
    sta_mvmt = [sta_mvmt, nanmean(mat(x).sta_mvmt,2)]; sta_rest = [sta_rest, nanmean(mat(x).sta_rest,2)];
    staZ_mvmt = [staZ_mvmt, nanmean(mat(x).staZ_mvmt,2)]; staZ_rest = [staZ_rest, nanmean(mat(x).staZ_rest,2)];
%     shuff_mvmt = [shuff_mvmt, (nanmean(mat(x).shuff_mvmt,2) - nanmean(mat(x).shuff_mvmt(:)))./nanstd(mat(x).shuff_mvmt(:))];
%     shuff_rest = [shuff_rest, (nanmean(mat(x).shuff_rest,2) - nanmean(mat(x).shuff_rest(:)))./nanstd(mat(x).shuff_rest(:))];
end

%% (4) Plot Average z-score
figure;
%shadederrbar(time, nanmean(shuff_mvmt,2), SEM(shuff_mvmt,2), 'k'); 
%shadederrbar(time, nanmean(shuff_rest,2), SEM(shuff_rest,2), 'k'); 
shadederrbar(time, nanmean(staZ_rest,2), SEM(staZ_rest,2), 'r'); hold on
shadederrbar(time, nanmean(staZ_mvmt,2), SEM(staZ_mvmt,2), 'g'); ylabel('ACh Fluorescence (z-score)');
% shadederrbar(time, nanmean(sta_rest,2), SEM(sta_rest,2), 'r'); hold on
% shadederrbar(time, nanmean(sta_mvmt,2), SEM(sta_mvmt,2), 'g'); ylabel('ACh Fluorescence (dF/F)');
xlabel('Latency to CIN burst (s)'); grid on; 
title(sprintf('Photometry aligned to CIN burst (n = %d units)',length(mat)));

%% (5) delta MOV vs REST change from t = 0
val_pkAdj = []; t_0 = find(time == 0);
for x = 1:size(staZ_mvmt,2)
    val_pkAdj(x,1) = max(staZ_mvmt(:,x)) - staZ_mvmt(t_0,x);
    val_pkAdj(x,2) = max(staZ_rest(:,x)) - staZ_rest(t_0,x);
end
val_diff = val_pkAdj(:,1)./val_pkAdj(:,2);
figure; hold on
plot([0.6:0.1:1.4],ones(9,1),'--k');
violinplot(val_diff(:)); grid on; ylabel('Fold Increase of Max STA - STA(t=0)'); xlim([0.6 1.4]);
title(sprintf('ACh aligned to CIN burst during MOV - REST\n mu = %1.3f, med = %1.3f (p = %1.3f, n = %d units)',mean(val_diff),median(val_diff),signrank(val_diff),length(val_diff)))

%% SUBPLOTS
fig = figure; plm = floor(sqrt(length(mat))); pln = ceil(length(mat)/plm);
for x = 1:length(mat)
    if isempty(mat(x).sta_mvmt); continue; end
    sp(x) = subplot(plm,pln,x); 
    shadederrbar(time, nanmean(mat(x).sta_mvmt,2), SEM(mat(x).sta_mvmt,2), 'g'); hold on
    shadederrbar(time, nanmean(mat(x).sta_rest,2), SEM(mat(x).sta_rest,2), 'r');
    title(sprintf('%s-%d',mat(x).rec, mat(x).n));
end

%%
mov = nan(length(time),length(align_m)); rest = mov;
for x = 5:length(align_m)
mov(:,x) = movmean(nanmean(align_m{x},2),1);
rest(:,x) = movmean(nanmean(align_r{x},2),1);
end

figure;
subplot(1,2,1);
shadederrbar(time, nanmean(rest,2), SEM(rest,2), 'r'); hold on
shadederrbar(time, nanmean(mov,2), SEM(mov,2), 'g');
xlabel('Latency to CIN burst (s)'); grid on;
ylabel('dF/F (z-score)');
title(sprintf('Photometry aligned to CIN st (n = %d)',size(mov,2)));

subplot(1,2,2);
m_max = max(mov) - mov(51,:);
r_max = max(rest) - rest(51,:);
violinplot([m_max; r_max]');
ylabel('mov - rest amplitude increase post burst');
title(sprintf('p = %1.4f',signrank([m_max - r_max])));
