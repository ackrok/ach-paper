%% (1) Load data
load('fullcin+beh_ACh')
sub = cinACh;  % Extract CINs from recordings with ACh
beh = behACh; % Extract ACh recordings

%% Proper Rest Threshold
velThres = 0.25;
Fs = 50;
minRestTime = 4*Fs; minRunTime = 1*Fs;
timeThres = 4*Fs; timeShift = 0.5*Fs;

for idx_b = 1:length(beh)
    vel = beh(idx_b).vel; vel = abs(vel);
    [onsetInd,offsetInd] = getOnsetOffset(-vel,-velThres,minRunTime,minRestTime,1);
    [onsetInd,offsetInd] = adjOnsetOffset(onsetInd,offsetInd,timeThres,vel);
    onsetInd = onsetInd+timeShift; offsetInd = offsetInd-timeShift;
    beh(idx_b).onRest = onsetInd; 
    beh(idx_b).offRest = offsetInd;
    
    if all(logical(~rem(beh(idx_b).on,1))); diffFs = 1; 
    else; diffFs = 50; 
        beh(idx_b).on = round(beh(idx_b).on.*diffFs); 
        beh(idx_b).off = round(beh(idx_b).off.*diffFs);
    end % Adjust for movement onset times being in samples or seconds
end

%% (2) PETH
Fs = 50; nShuff = 5;
bin = 0.02; window = [-1 1]; %CHANGE: window for PETH

mat = struct; % Initialize structure
h = waitbar(0, 'PETH: CIN spikes to FP peak times'); fprintf('Running ... ');
for x = 1:length(beh)
    fprintf('%s',beh(x).rec);
    % Extract spike times
    idx = find(strcmp({sub.rec},beh(x).rec));
    if isempty(idx); continue; end
    st = {sub(idx).st}; % Extract spike times of units from this recording
    fr = []; % Compute unit firing rate during movement and rest
    for y = 1:length(idx)
        fr(y,1) = 1/mean(diff(extractEventST(st{y},beh(x).on./Fs,beh(x).off./Fs,0))); % Event times during movement
        fr(y,2) = 1/mean(diff(extractEventST(st{y},beh(x).onRest./Fs,beh(x).offRest./Fs,0))); % Event times during rest
        fr(y,3) = 1/mean(diff(st{y}));
    end
    
    %% Compute fp peaks
    fp = beh(x).FP{1}; Fs = 50;
    fp = fp - nanmean(fp); % subtract mean of baseline recording to center around 0%
    thres_peak = 2; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
    width_peak = 4; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)
    locs = findIdxPause(fp,'peak',thres_peak,width_peak);
    locs = locs./50; % convert to seconds
    
    %% Extract event times during rest & movement
    ev_sub = cell(1,3); %ev_sub_forShuff = cell(1,2);
    ev_sub{1} = extractEventST(locs, beh(x).on./Fs, beh(x).off./Fs, 1); % Event times during movement
    ev_sub{2} = extractEventST(locs, beh(x).onRest./Fs, beh(x).offRest./Fs, 1); % Event times during rest
    ev_sub{3} = locs;
%     ev_sub_forShuff{1} = extractEventST(locs, beh(x).on/Fs, beh(x).off/Fs, 1); % Event times during movement
%     ev_sub_forShuff{2} = extractEventST(locs, beh(x).onRest/Fs, beh(x).offRest/Fs, 1); % Event times during rest

    %% PETH
    align = cell(1,length(ev_sub)); alignDelta = cell(1,length(ev_sub)); 
    alignZ = cell(1,length(ev_sub)); 
    shuff95 = alignZ; shuff50 = alignZ; shuff5 = shuff50;
    for z = 1:length(ev_sub)
        if isempty(ev_sub{z}); continue; end
        peth = getClusterPETH(st, ev_sub{z}, bin, window); % PETH: spike times aligned to fp peaks
        align{z} = peth.fr;
        for y = 1:length(idx)
            alignDelta{z} = [alignDelta{z}, (peth.fr(:,y)-fr(y,z))./fr(y,z)]; % Delta firing rate change
            stShuff = shuffleST(st{y}, nShuff);
            peth_shuff = getClusterPETH(stShuff, ev_sub{z}, bin, window); %PETH: shuffled spike times aligned to fp peaks
            mu = nanmean(nanmean(peth_shuff.fr,2)); sigma = nanmean(nanstd(peth_shuff.fr,[], 2)); 
            tmp_z = (peth.fr(:,y) - mu)./sigma; % z-score signal using shuffled mu, sigma
            alignZ{z} = [alignZ{z}, tmp_z]; % add to matrix
            prc = prctile(peth_shuff.fr,[5 50 95],2); %5th, 50th, 95th percentile of shuffled PETH
            shuff50{z} = [shuff50{z}, (prc(:,2)-fr(y,z))./fr(y,z)]; 
            shuff95{z} = [shuff95{z}, (prc(:,3)-fr(y,z))./fr(y,z)]; 
            shuff5{z} = [shuff5{z}, (prc(:,1)-fr(y,z))./fr(y,z)]; 
        end
    end
    
    %% STA
    sig = fp;
    sta_fp = cell(1,length(ev_sub));
    for z = 1:length(ev_sub)
        if isempty(ev_sub{z}); continue; end
        sta_fp{z} = getSTA(sig, ev_sub{z}, Fs, [window(1), window(2)]); % Align photometry to FP peaks
    end
    
    %% Load into output structure
    mat(x).rec = beh(x).rec; mat(x).FPnames{1} = beh(x).FPnames{1}; mat(x).n = [sub(idx).n]; mat(x).fr = fr;
    mat(x).align_mvmt = align{1}; mat(x).align_rest = align{2}; mat(x).align_full = align{3}; 
    mat(x).alignZ_mvmt = alignZ{1}; mat(x).alignZ_rest = alignZ{2}; mat(x).alignZ_full = alignZ{3};
    mat(x).alignDelta_mvmt = alignDelta{1}; mat(x).alignDelta_rest = alignDelta{2}; mat(x).alignDelta_full = alignDelta{3};
    mat(x).shuff5_mvmt = shuff5{1}; mat(x).shuff5_rest = shuff5{2}; mat(x).shuff5_full = shuff5{3};
    mat(x).shuff50_mvmt = shuff50{1}; mat(x).shuff50_rest = shuff50{2}; mat(x).shuff50_full = shuff50{3};
    mat(x).shuff95_mvmt = shuff95{1}; mat(x).shuff95_rest = shuff95{2}; mat(x).shuff95_full = shuff95{3};
    mat(x).sta_mvmt = sta_fp{1}; mat(x).sta_rest = sta_fp{2}; mat(x).sta_full = sta_fp{3};
    waitbar(x/length(beh),h); 
    fprintf(' ... ');
end
close(h); fprintf('\n Done: aligning spikes to ACh fp events \n');
time = peth.time; 

%% (3) Extract from output structure
align_mvmt = []; align_rest = [];
alignDelta_mvmt = []; alignDelta_rest = []; alignDelta_full = [];
%alignZ_mvmt = []; alignZ_rest = []; 
shuff50_rest = []; shuff95_rest = []; shuff5_rest = [];
for x = 1:length(mat)
    if isempty(mat(x).align_mvmt); continue; end
    align_mvmt = [align_mvmt, mat(x).align_mvmt]; % Concatenate all units
    align_rest = [align_rest, mat(x).align_rest];
    alignDelta_mvmt = [alignDelta_mvmt, mat(x).alignDelta_mvmt]; 
    alignDelta_rest = [alignDelta_rest, mat(x).alignDelta_rest];
    alignDelta_full = [alignDelta_full, mat(x).alignDelta_full];
%     alignZ_mvmt = [alignZ_mvmt, mat(x).alignZ_mvmt];
%     alignZ_rest = [alignZ_rest, mat(x).alignZ_rest];
    shuff95_rest = [shuff95_rest, mat(x).shuff95_rest];
    shuff50_rest = [shuff50_rest, mat(x).shuff50_rest];
    shuff5_rest = [shuff5_rest, mat(x).shuff5_rest];
end

%% (4) Plot Average z-score
figure;
% shadederrbar(time, nanmean(align_rest,2), SEM(align_rest,2), 'r'); hold on
% shadederrbar(time, nanmean(align_mvmt,2), SEM(align_mvmt,2), 'g'); ylabel('CIN Firing Rate (Hz)');
% shadederrbar(time, nanmean(alignZ_rest,2), SEM(alignZ_rest,2), 'r'); hold on
% shadederrbar(time, nanmean(alignZ_mvmt,2), SEM(alignZ_mvmt,2), 'g'); ylabel('CIN Firing Rate (z-score)');
shadederrbar(time, nanmean(shuff50_rest,2), nanmean(shuff5_rest,2), 'k'); hold on
shadederrbar(time, nanmean(alignDelta_rest,2), SEM(alignDelta_rest,2), 'r'); ylabel('CIN Firing Rate (deltaFR)'); hold on
% shadederrbar(time, nanmean(alignDelta_mvmt,2), SEM(alignDelta_mvmt,2), 'g'); ylabel('CIN Firing Rate (deltaFR)');
%shadederrbar(time, movmean(nanmean(alignDelta_rest(:,[5:end]),2),5), movmean(SEM(alignDelta_rest(:,[5:end]),2),5), 'b'); ylabel('CIN Firing Rate (deltaFR)');

% shadederrbar(time, nanmean(alignDelta_full,2), SEM(alignDelta_full,2), 'k'); ylabel('pMSN Firing Rate (deltaFR)');

xlabel('Latency to ACh Peak (s)'); grid on; 
title(sprintf('CIN spikes aligned to rest ACh Peaks (n = %d units)',size(align_rest,2)));

%% (4B) PLOT proportion of pair CCG > 95% CI
above95 = []; below5 = [];
for x = 1:length(mat)
    for y = 1:size(mat(x).alignDelta_rest,2)
        a = mat(x).alignDelta_rest(:,y); 
        b = mat(x).shuff95_rest(:,y);
        c = mat(x).shuff50_rest(:,y) - (mat(x).shuff95_rest(:,y) - mat(x).shuff50_rest(:,y));  % mat(x).shuff5_rest(:,y);
        above95 = [above95, a > b]; %binary vector where CCG passed 95% confidence interval
        below5 = [below5, a < c]; %binary vector where CCG below 5% confidence interval
    end
end

figure; hold on
bar(time, 100*sum(above95,2)/size(above95,2),'FaceColor','b','FaceAlpha',0.5);
bar(time, -100*sum(below5,2)/size(below5,2),'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag (s)'); ylabel('Prop of Pair CCG > 95% CI (%)'); 
ylim([-100 100]); xlim([-1 1]);
title(sprintf('Proportion of pair CCGs > 95p CI (n = %d pairs)',size(above95,2)))

%% (5) Representative ALIGN
figure; sm = 5;
x = find(strcmp({mat.rec},'IV043_rec02'));
yyaxis right; shadederrbar([-1:1/50:1], nanmean(mat(x).sta_rest,2), SEM(mat(x).sta_rest,2),'g'); ylabel('ACh (dF/F %)');
yyaxis left; 
%shadederrbar(time, movmean(mat(x).shuff50_rest(:,1),sm), movmean(mat(x).shuff95_rest(:,1)-mat(x).shuff50_rest(:,1),sm), 'k'); hold on
shadederrbar(time, movmean(mat(x).shuff50_rest(:,2),sm), movmean(mat(x).shuff95_rest(:,2)-mat(x).shuff50_rest(:,2),sm), 'k'); hold on
plot(time, movmean(mat(x).align_rest(:,2),sm), 'b','LineWidth',2); ylabel('CIN Firing Rate (Hz)');
xlabel('Latency to ACh Peak (s)'); 
title(sprintf('%s: CIN #%d - spikes aligned to rest ACh peaks',mat(x).rec,mat(x).n(2)));

%% (5b) Representative TRACE
% load('beh_wt_ACh+DA+PF.mat')
beh = behwt([1:19]); % Extract ACh recordings
idx_b = find(strcmp({beh.rec},'IV043_rec02'));
fp = beh(idx_b).FP{1}; fp_norm = (fp - min(fp)) / (max(fp) - min(fp)); % Min Max Normalization
[pks,locs] = findpeaks(fp_norm,'MinPeakProminence',0.05,'MinPeakDistance',0.5);
locs = beh(idx_b).time(locs);

figure; hold on; 
plot(beh(idx_b).time,movmean(fp,1),'g'); plot(locs,fp(round(locs*Fs)),'*k'); 
st = {sub(find(strcmp({sub.rec},'IV043_rec02'))).st};
for y = 1:length(st)
    plot(st{y}, ((y-1)/2)*ones(length(st{y}),1), '.b','MarkerSize',10);
end
xlabel('Time (s)'); ylabel('ACh Fluorescence (dF/F %)');
title(sprintf('%s - REST',beh(idx_b).rec),'Interpreter','none');
xlim([4715 4728]); ylim([-1 6.5]); yticks([0:2:6])

%% (6) DMS vs DLS
%% (6a) Extract 95% CI from mat
shuff50_DMS = []; shuff95_DMS = [];
for x = 1:length(matDMS)
shuff95_DMS = [shuff95_DMS, matDMS(x).shuff95_full];
shuff50_DMS = [shuff50_DMS, matDMS(x).shuff50_full];
end
shuff50_DLS = []; shuff95_DLS = [];
for x = 1:length(matDLS)
shuff95_DLS = [shuff95_DLS, matDLS(x).shuff95_full];
shuff50_DLS = [shuff50_DLS, matDLS(x).shuff50_full];    
end

%% (6b) PETH: CIN activity to ACh peaks -- DMS vs DLS
figure; hold on

% plot(time, movmean(nanmean(shuff95_DLS-shuff50_DLS,2),5), 'b');
% plot(time, movmean(nanmean(shuff95_DMS-shuff50_DMS,2),5), 'k');
% plot(time, -1*movmean(nanmean(shuff95_DLS-shuff50_DLS,2),5), 'b');
% plot(time, -1*movmean(nanmean(shuff95_DMS-shuff50_DMS,2),5), 'k');
xpatch = [time',fliplr(time')];
ypatch_dls = [movmean(nanmean(shuff95_DLS-shuff50_DLS,2),5)', fliplr(-1*movmean(nanmean(shuff95_DLS-shuff50_DLS,2),5))'];
ypatch_dms = [movmean(nanmean(shuff95_DMS-shuff50_DMS,2),5)', fliplr(-1*movmean(nanmean(shuff95_DMS-shuff50_DMS,2),5))'];
color_dls = char2rgb('b'); patchcolor_dls = color_dls+(1-color_dls)*.8;
color_dms = char2rgb('k'); patchcolor_dms = color_dms+(1-color_dms)*.8;
fill(xpatch,ypatch_dls,patchcolor_dls,'FaceAlpha',0.5,'EdgeColor',patchcolor_dls,'LineStyle','-');
fill(xpatch,ypatch_dms,patchcolor_dms,'FaceAlpha',0.5,'EdgeColor',patchcolor_dms,'LineStyle','-');

shadederrbar(time, movmean(nanmean([matDLS.alignDelta_full],2),5), movmean(SEM([matDLS.alignDelta_full],2),5), 'b'); hold on
shadederrbar(time, movmean(nanmean([matDMS.alignDelta_full],2),5), movmean(SEM([matDMS.alignDelta_full],2),5), 'k'); 
title('IV043: DLS (blue/n = 15 units) vs DMS (black/n = 10 units)');

% shadederrbar(time, movmean(nanmean([matDMS.alignDelta_rest],2),5), movmean(SEM([matDMS.alignDelta_rest],2),5), 'r'); hold on
% shadederrbar(time, movmean(nanmean([matDMS.alignDelta_mvmt],2),5), movmean(SEM([matDMS.alignDelta_mvmt],2),5), 'g'); 
%title(sprintf('IV043: DMS (n = 10 units) REST (red) vs MOV (green)'));

% shadederrbar(time, movmean(nanmean([matDLS.alignDelta_rest],2),5), movmean(SEM([matDLS.alignDelta_rest],2),5), 'r'); hold on
% shadederrbar(time, movmean(nanmean([matDLS.alignDelta_mvmt],2),5), movmean(SEM([matDLS.alignDelta_mvmt],2),5), 'g'); 
%title(sprintf('IV043: DLS (n = 15 units) REST (red) vs MOV (green)'));

ylabel('Firing Rate (deltaFR)'); xlabel('Latency to ACh Peak (s)'); grid on

%% (6c) STA: ACh photometry to ACh peaks
figure;
shadederrbar([-1:1/50:1], nanmean([matDMS.sta_full],2), SEM([matDMS.sta_full],2),'k'); hold on
shadederrbar([-1:1/50:1], nanmean([matDLS.sta_full],2), SEM([matDLS.sta_full],2),'b'); 
ylabel('ACh (dF/F %)'); xlabel('Latency to ACh Peak (s)'); grid on
title('IV043: DLS ACh photometry in DLS (blue) vs DMS (black) recordings');
