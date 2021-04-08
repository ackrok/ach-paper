%% (1) Load data
load('allACh_beh.mat'); beh = wt_ach(1:22); % Extract ACh recordings

%% Proper Rest Threshold
velThres = 0.25;
Fs = 50;
minRestTime = 4*Fs; minRunTime = 1*Fs;
timeThres = 4*Fs; timeShift = 0.5*Fs;
for x = 1:length(beh)
    for y = 1:size(beh(x).vel,2)
        vel = beh(x).vel(:,y); vel = abs(vel);
        [onsetInd,offsetInd] = getOnsetOffset(-vel,-velThres,minRunTime,minRestTime,1);
        [onsetInd,offsetInd] = adjOnsetOffset(onsetInd,offsetInd,timeThres,vel);
        onsetInd = onsetInd+timeShift; offsetInd = offsetInd-timeShift;
        beh(x).onRest{y} = onsetInd; beh(x).offRest{y} = offsetInd;
    end
end; fprintf('Done %s: new rest threshold = 0.25.\n',beh(1).FPnames{1});

%% (2) STA
Fs = 50; time = [-1:1/Fs:1]; nShuff = 50; % CHANGE: window for STA
peakProm = 0.2; peakDist = 0.5;

mat = struct; % Initialize structure
h = waitbar(0, 'STA: acceleration to FP peak times');
for x = 1:length(beh)
    sta_raw = cell(3,size(beh(x).fp,2)); sta_z = sta_raw; sta_shuff = sta_raw;
    for z = 1:size(beh(x).fp,2)
        %% Compute fp peaks
        fp = beh(x).fp(:,z); fp = (fp - min(fp)) / (max(fp) - min(fp)); % Min Max Normalization
        [pks,locs] = findpeaks(fp,'MinPeakProminence',peakProm,'MinPeakDistance',peakDist); 
        locs = beh(x).time(locs); %location of peaks from time signal
        
        %% Extract event times during rest & movement
        ev_sub = cell(1,3); %ev_sub_forShuff = cell(1,2);
        ev_sub{1} = extractEventST(locs, beh(x).on{z}/Fs, beh(x).off{z}/Fs, 0); % Event times during movement
        ev_sub{2} = extractEventST(locs, beh(x).onRest{z}/Fs, beh(x).offRest{z}/Fs, 0); % Event times during rest
        ev_sub{3} = locs;
%         ev_sub_forShuff{1} = extractEventST(locs, beh(x).on{z}/Fs, beh(x).off{z}/Fs, 1); % Event times during movement
%         ev_sub_forShuff{2} = extractEventST(locs, beh(x).onRest{z}/Fs, beh(x).offRest{z}/Fs, 1); % Event times during rest

        %% Select signal to align photometry peaks to
        %sig = [beh(x).vel(1,z); diff(movmean(beh(x).vel(:,z),10))]; lbl = 'Acceleration'; % Aligning acceleration signal
        sig = beh(x).fp(:,z); lbl = 'ACh Fluorescence'; % Aliging photometry signal

        %% 
        for y = 1:length(ev_sub) % Iterate over rest and movement
            if isempty(ev_sub{y}); continue; end
%             if length(ev_sub_forShuff{1}) < 2; continue; end
            sta_raw{y,z} = getSTA(sig, ev_sub{y}, Fs, [time(1), time(end)]);
%             st_new = shuffleST(ev_sub_forShuff{1}, nShuff); % Shuffle event times n times
%             mat_shuff = [];
%             for z = 1:nShuff
%                 tmp_shuff = getSTA(sig, st_new{z}, Fs, [time(1), time(end)]); % Aligning signal to shuffled event times
%                 mat_shuff = [mat_shuff, nanmean(tmp_shuff,2)]; % Mean across all events
%             end
%             mu = nanmean(mat_shuff(:)); sigma = nanmean(nanstd(mat_shuff,[],2)); % Mean, Std across number of shuffles
%             sta_z{y,z} = (sta_raw{y} - nanmean(mu))./nanmean(sigma); % Z-score using mu, sigma from shuffled STA
%             sta_shuff{y,z} = mat_shuff;
        end
    end
    %% Load into output structure
    mat(x).rec = beh(x).rec; mat(x).FPnames{1} = beh(x).FPnames{1};
    mat(x).sta_mvmt = []; mat(x).sta_rest = []; mat(x).sta_full = [];
    %mat(x).staZ_mvmt = []; mat(x).staZ_rest = [];
    for z = 1:size(beh(x).fp,2)
        mat(x).sta_mvmt = [mat(x).sta_mvmt, sta_raw{1,z}];
       	mat(x).sta_rest = [mat(x).sta_rest, sta_raw{2,z}];
        mat(x).sta_full = [mat(x).sta_full, sta_raw{3,z}];
%         mat(x).staZ_mvmt = [mat(x).staZ_mvmt, sta_z{1,z}]; 
%         mat(x).staZ_rest = [mat(x).staZ_rest, sta_z{2,z}]; 
%         mat(x).shuff_mvmt = [mat(x).shuff_mvmt, sta_shuff{1,z}]; 
%         mat(x).shuff_rest = [mat(x).shuff_rest, sta_shuff{2,z}]; 
    end
    waitbar(x/length(beh),h);
end
close(h); fprintf('Done: aligning acceleration to ACh fp events \n');

%% (3) Extract from output structure
sta_mvmt = []; sta_rest = []; 
for x = 1:length(mat)
    sta_mvmt = [sta_mvmt, nanmean(mat(x).sta_mvmt,2)];
    sta_rest = [sta_rest, nanmean(mat(x).sta_rest,2)];
end

% (4) Plot Average z-score
figure;
shadederrbar(time, nanmean(sta_rest,2), SEM(sta_rest,2), 'r'); hold on
shadederrbar(time, nanmean(sta_mvmt,2), SEM(sta_mvmt,2), 'g'); 
xlabel('Latency to ACh Peak (s)'); ylabel(sprintf('%s',lbl)); grid on; 
title(sprintf('%s aligned to ACh Peaks (n = %d recs)',lbl,size(sta_rest,2)));
