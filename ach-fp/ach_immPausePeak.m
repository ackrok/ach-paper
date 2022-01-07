s = struct;
y = 1;

thres_pause = -2;
width_pause = 7;
thres_peak = 2; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
width_peak = 7; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)

mag = cell(length(beh),2); 
amp = cell(length(beh),2); 
freq = [];
ach2ach = cell(length(beh),2);

%%
for x = 1:length(beh) % iterate over animal
    fp = beh(x).FP{1}; Fs = 50; % extract full photometry signal
    idx = [1:length(fp)]'; 
    if ~isempty(beh(x).reward)
        rewWindow = 50; % how many samples after reward delivery is the reward window
        idx_rew = extractEventST(idx, floor(beh(x).reward), floor(beh(x).reward)+rewWindow, 1); % identify sample during reward
    else; idx_rew = []; end
    idx_imm = extractEventST(idx, beh(x).onRest, beh(x).offRest, 1); % index: infusion + rest
    idx_imm = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    fp = fp - nanmean(fp); % subtract mean of recording to center around 0%

%% pause identification
    [pause_maxIdx, pause_good, pause_valmag, ~, pause_crossstart] = findIdxPause(fp(idx_imm),'pause',thres_pause,width_pause);
    [peak_maxIdx, peak_good, peak_valmag, ~, peak_crossstart] = findIdxPause(fp(idx_imm),'peak',thres_peak,width_peak);

%% ACh to ACh pause
    [ach2ach{x,1},sta_time] = getSTA(fp(idx_imm), pause_maxIdx/Fs, Fs, [-2 2]);
    ach2ach{x,2} = getSTA(fp(idx_imm), peak_maxIdx/Fs, Fs, [-2 2]);

%%
    mag{x,1} = pause_valmag(pause_good); % Pause magnitude
    mag{x,2} = peak_valmag(peak_good); % Peak magnitude
    amp{x,1} = diff(pause_crossstart(pause_good,:),1,2); % Duration of pause that exceeds threshold
    amp{x,2} = diff(peak_crossstart(peak_good,:),1,2); % Duration of pause that exceeds threshold
    freq(x,1) = (1/mean(diff(pause_maxIdx)))*Fs; % Frequency of maximum deflections
    freq(x,2) = (1/mean(diff(peak_maxIdx)))*Fs; % Frequency of maximum deflections

%%
end

%% by animal
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
amp_an = cell(nAn,2); dur_an = cell(nAn,2); freq_an = []; % Initialize cell arrays for by animal grouping
ach2ach_an = cell(2,1);

for x = 1:nAn
    for z = 1:2 % iterate over peaks, pauses
        ii = find(strcmp(tmp,uni{x}));
        ach2ach_an{z}(:,x) = nanmean([ach2ach{ii,z}],2);
        for y = 1:length(ii)
            amp_an{x,z} = [amp_an{x,z}; mag{ii(y),z}];
            dur_an{x,z} = [dur_an{x,z}; amp{ii(y),z}];
        end
        freq_an(x,z) = nanmean(freq(ii,z));
    end
end

mag_an_mu = [cellfun(@nanmean, amp_an)]; 
dur_an_mu = [cellfun(@nanmean, dur_an)]; 
dur_an_mu = dur_an_mu.*(1/Fs); % convert to seconds

%%
s(y).criteria = [thres_pause thres_peak; width_pause width_peak];
s(y).lbl = {'pause','peak'};
s(y).amp = mag_an_mu; 
s(y).dur = dur_an_mu; 
s(y).freq = freq_an;
s(y).ach2ach = ach2ach;
