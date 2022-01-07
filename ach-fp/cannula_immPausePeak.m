%%
thres_pause = -3;
width_pause = 5;
thres_peak = 3; % CHANGE -- threshold for peaks/pauses to exceed, after baseline
width_peak = 5; % CHANGE -- require that peaks/pauses are at least X samples long (1 sample = 20ms)

h = waitbar(0, 'find immobility pauses/peaks');
for y = 1:length(s) % iterate over infusion
    %%
    beh = s(y).s;
    mag = cell(length(beh),2); 
    dur = cell(length(beh),2); freq = [];
    ach2achpause = cell(length(beh),2);

    for x = 1:length(beh) % iterate over animal
        fp = beh(x).FP{1}; Fs = 50; % extract full photometry signal
        idx_inf = [s(y).win(x,1).*(50*60) : s(y).win(x,2).*(50*60)]'; % infusion window
        rewWindow = 50; % how many samples after reward delivery is the reward window
        idx_inf_rew = extractEventST(idx_inf, floor(beh(x).reward), floor(beh(x).reward)+rewWindow, 1); % identify sample during reward
        idx_inf_rest = extractEventST(idx_inf, beh(x).onRest, beh(x).offRest, 1); % index: infusion + rest
        idx_inf_rest = idx_inf_rest(~ismember(idx_inf_rest, idx_inf_rew)); % exclude reward, include rest
        fp = fp - s(y).base_mu_rest(x); % subtract mean of baseline recording to center around 0%
        % fp = fp - nanmean(fp);
        
    %% pause identification
        [pause_maxIdx, pause_good, pause_valmag, ~, pause_crossstart] = findIdxPause(fp(idx_inf_rest),'pause',thres_pause,width_pause);
        [peak_maxIdx, peak_good, peak_valmag, ~, peak_crossstart] = findIdxPause(fp(idx_inf_rew),'peak',thres_peak,width_peak);
        
    %% ACh to ACh pause
        [ach2achpause{x,1},sta_time] = getSTA(fp(idx_inf_rest), pause_maxIdx/Fs, Fs, [-2 2]);
        ach2achpause{x,2} = getSTA(fp(idx_inf_rew), peak_maxIdx/Fs, Fs, [-2 2]);
        
    %%
        mag{x,1} = pause_valmag(pause_good); % Pause magnitude
        mag{x,2} = peak_valmag(peak_good); % Peak magnitude
        dur{x,1} = diff(pause_crossstart(pause_good,:),1,2); % Duration of pause that exceeds threshold
        dur{x,2} = diff(peak_crossstart(peak_good,:),1,2); % Duration of pause that exceeds threshold
        freq(x,1) = (1/mean(diff(pause_maxIdx)))*Fs; % Frequency of maximum deflections
        freq(x,2) = (1/mean(diff(peak_maxIdx)))*Fs; % Frequency of maximum deflections
        
    waitbar(x/length(beh),h);
%%
    end
    s(y).criteria = [thres_pause thres_peak; width_pause width_peak];
    s(y).mag = mag; s(y).dur = dur; s(y).freq = freq;
    s(y).ach2achpause = ach2achpause;
    
end
close(h);