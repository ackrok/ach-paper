y = 1; %aCSF
flag = 'peak';
thres = [];
width = 5;

for x = 1:4

fp = beh(x).FP{1}; Fs = 50; % extract full photometry signal
idx_inf = [s(y).win(x,1).*(50*60) : s(y).win(x,2).*(50*60)]'; % infusion window
rewWindow = 50; % how many samples after reward delivery is the reward window
idx_inf_rew = extractEventST(idx_inf, floor(beh(x).reward), floor(beh(x).reward)+rewWindow, 1); % identify sample during reward
idx_inf_rest = extractEventST(idx_inf, beh(x).onRest, beh(x).offRest, 1); % index: infusion + rest
idx_inf_rest = idx_inf_rest(~ismember(idx_inf_rest, idx_inf_rew)); % exclude reward, include rest
fp = fp - s(y).base_mu_rest(x); % subtract mean of baseline recording to center around 0%
% fp = fp - nanmean(fp);
fpInputVec = fp;

%%
thresholds = [2:8]; % thresholds: defined as deviation from new baseline
mu = []; sigma = []; freq = [];
for t = 1:length(thresholds) % Iterate over possible thresholds
thres = thresholds(t);

switch flag
case 'pause'; cross = find(fpInputVec < thres); % Find all samples below threshold
case 'peak'; cross = find(fpInputVec > thres); % Find all samples below threshold
end
[crossGroups, crossStart] = consecutive_vec2cell(cross); % Extract grouped crossings
valMag = nan(length(crossGroups),1); % Initialize vector for amplitudes
idxGood = []; idxMax = [];
for c = 1:length(crossGroups) % Iterate over discrete pause/peaks
fp_vec = [fpInputVec(crossGroups{c})]; % Vector for photometry values that cross during this pause/peak
if length(crossGroups{c}) > width % Only find amplitudes for peaks that are X samples long (1 sample = 20ms if Fs = 50Hz)
    switch flag
        case 'pause'; [valMag(c),ii] = min(fp_vec); % Pause magnitude
        case 'peak'; [valMag(c),ii] = max(fp_vec); % Peaks magnitude
    end
    idxGood = [idxGood; c]; % Index of crosses that satisfy parameters (width)
    idxMax = [idxMax; crossGroups{c}(ii)]; % Index of max point in cross
end
end

    mu(t) = nanmean(valMag(idxGood));
    sigma(t) = nanstd(valMag(idxGood));
    freq(t) = (1/mean(diff(idxMax)))*Fs;
end

%%
[align, time, ev] = plot_fp2event(beh, [-6 2], 0);
idx_inf = [s(y).win(x,:)*60]; % infusion window
ii = find(ev{x} > idx_inf(1) & ev{x} < idx_inf(end)); % index, wrt ev vector, where events are in infusion window
a_rew = align{x,1}(:,ii); % extract aligned traces that are in infusion window
a_rew = a_rew - s(y).base_mu_rest(x); % subtract mean of baseline recording to center around 0%
switch flag
    case 'pause'; m = min(a_rew([find(time == 0):find(time == 1)],:)); % Pause magnitude
    case 'peak'; m = max(a_rew([find(time == 0):find(time == 1)],:)); % Peaks magnitude
end
%%
figure; errorbar(mu, sigma); hold on; 
errorbar(0,nanmean(m),nanstd(m)); xlim([-1 10])
title(sprintf('%s',strtok(beh(x).rec,'_')));
end