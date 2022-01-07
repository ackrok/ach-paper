thresholds = [0:1:30]; % High Thresholds to test: +10 to +30dF/F
tpr_fp = []; fpr_fp = []; % Initialize output matrix for True-Positive and False-Positive Rates

for x = 1:length(beh) % Iterate over all units
    tpr_tmp = []; fpr_tmp = [];
    for x2 = 1:length(beh(x).on) % Iterate over sub-recordings
%%
% tic
fp = beh(x).fp_prc1(:,x2); % Extract photometry signal from data structure
ev_on = beh(x).on{x2}; % Extract event onset times from data structure
ev_off = beh(x).off{x2}; % End of event
bin = 1; % Set bin size
nBins = length(fp)/bin; % Total number of bins
sampBins = reshape([1:length(fp)], bin, nBins); % Reshape sample vector to form bins

binEv = zeros(nBins,1); % Initialize binarized vector for events
for z = 1:length(ev_on)
   binEv(ev_on(z):ev_off(z)) = 1; 
end

for y = 1:length(thresholds)
    thres = thresholds(y);
    binPos = zeros(nBins,1);
    for z = 1:nBins
        if max(fp(sampBins(:,z))) > thres; binPos(z) = 1; end
    end
    crossYes = length(find(binPos)); % Total number of threshold crossings
    crossNo = nBins - crossYes;
    crossIdx = find(fp > thres); % Index of threshold crossings
%     posBins = cross > sampBins(1,:) & cross < sampBins(bin,:); % Test positive when crossing occurs within bin
%     posBins = sum(posBins, 1); posBins( posBins > 1 ) = 1; % Sum over all crossings, convert to binary

    overlap_YesYes = binPos(:) + binEv(:); % Identify overlap between threshold crossings and event bins
    crossYes_evYes = length(find(overlap_YesYes == 2)); % Number of threshold crossings where event is also occuring
    crossNo_evNo = length(find((~binPos(:) + ~binEv(:)) == 2)); % True negative
    crossNo_evYes = length(find(binEv)) - crossYes_evYes;
    crossYes_evNo = crossYes - crossYes_evYes;
    
    TP = crossYes_evYes; % cross yes, event yes = # correct crossings
    TN = crossNo_evNo;   % cross no, event no
    FN = crossNo_evYes;  % cross no, event yes =  # true - # correct crossings
    FP = crossYes_evNo;  % cross yes, event no = # crossings - # correct crossings

    tpr_tmp(y,x2) = TP/(TP + FN); % True Positive Rate =  TP / (TP + FN)
    fpr_tmp(y,x2) = FP/(FP + TN); % False Positive Rate =  FP / (FP + TN)
end
    end
    tpr_fp(:,x) = nanmean(tpr_tmp,2); % Average over all sub-recordings
    fpr_fp(:,x) = nanmean(fpr_tmp,2); % Average over all sub-recordings
toc
end
%%
figure; 
% plot(fpr, tpr); hold on; plot([0:1],[0:1],':k'); xlim([0 1]); ylim([0 1]);
for x = 1:length(beh)
    sp(x) = subplot(4,6,x); hold on
    plot(fpr_fp(:,x), tpr_fp(:,x)); plot([0:1],[0:1],':k');
     xlabel('False Positive Rate'); ylabel('True Positive Rate');
    title(sprintf('%s',beh(x).rec));
end

%% Plot all posBins
crossBins = find(binPos);
figure; hold on
for y = 1:length(crossBins)
    plot(beh(x).time(sampBins(:,crossBins(y))), beh(x).FP{1}(sampBins(:,crossBins(y))));
end
