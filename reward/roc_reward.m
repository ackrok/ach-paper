thresholds = [-5:0.1:2]; % Low Thresholds to test: -5 to +2 dF/F
% thresholds = [2:1:30]; % High Thresholds to test: +10 to +30dF/F

tpr = []; fpr = []; % Initialize output matrix for True-Positive and False-Positive Rates

for x = 1:length(beh)
%%
tic
fp = beh(x).FP{1}(1:round(beh(x).reward+beh(x).Fs,-1)); % Extract photometry signal from data structure
ev = beh(x).reward; % Extract reward delivery times from data structure
bin = beh(x).Fs/50; % Set bin size
nBins = length(fp)/bin; % Total number of bins
sampBins = reshape([1:length(fp)], bin, nBins); % Reshape sample vector to form bins

for y = 1:length(thresholds)
    thres = thresholds(y);
    posBins = [];
    for z = 1:nBins
        if min(fp(sampBins(:,z))) < thres; posBins(z) = 1; end
    end
     cross = find(fp < thres); % Find all samples below threshold
%     posBins = cross > sampBins(1,:) & cross < sampBins(bin,:); % Test positive when crossing occurs within bin
%     posBins = sum(posBins, 1); posBins( posBins > 1 ) = 1; % Sum over all crossings, convert to binary
    crossYes = length(find(posBins)); 
    posRew = cross > ev & cross < ev+bin; % True positive if crossing occurs within reward bin
    posRew = sum(posRew, 1); % Sum over all crossings
    
    cross = find(fp < thres); 
    
    TP = length(find(posRew));  % cross yes, reward yes = # correct crossings
    FN = length(ev) - TP;       % cross no, reward yes =  # true - # correct crossings
    FP = crossYes - TP;         % cross yes, reward no = # crossings - # correct crossings
    TN = nBins - (TP + FN + FP); % cross no, reward no
    tpr(y,x) = TP/(TP + FN); % True Positive Rate =  TP / (TP + FN)
    fpr(y,x) = FP/(FP + TN); % False Positive Rate =  FP / (FP + TN)
end
toc
end
%%
figure; 
% plot(fpr, tpr); hold on; plot([0:1],[0:1],':k'); xlim([0 1]); ylim([0 1]);
for x = 1:length(beh)
    sp(x) = subplot(3,5,x); hold on
    plot(fpr(:,x), tpr(:,x)); plot([0:1],[0:1],':k');
     xlabel('False Positive Rate'); ylabel('True Positive Rate');
    title(sprintf('%s-%s',beh(x).rec,beh(x).site));
end

%% Plot all posBins
crossBins = find(posBins);
figure; hold on
for y = 1:length(crossBins)
    plot(beh(x).time(sampBins(:,crossBins(y))), beh(x).FP{1}(sampBins(:,crossBins(y))));
end
