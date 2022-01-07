thresholds = [0:0.1:20]; % Thresholds to test: graded values from min_fr to max_fr of this unit
tpr_cin = []; fpr_cin = []; % Initialize output matrix for True-Positive and False-Positive Rates

h = waitbar(0,'ROC: cin fr to predict locomotion');
for x = 1:length(sub) % Iterate over all units
%%
% tic
st = sub(x).st; % Spike times, in seconds, of this unit
inst_fr_bin = 1; % Bin size, in secondss, for inst firing rate
inst = get_inst_fr(st, floor(st(end)), inst_fr_bin, 30000); % Instantaneous firing rate

b = find(strcmp({beh.rec},sub(x).rec)); % Identify matching behavioral recording
if isempty(b); continue; end
ev_on = beh(b).on; % Extract event onset times from data structure
ev_off = beh(b).off; % End of event
if isempty(ev_on); continue; end
bin = 1; % Set bin size
nBins = length(inst); % Total number of bins
sampBins = reshape([1:length(inst)], bin, nBins); % Reshape sample vector to form bins

binEv = zeros(nBins,1); % Initialize binarized vector for events
for z = 1:length(ev_on)
   binEv(round(ev_on(z)):round(ev_off(z))) = 1; 
end

for y = 1:length(thresholds)
    thres = thresholds(y);
    binPos = zeros(nBins,1);
    for z = 1:nBins
        if max(inst(sampBins(:,z))) > thres; binPos(z) = 1; end
    end
    crossYes = length(find(binPos)); % Total number of threshold crossings
    crossNo = nBins - crossYes;
    crossIdx = find(inst > thres); % Index of threshold crossings
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

    tpr_cin(y,x) = TP/(TP + FN); % True Positive Rate =  TP / (TP + FN)
    fpr_cin(y,x) = FP/(FP + TN); % False Positive Rate =  FP / (FP + TN)
end
waitbar(x/length(sub),h)
end; close(h);
%%
figure; 
% plot(fpr, tpr); hold on; plot([0:1],[0:1],':k'); xlim([0 1]); ylim([0 1]);
for x = 1:length(sub)
    sp(x) = subplot(10,10,x); hold on
    plot(fpr(:,x), tpr(:,x)); plot([0:1],[0:1],':k');
     xlabel('False Positive Rate'); ylabel('True Positive Rate');
    title(sprintf('%s-#%d',sub(x).rec,sub(x).n));
end

%% Plot all posBins
crossBins = find(binPos);
figure; hold on
for y = 1:length(crossBins)
    plot(beh(x).time(sampBins(:,crossBins(y))), beh(x).FP{1}(sampBins(:,crossBins(y))));
end
%%
figure;
shadederrbar(nanmean(fpr_cin,2), nanmean(tpr_cin,2), SEM(tpr_cin,2),'b'); hold on;
shadederrbar(nanmean(fpr_fp,2), nanmean(tpr_fp,2), SEM(tpr_fp,2),'g'); hold on;
plot([0:1],[0:1],':k');