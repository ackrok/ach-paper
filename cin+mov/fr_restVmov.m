load('fullWT_2-10-21.mat')
sub = WT(find([WT.label] == 2)); %n = 125 CINs
beh = behwt;

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
    beh(idx_b).onRest = onsetInd; beh(idx_b).offRest = offsetInd;
end
%% Extract Firing Rates for ALL CINs
fr = []; % fr_mvmt = fr(:,1); fr_rest = fr(:,2); 
for x = 1:length(beh)
    idx = find(strcmp({sub.rec},beh(x).rec));
    if isempty(idx); continue; end
    times = cell(2,2); times{1,1} = beh(x).on/50; times{2,1} = beh(x).off/50; 
                times{1,2} = beh(x).onRest/50; times{2,2} = beh(x).offRest/50; 
    for y = 1:length(idx) % Repeat over all CINs in this recording
        st = sub(idx(y)).st; % Extract spike times for this CIN
        tmp_fr = [];
        for a = 1:2 % Repeat over movement, rest
            tmp_isi = []; % Clear/initialize ISI vector
            for z = 1:length(times{1,a}) % Iterate over all mov/rest periods
                range = [times{1,a}(z), times{2,a}(z)];  % Sample range
                check = st >= range(1) & st <= range(2); % Check if values within range
                stTmp = st;
                stTmp(~check) = nan; % If no values within range, return NaN
                stTmp = stTmp(~isnan(stTmp)); % Remove NaNs from spike time vector
                tmp_isi = [tmp_isi; diff(stTmp)]; % concatenate ISI between spike times within event range
            end
            tmp_fr(a) = 1/nanmean(tmp_isi); % Firing rate = 1/mean(diff(st))
        end
        fr = [fr; tmp_fr]; % Concatenate to output matrix
    end
end
fprintf('Done Extracting REST + MOV Firing Rates.\n');

%% PLOT firing rate scatter/bar plot
figure; hold on
plot([1*ones(1,length(fr));2*ones(1,length(fr))],[fr(:,1)';fr(:,2)'],':k');
for x = 1:2; clr = {'g','r'};
    scatter(x*ones(length(fr(:,x)),1),fr(:,x),clr{x});
    bar(x,nanmean(fr(:,x)),clr{x},'FaceAlpha',0.1); 
    errorbar(x,nanmean(fr(:,x)),SEM(fr(:,x),1),'k');
end
xticks([1:2]); xticklabels({'MOV','REST'}); xlim([0.25 2.75]); ylabel('Firing Rate (Hz)');
title(sprintf('Firing Rate MOV - REST (n = %d units)',length(fr)));

%% PLOT difference violin plot
figure; violinplot(fr(:,1)-fr(:,2));
xticks([1]); xticklabels({'MOV - REST'}); xlim([0.6 1.4]); grid on
ylabel('Firing Rate Difference (Hz)')
title(sprintf('Firing Rate MOV - REST: mu = %1.3f (p = %1.3f)',nanmean(fr(:,1)-fr(:,2)),signrank(fr(:,1)-fr(:,2))));

%%
diff_pos = length(find(fr(:,1)-fr(:,2) > 0)); 
diff_neg = length(find(fr(:,1)-fr(:,2) < 0)); 
