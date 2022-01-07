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

%% CCG rest and mvmt spikes
% mat = struct; %Initialize structure to save CCG output data into
% uni = unique({sub.rec}); %Find unique recording IDs across all units
% Fs = 50; %Sampling frequency for behavioral data
% 
% for x = 1:length(uni)
%     idx = find(strcmp({sub.rec},uni{x})); %Find all units that match this unique recording ID
%     idx_b = find(strcmp({beh.rec},uni{x})); %Find behavior data that matches this unique recording ID
%     if length(idx) < 2; continue; end %If there are <2 units (thus no possible unit pairs), continue to next unique recording ID
%     ccgst = rununitccg({sub(idx).st},beh(idx_b).on/50,beh(idx_b).off/50,beh(idx_b).onRest/50,beh(idx_b).offRest/50); 
% 
%     mat(x).rec = uni{x};
%     
%     mat(x).fr_mvmt = [];
%     mat(x).ccg_mvmt = [ccgst.pairs.ccg_mvmt]; 
%     mat(x).ccgZ_mvmt = [ccgst.pairs.ccgZ_mvmt];
%     mat(x).shuffPrc_mvmt = cell(length(ccgst.pairs),1); %5th, 50th, 95th percentile of shuffled CCG's
%     
%     mat(x).fr_rest = [];
%     mat(x).ccg_rest = [ccgst.pairs.ccg_rest]; 
%     mat(x).ccgZ_rest = [ccgst.pairs.ccgZ_rest];
%     mat(x).shuffPrc_rest = cell(length(ccgst.pairs),1); %5th, 50th, 95th percentile of shuffled CCG's
%     
%     mat(x).dist = [];
%     for y = 1:length(ccgst.pairs)
%         mat(x).fr_mvmt(y) = 1/mean(diff(extractEventST(ccgst.times(ccgst.pairs(y).m).full,beh(idx_b).on/50,beh(idx_b).off/50,1)));
%         mat(x).fr_rest(y) = 1/mean(diff(extractEventST(ccgst.times(ccgst.pairs(y).m).full,beh(idx_b).onRest/50,beh(idx_b).offRest/50,1)));
%         mat(x).shuffPrc_mvmt{y} = prctile(ccgst.pairs(y).shuff_mvmt,[5 50 95],2); %5th, 50th, 95th percentile of shuffled CCG's
%         mat(x).shuffPrc_rest{y} = prctile(ccgst.pairs(y).shuff_rest,[5 50 95],2); %5th, 50th, 95th percentile of shuffled CCG's
%         xc_diff = sub(idx(ccgst.pairs(y).n)).coor(1) - sub(idx(ccgst.pairs(y).m)).coor(1); %Distance in x- or y-dimension
%         zc_diff = sub(idx(ccgst.pairs(y).n)).coor(2) - sub(idx(ccgst.pairs(y).m)).coor(2); %Distnace in z-dimension
%         mat(x).dist(y) = hypot(xc_diff,zc_diff);
%     end
% end
% time = ccgst.lag;

%% EXTRACT FROM MAT
ccgZ_mov = []; 
ccgDelta_mov = []; ccgDelta_95mov = []; ccgDelta_50mov = [];
ccgDelta_rest = []; ccgDelta_rest95 = []; ccgDelta_rest50 = [];

for x = 1:length(mat)
    if isempty(mat(x).ccgZ_mvmt); continue; end
    ccgZ_mov = [ccgZ_mov, mat(x).ccgZ_mvmt]; 
    
    for y = 1:length(mat(x).fr_mvmt)
        tmp = (mat(x).ccg_mvmt(:,y) - mat(x).fr_mvmt(y))./mat(x).fr_mvmt(y); %deltaFR = (rate - mu(rate))/mu(rate)
        ccgDelta_mov = [ccgDelta_mov, tmp];
        tmp_95 = (mat(x).shuffPrc_mvmt{y}(:,3) - mat(x).fr_mvmt(y))./mat(x).fr_mvmt(y); 
        tmp_50 = (mat(x).shuffPrc_mvmt{y}(:,2) - mat(x).fr_mvmt(y))./mat(x).fr_mvmt(y);
        ccgDelta_95mov = [ccgDelta_95mov, tmp_95];
        ccgDelta_50mov = [ccgDelta_50mov, tmp_50];
        
        tmp = (mat(x).ccg_rest(:,y) - mat(x).fr_rest(y))./mat(x).fr_rest(y); %deltaFR = (rate - mu(rate))/mu(rate)
        ccgDelta_rest = [ccgDelta_rest, tmp];
        tmp_95 = (mat(x).shuffPrc_rest{y}(:,3) - mat(x).fr_rest(y))./mat(x).fr_rest(y); 
        tmp_50 = (mat(x).shuffPrc_rest{y}(:,2) - mat(x).fr_rest(y))./mat(x).fr_rest(y);
        ccgDelta_rest95 = [ccgDelta_rest95, tmp_95];
        ccgDelta_rest50 = [ccgDelta_rest50, tmp_50];
    end
end

%% PLOT Z-SCORE
figure;
for x = 1:length(mat)
    if isempty(mat(x).ccgZ_mvmt); continue; end
    sp(x) = subplot(5,6,x); plot(time, [mat(x).ccgZ_mvmt], ':'); end
linkaxes(sp,'y');

figure; 
shadederrbar(time, nanmean(ccgZ_mov,2), SEM(ccgZ_mov,2), 'b'); 
xlabel('Lag (s)'); ylabel('CCG (z-score)'); xlim([-1 1]); grid on; 
title(sprintf('CCG mvmt spikes (n = %d pairs)',size(ccgZ_mov,2)));

%% PLOT DELTA-FIRING RATE
figure; hold on
shadederrbar(time, nanmean(ccgDelta_50mov,2), nanmean(ccgDelta_95mov,2), 'g'); hold on
shadederrbar(time, nanmean(ccgDelta_rest50,2), nanmean(ccgDelta_rest95,2), 'r'); hold on
% plot(time, nanmean(ccgDelta_rest95,2), '--r'); plot(time, nanmean(ccgDelta_95mov,2), '--g'); plot(time, zeros(length(time),1), ':k');

shadederrbar(time, nanmean(ccgDelta_rest,2), SEM(ccgDelta_rest,2), 'r'); 
shadederrbar(time, nanmean(ccgDelta_mov,2), SEM(ccgDelta_mov,2), 'g'); 
xlabel('Lag (s)'); ylabel('CCG (deltaFR)'); xlim([-1 1]); grid on; 
title(sprintf('CCG mvmt vs rest spikes (n = %d pairs)',size(ccgDelta_mov,2)));

figure; violinplot(max(ccgDelta_mov,[],1)-max(ccgDelta_rest,[],1)); 
xlim([0.6 1.4]); ylim([-0.5 1]); ylabel('MOV - REST deltaFR'); grid on
title(sprintf('MOV - REST deltaFR (n = %d pairs)',length(max(ccgDelta_mov,[],1))));

%% PLOT proportion of pair CCG > 95% CI
above95_mov = []; above95_rest = [];
for x = 1:length(mat)
    for y = 1:length(mat(x).fr_mvmt)
        a = mat(x).ccg_mvmt(:,y); b = mat(x).shuffPrc_mvmt{y};
        above95_mov = [above95_mov, a > b(:,3)]; %binary vector where MOV CCG passed 95% confidence interval
        a = mat(x).ccg_rest(:,y); b = mat(x).shuffPrc_rest{y};
        above95_rest = [above95_rest, a > b(:,3)]; %binary vector where REST CCG passed 95% confidence interval
    end
end

figure; hold on
bar(time, 100*sum(above95_rest,2)/size(above95_rest,2),'FaceColor','r','FaceAlpha',0.5,'DisplayName','REST');
bar(time, 100*sum(above95_mov,2)/size(above95_mov,2),'FaceColor','g','FaceAlpha',0.5,'DisplayName','MOV');
xlabel('Lag (s)'); ylabel('Prop of Pair CCG > 95% CI (%)'); legend
title(sprintf('Proportion of pair CCGs > 95p CI (n = %d pairs)',size(above95_rest,2)))

%% PLOT versus DISTANCE
unitdist = [mat.dist]';
max_ccg = max(ccgDelta_mov,[],1)';

tbl = table(unitdist,max_ccg,'VariableNames',{'dist','maxCCG'}); %Create a table
mdl = fitlm(tbl,'maxCCG ~ dist'); %Fit a linear regression model with max(CCG) as response variable, unitDistance as predictor variable
ci = coefCI(mdl); %Find confidence intervals for the coefficients of the model
x_hat = [1:max(unitdist)]; 
y_hat = mdl.Coefficients{1,1} + mdl.Coefficients{2,1}*x_hat; %Model fit y
yerrl = ci(1,1) + ci(2,1)*x_hat; yerru = ci(1,2) + ci(2,2)*x_hat; %Lower and upper confidence intervals
xpatch = [x_hat,fliplr(x_hat)]; ypatch = [yerru,fliplr(yerrl)]; %Generate patch for plotting shaded confidence interval

figure; hold on
plot(x_hat,nanmean(ccgDelta_95mov(:))*ones(length(x_hat),1),'--k'); %Dashed line for 95% confidence interval
plot(unitdist,max_ccg,'.b','MarkerSize',10);
fill(xpatch,ypatch,[0 0 0],'FaceAlpha',0.25,'EdgeAlpha',0); %Plotting shaded confidence interval
plot(x_hat,y_hat,'-r'); %Plotting model fit line
xlabel('Unit Distance (um)'); ylabel('Max CCG mov spikes (deltaFR)');
title(sprintf('CCG peak ~ Unit Distance (n = %d pairs)',size(ccgZ_mov,2)));
title(sprintf('CCG peak ~ Unit Distance (n = %d pairs)\nR-squared = %1.3f',size(ccgZ_mov,2),mdl.Rsquared.Ordinary));
xlim([-50 1150]); 

%% PLOT IV069_rec01 ONLY
x = 21; a = [12,23]; sm = 5;
fig = figure; fig.Position(3) = 1000;
for b = 1:2
    sp(b) = subplot(1,2,b); y = a(b);
    shadederrbar(time, movmean(mat(x).shuffPrc_mvmt{y}(:,2),sm), movmean(mat(x).shuffPrc_mvmt{y}(:,3)-mat(x).shuffPrc_mvmt{y}(:,2),sm), 'k'); hold on
    plot(time, movmean(mat(x).ccg_mvmt(:,y),5),'b');
    xlabel('Lag (s)'); ylabel('CCG (Firing Rate)'); xlim([-1 1]); grid on; ylim([3.5 7]);
    title(sprintf('IV069-rec01: run(#6)||ref(#%d)',b*2));
end

fig = figure; fig.Position(3) = 1000;
for b = 1:2
    sp(b) = subplot(1,2,b); y = a(b);
    tmp = (mat(x).ccg_mvmt(:,y) - mat(x).fr_mvmt(y))./mat(x).fr_mvmt(y);
    plot(time, movmean(tmp,5),'b');
    xlabel('Lag (s)'); ylabel('CCG (delta FiringRate)'); xlim([-1 1]); grid on; ylim([-0.25 0.5]);
    title(sprintf('IV069-rec01: run(#6)||ref(#%d)',b*2));
end

%% PLOT IV069_rec01 ONLY
x = 21; a = [106:141];
figure;
shadederrbar(time, nanmean(ccgDelta_50mov(:,a),2), nanmean(ccgDelta_95mov(:,a),2), 'k'); hold on
%shadederrbar(time, nanmean(ccgDelta_rest(:,a),2), SEM(ccgDelta_rest(:,a),2), 'b');
plot(time, ccgDelta_mov(:,a), ':b');
xlabel('Lag (s)'); ylabel('CCG (deltaFR)'); xlim([-1 1]); grid on; 
title('IV069-rec01: CCG mvmt spikes (n = 36 pairs)');

unitdist = [mat(x).dist]';
max_ccg = max(ccgDelta_mov(:,a),[],1)';
figure; hold on
plot(x_hat,nanmean(nanmean(ccgDelta_95mov(:,a)))*ones(length(x_hat),1),'--k'); %Dashed line for 95% confidence interval
plot(unitdist, max_ccg,'.b','MarkerSize',10);
xlabel('Unit Distance (um)'); ylabel('Max CCG rest spikes (z-score)'); 
xlim([-50 1150]); 
title('IV069-rec01: CCG peak ~ Unit Distance (n = 36 pairs)');

tbl = table(unitdist,max_ccg,'VariableNames',{'dist','maxCCG'}); %Create a table
mdl = fitlm(tbl,'maxCCG ~ dist'); %Fit a linear regression model with max(CCG) as response variable, unitDistance as predictor variable
ci = coefCI(mdl); %Find confidence intervals for the coefficients of the model
x_hat = [1:1000]; 
y_hat = mdl.Coefficients{1,1} + mdl.Coefficients{2,1}*x_hat; %Model fit y
yerrl = ci(1,1) + ci(2,1)*x_hat; yerru = ci(1,2) + ci(2,2)*x_hat; %Lower and upper confidence intervals
xpatch = [x_hat,fliplr(x_hat)]; ypatch = [yerru,fliplr(yerrl)]; %Generate patch for plotting shaded confidence interval
fill(xpatch,ypatch,[0 0 0],'FaceAlpha',0.25,'EdgeAlpha',0); %Plotting shaded confidence interval
plot(x_hat,y_hat,'-r'); %Plotting model fit line
title(sprintf('IV069-rec01: CCG peak ~ Unit Distance (n = 36 pairs)\nR-squared = %1.3f',mdl.Rsquared.Ordinary));
