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
mat = struct; %Initialize structure to save CCG output data into
uni = unique({sub.rec}); %Find unique recording IDs across all units
Fs = 50; %Sampling frequency for behavioral data
diffFs = 50;

for x = 1
    fprintf('%s \n',uni{x});
    idx = find(strcmp({sub.rec},uni{x})); %Find all units that match this unique recording ID
    idx_b = find(strcmp({beh.rec},uni{x})); %Find behavior data that matches this unique recording ID
    if isempty(idx_b); continue; end
    if length(idx) < 2; continue; end %If there are <2 units (thus no possible unit pairs), continue to next unique recording ID
    ccgst = rununitccg({sub(idx).st},beh(idx_b).on/diffFs,beh(idx_b).off/diffFs,beh(idx_b).onRest/diffFs,beh(idx_b).offRest/diffFs); 

    mat(x).rec = uni{x};
    
    mat(x).n = [ccgst.pairs.n]; mat(x).m = [ccgst.pairs.m];
    mat(x).fr = [];
    mat(x).ccg = [ccgst.pairs.ccg];
    mat(x).ccgZ = [ccgst.pairs.ccgZ];
    mat(x).shuffPrc = cell(length(ccgst.pairs),1); %5th, 50th, 95th percentile of shuffled CCG's
    
    mat(x).fr_mvmt = [];
    mat(x).ccg_mvmt = [ccgst.pairs.ccg_mvmt]; 
    mat(x).ccgZ_mvmt = [ccgst.pairs.ccgZ_mvmt];
    mat(x).shuffPrc_mvmt = cell(length(ccgst.pairs),1); %5th, 50th, 95th percentile of shuffled CCG's
    
    mat(x).fr_rest = [];
    mat(x).ccg_rest = [ccgst.pairs.ccg_rest]; 
    mat(x).ccgZ_rest = [ccgst.pairs.ccgZ_rest];
    mat(x).shuffPrc_rest = cell(length(ccgst.pairs),1); %5th, 50th, 95th percentile of shuffled CCG's
    
    mat(x).dist = [];
    for y = 1:length(ccgst.pairs)
        mat(x).fr(y) = 1/mean(diff(ccgst.times(ccgst.pairs(y).m).full));
        mat(x).fr_mvmt(y) = 1/mean(diff(extractEventST(ccgst.times(ccgst.pairs(y).m).full,beh(idx_b).on/diffFs,beh(idx_b).off/diffFs,0)));
        mat(x).fr_rest(y) = 1/mean(diff(extractEventST(ccgst.times(ccgst.pairs(y).m).full,beh(idx_b).onRest/diffFs,beh(idx_b).offRest/diffFs,0)));
        mat(x).shuffPrc{y} = prctile(ccgst.pairs(y).shuff,[5 50 95],2); %5th, 50th, 95th percentile of shuffled CCG's
        mat(x).shuffPrc_mvmt{y} = prctile(ccgst.pairs(y).shuff_mvmt,[5 50 95],2); %5th, 50th, 95th percentile of shuffled CCG's
        mat(x).shuffPrc_rest{y} = prctile(ccgst.pairs(y).shuff_rest,[5 50 95],2); %5th, 50th, 95th percentile of shuffled CCG's
%         xc_diff = sub(idx(ccgst.pairs(y).n)).coor(1) - sub(idx(ccgst.pairs(y).m)).coor(1); %Distance in x- or y-dimension
%         zc_diff = sub(idx(ccgst.pairs(y).n)).coor(2) - sub(idx(ccgst.pairs(y).m)).coor(2); %Distnace in z-dimension
%         mat(x).dist(y) = hypot(xc_diff,zc_diff);
    end
end
time = ccgst.lag;

%% EXTRACT FROM MAT
time = [-2:0.01:2];
ccg_full = []; ccg_mvmt = []; ccg_rest = [];
ccgDelta = []; ccg95 = []; ccg50 = [];
ccgDelta_mvmt = []; ccgDelta_95mvmt = []; ccgDelta_50mvmt = [];
ccgDelta_rest = []; ccgDelta_95rest = []; ccgDelta_50rest = [];

for x = 1:length(mat)
    if isempty(mat(x).ccgZ_rest); continue; end
    ccg_full = [ccg_full, mat(x).ccg];
    ccg_mvmt = [ccg_mvmt, mat(x).ccg_mvmt]; 
    ccg_rest = [ccg_rest, mat(x).ccg_rest];
    for y = 1:length(mat(x).fr_rest)
        tmp = (mat(x).ccg(:,y) - mat(x).fr(y))./mat(x).fr(y); %deltaFR = (rate - mu(rate))/mu(rate)
        ccgDelta = [ccgDelta, tmp];
        tmp = (mat(x).ccg_rest(:,y) - mat(x).fr_rest(y))./mat(x).fr_rest(y); %deltaFR = (rate - mu(rate))/mu(rate)
        ccgDelta_rest = [ccgDelta_rest, tmp];
        tmp = (mat(x).ccg_mvmt(:,y) - mat(x).fr_mvmt(y))./mat(x).fr_mvmt(y); %deltaFR = (rate - mu(rate))/mu(rate)
        ccgDelta_mvmt = [ccgDelta_mvmt, tmp];
        
        tmp_95 = (mat(x).shuffPrc{y}(:,3) - mat(x).fr(y))./mat(x).fr(y); 
        tmp_50 = (mat(x).shuffPrc{y}(:,2) - mat(x).fr(y))./mat(x).fr(y);
        ccg95 = [ccg95, tmp_95];
        ccg50 = [ccg50, tmp_50];
        
        tmp_95 = (mat(x).shuffPrc_rest{y}(:,3) - mat(x).fr_rest(y))./mat(x).fr_rest(y); 
        tmp_50 = (mat(x).shuffPrc_rest{y}(:,2) - mat(x).fr_rest(y))./mat(x).fr_rest(y);
        ccgDelta_95rest = [ccgDelta_95rest, tmp_95];
        ccgDelta_50rest = [ccgDelta_50rest, tmp_50];
        
        tmp_95 = (mat(x).shuffPrc_mvmt{y}(:,3) - mat(x).fr_mvmt(y))./mat(x).fr_mvmt(y); 
        tmp_50 = (mat(x).shuffPrc_mvmt{y}(:,2) - mat(x).fr_mvmt(y))./mat(x).fr_mvmt(y);
        ccgDelta_95mvmt = [ccgDelta_95mvmt, tmp_95];
        ccgDelta_50mvmt = [ccgDelta_50mvmt, tmp_50];
    end
end

%% PLOT Z-SCORE
figure;
for x = 1:length(mat)
    if isempty(mat(x).ccgZ_rest); continue; end
    sp(x) = subplot(3,3,x); plot(time, [mat(x).ccgZ_rest], ':'); end
linkaxes(sp,'y');

figure; 
shadederrbar(time, nanmean(ccgZ,2), SEM(ccgZ,2), 'b'); 
xlabel('Lag (s)'); ylabel('CCG (z-score)'); xlim([-1 1]); grid on; 
title(sprintf('CCG rest spikes (n = %d pairs)',size(ccgZ,2)));

%% PLOT DELTA-FIRING RATE
figure; 
sm = 5;
shadederrbar(time, movmean(nanmean(ccgDelta_50rest,2),sm), movmean(nanmean(ccgDelta_95rest,2),sm), 'r'); hold on
shadederrbar(time, movmean(nanmean(ccgDelta_rest,2),sm), movmean(SEM(ccgDelta_rest,2),sm), 'b'); 
shadederrbar(time, movmean(nanmean(ccgDelta_50mvmt,2),sm), movmean(nanmean(ccgDelta_95mvmt,2),sm), 'g'); hold on
shadederrbar(time, movmean(nanmean(ccgDelta_mvmt,2),sm), movmean(SEM(ccgDelta_mvmt,2),sm), 'b'); 

%plot(time, movmean(ccgDelta_mvmt,sm,1), ':b');
xlabel('Lag (s)'); ylabel('CCG (deltaFR)'); xlim([-1 1]); grid on; 
title(sprintf('CCG spikes REST (n = %d pairs)',size(ccgDelta_rest,2)));

% figure; violinplot(max(ccgDelta,[],1)); xlim([0.6 1.4]); ylim([-0.1 1.5]);
% ylabel('Max deltaFR'); grid on
% title(sprintf('Max deltaFR (n = %d pairs)',length(max(ccgDelta,[],1))));

%% HEATMAP
figure;
a = ccgDelta_rest;
b = a(time == 0,:); % find delta @ lag = 0
[c, ii] = sort(b); % sort in ascending order
h = heatmap(a([find(time == -0.5):find(time == 0.5)],ii)');
h.Title = 'REST CIN pair CCG';
h.XLabel = 'Lag from reference spike (s)'; 
h.YLabel = 'Pair Number';
h.Colormap = jet; h.ColorLimits = [-0.4 1.4]; h.GridVisible = 'off';


%% PLOT versus DISTANCE
x_var = [mat.dist]'; %x_var = [mat.fr]';
max_ccg = max(ccgDelta,[],1)';
%%
tbl = table(x_var,max_ccg,'VariableNames',{'dist','maxCCG'}); %Create a table
mdl = fitlm(tbl,'maxCCG ~ dist'); %Fit a linear regression model with max(CCG) as response variable, unitDistance as predictor variable
ci = coefCI(mdl); %Find confidence intervals for the coefficients of the model
x_hat = [1:max(x_var)]; 
y_hat = mdl.Coefficients{1,1} + mdl.Coefficients{2,1}*x_hat; %Model fit y
yerrl = ci(1,1) + ci(2,1)*x_hat; yerru = ci(1,2) + ci(2,2)*x_hat; %Lower and upper confidence intervals
xpatch = [x_hat,fliplr(x_hat)]; ypatch = [yerru,fliplr(yerrl)]; %Generate patch for plotting shaded confidence interval

figure; hold on
plot(x_hat,nanmean(ccg95(:))*ones(length(x_hat),1),'--k'); %Dashed line for 95% confidence interval
plot(x_var,max_ccg,'.b','MarkerSize',10);
fill(xpatch,ypatch,[0 0 0],'FaceAlpha',0.25,'EdgeAlpha',0); %Plotting shaded confidence interval
plot(x_hat,y_hat,'-r'); %Plotting model fit line
xlabel('Unit Distance (um)'); ylabel('Max CCG rest spikes (deltaFR)');
title(sprintf('CCG peak ~ Unit Distance (n = %d pairs)',size(ccgDelta,2)));
title(sprintf('CCG peak ~ Unit Distance (n = %d pairs)\nR-squared = %1.3f',size(ccgDelta,2),mdl.Rsquared.Ordinary));
xlim([-50 1150]); 

%% PLOT proportion of pair CCG > 95% CI
above95 = []; below5 = [];
for x = 1:length(mat)
    for y = 1:size(mat(x).ccg_mvmt,2)
        a = mat(x).ccg_mvmt(:,y); b = mat(x).shuffPrc_mvmt{y};
        above95 = [above95, a > b(:,3)]; %binary vector where CCG passed 95% confidence interval
        below5 = [below5, a < b(:,1)]; %binary vector where CCG below 5% confidence interval
    end
end

figure; hold on
bar(time, 100*sum(above95,2)/size(above95,2),'FaceColor','b','FaceAlpha',0.5);
bar(time, -100*sum(below5,2)/size(below5,2),'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag (s)'); ylabel('Prop of Pair CCG > 95% CI (%)'); 
xlim([-1 1]); ylim([-50 100]);
title(sprintf('MVMT (n = %d pairs) %1.2f p > 95CI',size(above95,2),100*sum(above95(201,:),2)/size(above95,2)))

%% PLOT IV069_rec01 ONLY
x = 21; a = [12,23]; sm = 5;
fig = figure; fig.Position(3) = 1000;
for b = 1:2
    sp(b) = subplot(1,2,b); y = a(b);
    shadederrbar(ccgst.lag, movmean(mat(x).shuffPrc_rest{y}(:,2),sm), movmean(mat(x).shuffPrc_rest{y}(:,3)-mat(x).shuffPrc_rest{y}(:,2),sm), 'k'); hold on
    plot(ccgst.lag, movmean(mat(x).ccg_rest(:,y),5),'b');
    xlabel('Lag (s)'); ylabel('CCG (Firing Rate)'); xlim([-1 1]); grid on; ylim([3.5 7]);
    title(sprintf('IV069-rec01: run(#6)||ref(#%d)',b*2));
end

fig = figure; fig.Position(3) = 1000;
for b = 1:2
    sp(b) = subplot(1,2,b); y = a(b);
    tmp = (mat(x).ccg_rest(:,y) - mat(x).fr_rest(y))./mat(x).fr_rest(y);
    plot(ccgst.lag, movmean(tmp,5),'b');
    xlabel('Lag (s)'); ylabel('CCG (delta FiringRate)'); xlim([-1 1]); grid on; ylim([-0.25 0.5]);
    title(sprintf('IV069-rec01: run(#6)||ref(#%d)',b*2));
end

%% PLOT IV069_rec01 ONLY
x = 21; a = [106:141];

figure;
shadederrbar(time, nanmean(ccgDelta_50rest(:,a),2), nanmean(ccgDelta_95rest(:,a),2), 'k'); hold on
%shadederrbar(ccgst.lag, nanmean(ccgDelta_rest(:,a),2), SEM(ccgDelta_rest(:,a),2), 'b');
plot(ccgst.lag, ccgDelta(:,a), ':b');
xlabel('Lag (s)'); ylabel('CCG (deltaFR)'); xlim([-1 1]); grid on; 
title('IV069-rec01: CCG rest spikes (n = 36 pairs)');

x_var = [mat(x).dist]';
max_ccg = max(ccgDelta(:,a),[],1)';
figure; hold on
plot(x_hat,nanmean(nanmean(ccgDelta_95(:,a)))*ones(length(x_hat),1),'--k'); %Dashed line for 95% confidence interval
plot(x_var, max_ccg,'.b','MarkerSize',10);
xlabel('Unit Distance (um)'); ylabel('Max CCG rest spikes (z-score)'); 
xlim([-50 1150]); 
title('IV069-rec01: CCG peak ~ Unit Distance (n = 36 pairs)');

tbl = table(x_var,max_ccg,'VariableNames',{'dist','maxCCG'}); %Create a table
mdl = fitlm(tbl,'maxCCG ~ dist'); %Fit a linear regression model with max(CCG) as response variable, unitDistance as predictor variable
ci = coefCI(mdl); %Find confidence intervals for the coefficients of the model
x_hat = [1:1000]; 
y_hat = mdl.Coefficients{1,1} + mdl.Coefficients{2,1}*x_hat; %Model fit y
yerrl = ci(1,1) + ci(2,1)*x_hat; yerru = ci(1,2) + ci(2,2)*x_hat; %Lower and upper confidence intervals
xpatch = [x_hat,fliplr(x_hat)]; ypatch = [yerru,fliplr(yerrl)]; %Generate patch for plotting shaded confidence interval
fill(xpatch,ypatch,[0 0 0],'FaceAlpha',0.25,'EdgeAlpha',0); %Plotting shaded confidence interval
plot(x_hat,y_hat,'-r'); %Plotting model fit line
title(sprintf('IV069-rec01: CCG peak ~ Unit Distance (n = 36 pairs)\nR-squared = %1.3f',mdl.Rsquared.Ordinary));
