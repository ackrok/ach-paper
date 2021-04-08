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

%% CCG rest spikes - ATTEMPT 1
mat = struct; %Initialize structure to save CCG output data into
uni = unique({sub.rec}); %Find unique recording IDs across all units
Fs = 50; %Sampling frequency for behavioral data
for x = 1:length(uni)
    idx = find(strcmp({sub.rec},uni{x})); %Find all units that match this unique recording ID
    idx_b = find(strcmp({beh.rec},uni{x})); %Find behavior data that matches this unique recording ID
    if length(idx) < 2; continue; end %If there are <2 units (thus no possible unit pairs), continue to next unique recording ID
    st_rest = cell(length(idx),1);
    for y = 1:length(idx)
        st = sub(idx(y)).st; %Extract spike times for unit
        %if isempty(sub(idx(y)).burst); continue; end
        %st = sub(idx(y)).burst.Windows(:,1); %Extract burst onset times for unit
        for z = 1:length(beh(idx_b).onRest)
            logicalIndexes = st <= beh(idx_b).offRest(z)/Fs & st >= beh(idx_b).onRest(z)/Fs;
            st_rest{y} = [st_rest{y}; st(logicalIndexes)];
        end
    end
    ccgst = rununitccg(st_rest);

    mat(x).rec = uni{x};
    mat(x).fr = [sub(idx).fr];
    mat(x).ccg_rest = [ccgst.pairs.ccg]; 
    mat(x).ccgZ_rest = [ccgst.pairs.ccgZ];
    mat(x).dist = [];
    mat(x).mu = []; mat(x).sigma = []; 
    for y = 1:length(ccgst.pairs)
        xc_diff = sub(idx(ccgst.pairs(y).n)).coor(1) - sub(idx(ccgst.pairs(y).m)).coor(1); %Distance in x- or y-dimension
        zc_diff = sub(idx(ccgst.pairs(y).n)).coor(2) - sub(idx(ccgst.pairs(y).m)).coor(2); %Distnace in z-dimension
        mat(x).dist(y) = hypot(xc_diff,zc_diff);
        mat(x).mu(y) = nanmean(nanmean([ccgst.pairs(y).shuff]));
        mat(x).sigma(y) = nanmean(nanstd([ccgst.pairs(y).shuff],[],2));
    end
    mat(x).std = nanstd([ccgst.pairs.shuff],[],2);
end
time = ccgst.lag;

%% CCG rest spikes - ATTEMPT 2
mat = struct; %Initialize structure to save CCG output data into
uni = unique({sub.rec}); %Find unique recording IDs across all units
Fs = 50; %Sampling frequency for behavioral data
for x = 1:length(uni)
    idx = find(strcmp({sub.rec},uni{x})); %Find all units that match this unique recording ID
    idx_b = find(strcmp({beh.rec},uni{x})); %Find behavior data that matches this unique recording ID
    if length(idx) < 2; continue; end %If there are <2 units (thus no possible unit pairs), continue to next unique recording ID
    st = {sub(idx).st}; %Extract spike times for units in this recording into a cell array
    ccgst = rununitccg(st,beh(idx_b).on/Fs,beh(idx_b).off/Fs,beh(idx_b).onRest/Fs,beh(idx_b).offRest/Fs);
    mat(x).rec = uni{x};
    mat(x).fr = [sub(idx).fr];
    mat(x).st_rest = {ccgst.times.rest};
    mat(x).st_mov = {ccgst.times.mvmt};
    mat(x).ccg_rest = [ccgst.pairs.ccg_rest]; 
    mat(x).ccgZ_rest = [ccgst.pairs.ccgZ_rest];
    mat(x).ccg_mov = [ccgst.pairs.ccg_mvmt]; 
    mat(x).ccgZ_mov = [ccgst.pairs.ccgZ_mvmt];
    mat(x).dist = [];
    mat(x).mu = []; mat(x).sigma = []; 
    for y = 1:length(ccgst.pairs)
        xc_diff = sub(idx(ccgst.pairs(y).n)).coor(1) - sub(idx(ccgst.pairs(y).m)).coor(1); %Distance in x- or y-dimension
        zc_diff = sub(idx(ccgst.pairs(y).n)).coor(2) - sub(idx(ccgst.pairs(y).m)).coor(2); %Distnace in z-dimension
        mat(x).dist(y) = hypot(xc_diff,zc_diff);
        mat(x).mu(y) = nanmean(nanmean([ccgst.pairs(y).shuff]));
        mat(x).sigma(y) = nanmean(nanstd([ccgst.pairs(y).shuff],[],2));
    end
    mat(x).std = nanstd([ccgst.pairs.shuff],[],2);
end
time = ccgst.lag;

%%
ccgZ_rest = []; 
for x = 1:length(mat)
    if isempty(mat(x).ccgZ_rest); continue; end
    ccgZ_rest = [ccgZ_rest, mat(x).ccgZ_rest]; end
shuff = [mat.std]; 

%% PLOT ALL
figure;
for x = 1:length(mat)
    if isempty(mat(x).ccgZ_rest); continue; end
    sp(x) = subplot(5,6,x); plot(time, [mat(x).ccgZ_rest], ':'); end
linkaxes(sp,'y');

%% PLOT AVERAGE
figure; 
%shadederrbar(ccgst.lag, nanmean(shuff,2), nanstd(shuff,[],2), 'k'); 
shadederrbar(ccgst.lag, nanmean(ccgZ_rest,2), SEM(ccgZ_rest,2), 'b'); 
xlabel('Lag (s)'); ylabel('CCG (z-score)'); xlim([-1 1]); grid on; %ylim([-8 22]);
title(sprintf('CCG rest spikes (n = %d pairs)',size(ccgZ_rest,2)));

%% PLOT versus DISTANCE
unitdist = [mat.dist]';
max_ccg = max(ccgZ_rest,[],1)' - nanmean(ccgZ_rest([1:151],:),1)';
tbl = table(unitdist,max_ccg,'VariableNames',{'dist','maxCCG'}); %Create a table
mdl = fitlm(tbl,'maxCCG ~ dist'); %Fit a linear regression model with max(CCG) as response variable, unitDistance as predictor variable
ci = coefCI(mdl); %Find confidence intervals for the coefficients of the model
x = [1:max(unitdist)]; 
y_hat = mdl.Coefficients{1,1} + mdl.Coefficients{2,1}*x; %Model fit y
yerrl = ci(1,1) + ci(2,1)*x; yerru = ci(1,2) + ci(2,2)*x; %Lower and upper confidence intervals
xpatch = [x,fliplr(x)]; ypatch = [yerru,fliplr(yerrl)]; %Generate patch for plotting shaded confidence interval

figure; hold on
plot(x,std(base)*ones(length(x),1),'--k'); %Dashed line for noise
plot(unitdist,max_ccg,'.b','MarkerSize',10);
fill(xpatch,ypatch,[0 0 0],'FaceAlpha',0.25,'EdgeAlpha',0); %Plotting shaded confidence interval
plot(x,y_hat,'-r'); %Plotting model fit line
xlabel('Unit Distance (um)'); ylabel('Max CCG rest spikes (z-score)');
title(sprintf('CCG peak ~ Unit Distance (n = %d pairs)',size(ccgZ_rest,2)));
title(sprintf('CCG peak ~ Unit Distance (n = %d pairs)\nR-squared = %1.3f',size(ccgZ_rest,2),mdl.Rsquared.Ordinary));
xlim([-50 1150]); 

%% PLOT IV069_rec01 ONLY
figure;
%plot(ccgst.lag, mat(29).ccgZ, '--b');
shadederrbar(ccgst.lag, nanmean(nanmean(ccgZ_rest([1:100],[124:159]),2))*ones(1,401), std(nanmean(ccgZ_rest([1:100],[124:159])))*ones(1,401), 'k'); hold on
plot(ccgst.lag, ccgZ_rest(:,[124:159]), ':b');
xlabel('Lag (s)'); ylabel('CCG (z-score)'); xlim([-1 1]); grid on
title('IV069-rec01: CCG rest spikes (n = 36 pairs)');

unitdist = [mat(29).dist]';
max_ccg = max(ccgZ_rest(:,[124:159]),[],1)' - nanmean(ccgZ_rest([1:151],[124:159]),1)';
figure; hold on
plot(x,std(base)*ones(length(x),1),'--k'); %Dashed line for noise
plot(unitdist, max_ccg,'.b','MarkerSize',10);
xlabel('Unit Distance (um)'); ylabel('Max CCG rest spikes (z-score)'); 
xlim([-50 1150]); 
title('IV069-rec01: CCG peak ~ Unit Distance (n = 36 pairs)');
