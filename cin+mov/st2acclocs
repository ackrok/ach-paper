%% Input Variables
beh = behwt(10:12); sub = cinwt(20:29); %DMS recordings
beh = behwt([1:9,13:28,36:40]); sub = cinwt([1:19,30:98]); %DLS recordings

%% Extract Peaks
if all(logical(~rem(beh(1).on,1))); diffFs = 1; else; diffFs = 50; end
for x = 1:length(beh)
	locs = [];
	for y = 1:length(beh(x).on)
        vel_seg = beh(x).vel(diffFs*beh(x).on(y):diffFs*beh(x).off(y));
        acc_seg = [vel_seg(1); diff(movmean(vel_seg,10))];
        [~,locs_seg] = findpeaks(acc_seg,'MinPeakProminence',0.5,'MinPeakDistance',1);
        time_seg = beh(x).time(diffFs*beh(x).on(y):diffFs*beh(x).off(y));
        locs_seg = time_seg(locs_seg);
        locs = [locs; locs_seg(:)];
    end
	beh(x).acc_locs = locs; clc
end; clc

%% Align CIN spikes to acceleration peaks
mat = struct; % Initialize output structure
gen = struct; gen.bin = 0.01; gen.window = [-1 1]; nShuff = 10; %CHANGE: window for PETH 
h = waitbar(0, 'PETH: CIN spikes to acceleration peaks');
for x = 1:length(beh)
    idx = find(strcmp({sub.rec},beh(x).rec)); % Find matching units from beh recording
    if isempty(idx); continue; end
    st = {sub(idx).st}; % Unit spike times
    events = beh(x).acc_locs; % Acceleration peaks
    peth = getClusterPETH(st,events,gen); % PETH: CIN spikes to acc peaks
    tmp_peth = peth.fr;
    tmpZ = []; tmp50 = []; tmp95 = []; fr = [];
    for y = 1:length(st)
        % st_new = shuffleST(st{y},nShuff);
        % st_new = shiftST(st{y},nShuff,1000/nShuff); 
        st_new = poissonSpikeGen(sub(idx(y)).fr, beh(x).time(end), nShuff);
        peth = getClusterPETH(st_new,events,gen);
        mu = nanmean(peth.fr(:)); sigma = nanstd(peth.fr(:)); % mean, sigma of shuffled PETH
        tmpZ(:,y) = (tmp_peth(:,y) - mu)./sigma;
        tmpPrc = prctile(peth.fr,[5 50 95],2); %5th, 50th, 95th percentile of shuffled PETH
        tmp50(:,y) = tmpPrc(:,2); tmp95(:,y) = tmpPrc(:,3);
        
        isiSub = [];
        for z = 1:length(beh(x).on)
            logicalIndexes = st{y} < beh(x).off(z)/Fs & st{y} > beh(x).on(z)/Fs;
            isiSub = [isiSub; diff(st{y}(logicalIndexes))];
        end
        fr(y) = 1/mean(isiSub);
    end
    mat(x).rec = beh(x).rec;
    mat(x).n = [sub(idx).n];
    mat(x).fr = fr;
    mat(x).peth = tmp_peth;
    mat(x).pethZ = tmpZ;
    mat(x).prc50 = tmp50;
    mat(x).prc95 = tmp95;
    waitbar(x/length(beh),h);
end; fprintf('Done! \n'); close(h);
time = peth.time;

%% MEAN deltaFR
sm = 1;
pethMat = [mat.peth]; prc50 = [mat.prc50]; prc95 = [mat.prc95]; fr = [mat.fr];
deltaMat = []; delta50 = []; delta95 = []; pethZ = [mat.pethZ];
for x = 1:size(pethMat,2)
    deltaMat = [deltaMat,movmean([(pethMat(:,x) - fr(x))./fr(x)], sm)];
    delta50 = [delta50,movmean([(prc50(:,x) - fr(x))./fr(x)], sm)];
    delta95 = [delta95,movmean([(prc95(:,x) - fr(x))./fr(x)], sm)];
end

figure;
shadederrbar(time, nanmean(delta50,2), nanmean(delta95-delta50,2), 'k'); hold on
shadederrbar(time, nanmean(deltaMat,2), SEM(deltaMat,2), 'b'); 
xlabel('Latency (s)'); grid on; xlim([-1 1]);
ylabel('deltaFR');
title(sprintf('CIN st aligned to acc peak (n = %d units)',size(deltaMat,2)));

%% PLOT proportion > 95% CI
above95 = []; pethMat = [mat.peth]; prc95 = [mat.prc95];
for x = 1:size(pethMat,2)
    above95(:,x) = pethMat(:,x) > prc95(:,x); %binary vector where CCG passed 95% confidence interval
end

figure; hold on
bar(time, 100*sum(above95,2)/size(above95,2),'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag (s)'); ylabel('Prop > 95% CI (%)');
title(sprintf('Proportion > 95p CI (n = %d units)',size(above95,2)))

%% SUBPLOT FR
sm = 10;
figure; plm = floor(sqrt(size(pethMat,2))); pln = ceil(size(pethMat,2)/plm);
for x = 1:size(pethMat,2)
    sp(x) = subplot(plm,pln,x);
    shadederrbar(time, movmean(prc50(:,x),sm), movmean(prc95(:,x)-prc50(:,x),sm), 'k'); hold on
    plot(time, movmean(pethMat(:,x),sm), 'b'); ylabel('FR (Hz)');
    xlabel('Latency (s)'); grid on; xlim([-1 1]);
    title(sprintf('%s-#%d',sub(x).rec,sub(x).n));
end
