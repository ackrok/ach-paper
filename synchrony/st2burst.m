%% (1)
lbl = 'burst'; %CHANGE: 'burst','pause',spike'

mat = struct; %Initialize structure to save output data into
uni = unique({sub.rec}); %Find unique recording IDs across all units
gen = struct; gen.bin = 0.01; gen.window = [-0.5 1];
Fs = 50; %Sampling frequency for behavioral data

for x = 1:length(uni)
    idx = find(strcmp({sub.rec},uni{x})); %Find all units that match this unique recording ID
    b = find(strcmp({beh.rec},uni{x})); %Find behavior data that matches this unique recording ID
    if length(idx) < 2; continue; end %If there are <2 units (thus no possible unit pairs), continue to next unique recording ID
    peth_save = cell(1,2); pethZ_save = cell(1,2); 
    prc50_save = cell(1,2); prc95_save = cell(1,2);
    fr_mvmt = []; fr_rest = [];
    h = waitbar(0,'Aligning spikes to events');
    for y = 1:length(idx)
        switch lbl
            case 'burst'
                try events = sub(idx(1)).burst.Windows(:,1); 
                catch continue; end %Extract event time points for this unit
            case 'pause'
                try events = sub(idx(1)).pause.AllSpikes; 
                catch continue; end 
            case 'spike'; events = sub(idx(1)).st;
        end
        st = {sub(idx(2:end)).st}; %Spike times to align from remaining units
        
        %% MOV vs REST
        evSub = cell(1,2); %Extract for 1st unit
        evSub{1} = extractEventST(events,beh(b).on/Fs,beh(b).off/Fs,0);
        evSub{2} = extractEventST(events,beh(b).onRest/Fs,beh(b).offRest/Fs,0);
        
        for n = 1:length(evSub)
            peth = getClusterPETH(st,evSub{n},gen); %PETH: aligning spikes to events from 1 unit
            peth_save{n} = [peth_save{n}, peth.fr];
            
            peth_z = []; tmp50 = []; tmp95 = [];
            for z = 1:length(st)
                st_shuff = shuffleST(st{z}, nShuff); 
                peth_shuff = getClusterPETH(st_shuff,evSub{n},gen); %PETH with shuffled spike times for a unit
                mu = nanmean(peth_shuff.fr(:)); sigma = nanstd(peth_shuff.fr(:)); 
                peth_z(:,z) = (peth.fr(:,z) - mu)./sigma;
                tmp50(:,z) = prctile(peth_shuff.fr,50,2); %50th percentile of shuffled PETH
                tmp95(:,z) = prctile(peth_shuff.fr,95,2); %95th percentile of shuffled PETH
            end
            pethZ_save{n} = [pethZ_save{n}, peth_z];
            prc50_save{n} = [prc50_save{n}, tmp50]; prc95_save{n} = [prc95_save{n}, tmp95];
        end
        
        for z = 1:length(st)
            fr_mvmt = [fr_mvmt; 1/mean(diff(extractEventST(st{z},beh(b).on/Fs,beh(b).off/Fs,1)))];
            fr_rest = [fr_rest; 1/mean(diff(extractEventST(st{z},beh(b).onRest/Fs,beh(b).offRest/Fs,1)))];
        end
        waitbar(y/length(idx),h); 
    end; close(h)
    %%
    mat(x).rec = uni{x};
    
    mat(x).fr_mvmt = fr_mvmt;
    mat(x).peth_mvmt = peth_save{1}; 
    mat(x).pethZ_mvmt = pethZ_save{1};
    mat(x).prc50_mvmt = prc50_save{1}; %50th percentile of shuffled PETH
    mat(x).prc95_mvmt = prc95_save{1}; %95th percentile of shuffled PETH
    
    mat(x).fr_rest = fr_rest;
    mat(x).peth_rest = peth_save{2}; 
    mat(x).pethZ_rest = pethZ_save{2};
    mat(x).prc50_rest = prc50_save{2}; %50th percentile of shuffled PETH
    mat(x).prc95_rest = prc95_save{2}; %95th percentile of shuffled PETH
end
time = peth.time;
fprintf('Done!\n');

%% (2) EXTRACT FROM MAT
delta = []; delta95 = []; delta50 = [];

for x = 1:length(mat)
    if isempty(mat(x).peth_rest); continue; end
    for y = 1:length(mat(x).fr_rest)
        tmp = (mat(x).peth_rest(:,y) - mat(x).fr_rest(y))./mat(x).fr_rest(y); %deltaFR = (rate - mu(rate))/mu(rate)
        delta = [delta, tmp];
        tmp_95 = (mat(x).prc95_rest(:,y) - mat(x).fr_rest(y))./mat(x).fr_rest(y); 
        tmp_50 = (mat(x).prc50_rest(:,y) - mat(x).fr_rest(y))./mat(x).fr_rest(y);
        delta95 = [delta95, tmp_95];
        delta50 = [delta50, tmp_50];
    end
end

%% (3) PLOT DELTA-FIRING RATE
figure; 
shadederrbar(time, nanmean(delta50,2), nanmean(delta95-delta50,2), 'k'); hold on
shadederrbar(time, nanmean(delta,2), SEM(delta,2), 'b'); 
xlabel('Lag (s)'); ylabel('PETH (deltaFR)'); xlim([-0.5 1]); grid on; 
title(sprintf('CIN spikes to rest %s (n = %d pairs)',lbl,size(delta,2)));
% 
% figure; violinplot(max(delta,[],1)); xlim([0.6 1.4]); ylim([-0.1 1.5]);
% ylabel('Max deltaFR'); grid on
% title(sprintf('Max deltaFR (n = %d pairs)',length(max(ccgDelta,[],1))));

%% Plotting average CIN to CIN burst alignment
figure;
%shadederrbar(peth.time, zeros(length(peth.time),1), mean(sem_all).*ones(length(peth.time),1),'k');
%shadederrbar(peth.time, nanmean(mat,2), SEM(mat,2), 'b');
shadederrbar(peth.time, nanmean(mat_b{2},2), SEM(mat_b{2},2), 'r');
shadederrbar(peth.time, nanmean(mat_b{1},2), SEM(mat_b{1},2), 'g');
xlabel(sprintf('Latency from %s',lbl{1})); ylabel('Firing Rate (z-score)'); xlim([-0.5 1])
grid on;
title(sprintf('CIN activity aligned to %s (n = %d pairs)',lbl{1},size(mat,2)));


%% Align PHOTOMETRY to bursts of CINs in same recordings
lbl = {'CIN burst'};
Fs = 50; time = [-1:1/Fs:1]; 

mat_a = cell(length(sub),1); mat_b = cell(1);
h = waitbar(0,'align photometry to CIN bursts');

for x = 1:length(sub)
    if lbl{1} == 'CIN burst'
        try events = sub(x).burst.Windows(:,1); catch continue; end %Extract event time points for this unit
    elseif lbl{1} == 'CIN pause'
        try events = sub(x).pause.Windows(:,1); catch continue; end 
    else; continue; 
    end
    b = find(strcmp({beh.rec},sub(x).rec));
    fp = beh(b).FP{1};
    
    evSub = {events};
%     evSub = cell(1,2);
%     for a = 1:length(beh(b).on)
%         logicalIndexes = events <= beh(b).off(a)/Fs & events >= beh(b).on(a)/Fs;
%         evSub{1} = [evSub{1}; events(logicalIndexes)];
%     end
%     for a = 1:length(beh(b).onRest)
%         logicalIndexes = events <= beh(b).offRest(a)/Fs & events >= beh(b).onRest(a)/Fs;
%         evSub{2} = [evSub{2}; events(logicalIndexes)];
%     end
    
    evShuff = shuffleST(events, 1); evShuff = evShuff{1}; %shuffle spike times n times
    mat = getSTA(fp, evShuff, Fs, [time(1), time(end)]);
    mu = nanmean(nanmean(mat,2)); sigma = nanmean(nanstd(mat,[],2)); %mu, sigma from shuffled events
    
    for y = 1:length(evSub)
        mat = getSTA(fp, evSub{y}, Fs, [time(1), time(end)]);
        for z = 1:size(mat,2)
            mat(:,z) = (mat(:,z) - mu)./sigma; %z-score using mu, sigma from shuffled events
            mat(:,z) = movmean(mat(:,z),5);
        end
        mat_a{x,y} = mat;
        mat_b{y} = [mat_b{y}, mat];
    end
    waitbar(x/length(sub),h);
end
close(h);

%% Plotting comparison
figure; clr = {'g','r'};
for x = 1:length(mat_b)
    shadederrbar(time, nanmean(mat_b{x},2), SEM(mat_b{x},2), clr{x}); hold on
end
xlabel(sprintf('Latency from %s',lbl{1})); ylabel(sprintf('%s (z-score)',beh(1).FPnames{1})); 
xlim([-0.5 1.0]); grid on; grid minor
title(sprintf('Photometry aligned to %s (n = %d)',lbl{1},length(sub)));

%% Plotting average FP to CIN burst alignment
figure; clr = {'g','r'}; plm = floor(sqrt(length(sub))); pln = ceil(length(sub)/plm);
for x = 1:length(sub)
    if isempty(mat_a{x,1}); continue; end
    sp(x) = subplot(plm,pln,x);
    for y = 1:size(mat_a,2)
        shadederrbar(time, nanmean(mat_a{x,y},2), SEM(mat_a{x,y},2), clr{y}); hold on; 
    end
	xlabel(sprintf('Latency from %s',lbl{1})); ylabel(sprintf('%s (z-score)',beh(1).FPnames{1})); 
    xlim([-0.5 1.0]); grid on; grid minor
	title(sprintf('%s-%d',sub(x).rec,sub(x).n));
end
linkaxes(sp,'y');
 
