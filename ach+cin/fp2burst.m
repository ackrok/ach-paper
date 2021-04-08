%% Align FR CINs to bursts of CINs in same recordings
lbl = {'CIN pause'};
uni = unique({sub.rec});
peth_gen = struct; peth_gen.bin = 0.01; peth_gen.window = [-0.5 1];
mat_a = cell(length(uni),2); mat_b = cell(1,2);
%mat = []; mu_all = []; sem_all = []; nShuff = 5;

h = waitbar(0, 'PETH: CIN activity to CIN bursts (pairs)');
for x = 1:length(uni)
	idx = find(strcmp({sub.rec},uni{x})); %Find all units within same recording
    if length(idx) < 2; continue; end %Exit loop if don't have enough units to make pairs
    %b = find(strcmp({beh.rec},uni{x})); %Find matching behavior data

    for y = 1:length(idx) %Repeat for each unit
        if lbl{1} == 'CIN burst'
            try events = sub(idx(1)).burst.Windows(:,1); catch continue; end %Extract event time points for this unit
        elseif lbl{1} == 'CIN pause'
            try events = sub(idx(1)).pause.AllSpikes; catch continue; end 
        else; continue; 
        end
        st = {sub(idx(2:end)).st}; 
        %% MOV vs REST
        evSub = cell(1,2); %Extract for 1st unit
        b = find(strcmp({beh.rec},sub(x).rec));
        for a = 1:length(beh(b).on) %Extract events that occur within movement periods
            logicalIndexes = events <= beh(b).off(a) & events >= beh(b).on(a);
            evSub{1} = [evSub{1}; events(logicalIndexes)];
        end
        for a = 1:length(beh(b).onRest) %Extract events that occur during rest periods
            logicalIndexes = events <= beh(b).offRest(a) & events >= beh(b).onRest(a);
            evSub{2} = [evSub{2}; events(logicalIndexes)];
        end
        for n = 1:length(evSub)
            peth = getClusterPETH(st,evSub{n},peth_gen); %PETH: aligning activity of units to bursts of 1 unit
            peth_z = [];
            for z = 1:length(st)
                st_shuff = shuffleST(st{z}, nShuff); 
                peth_shuff = getClusterPETH(st_shuff,evSub{n},peth_gen); %PETH with shuffled spike times for a unit
                mu = nanmean(nanmean(peth_shuff.fr,2)); sigma = nanmean(nanstd(peth_shuff.fr,[],2)); 
                peth_z(:,z) = (peth.fr(:,z) - mu)./sigma;
                peth_z(:,z) = movmean(peth_z(:,z),5);
            end
            mat_a{x,n} = [mat_a{x,n}, peth_z];
            mat_b{n} = [mat_b{n}, peth_z];
        end
        %% SHUFFLE
%         peth = getClusterPETH(st,events,peth_gen);
%         peth_z = []; 
%         for z = 1:length(st)
%             st_shuff = shuffleST(st{z}, nShuff); 
%             peth_shuff = getClusterPETH(st_shuff,events,peth_gen); %PETH with shuffled spike times for a unit
%             mu = nanmean(nanmean(peth_shuff.fr,2)); sigma = nanmean(nanstd(peth_shuff.fr,[],2)); 
%             peth_z(:,z) = (peth.fr(:,z)-mu)./sigma; %z-score using mean, std of PETH with shuffled spike times
%             peth_z(:,z) = movmean(peth_z(:,z),5);
%             mu_all = [mu_all; mu]; sem_all = [sem_all; sigma/sqrt(nShuff)]; 
%         end
%         mat = [mat, peth_z]; %Save into matrix
        %%
        idx = circshift(idx, 1); %Shift index forward to repeat, using next unit as reference
    end
    waitbar(x/length(uni),h);
end
close(h); fprintf('Done! \n');

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
    sp(x) = subplot(plm,pln,x);
    for y = 1:size(mat_a,2)
        shadederrbar(time, nanmean(mat_a{x,y},2), SEM(mat_a{x,y},2), clr{y}); hold on; 
    end
	xlabel(sprintf('Latency from %s',lbl{1})); ylabel(sprintf('%s (z-score)',beh(1).FPnames{1})); 
    xlim([-0.5 1.0]); grid on; grid minor
	title(sprintf('%s-%d',sub(x).rec,sub(x).n));
end
linkaxes(sp,'y');
 
