%% Input Variables
beh = behwt(10:12); sub = cinwt(20:29); %DMS recordings
beh = behwt([1:9,13:28,36:40]); sub = cinwt([1:19,30:98]); %DLS recordings

%% Align acceleration to CIN unit spikes
Fs = 50; time = [-1:1/Fs:1]; %CHANGE: window for STA
align = {length(sub),2}; nShuff = 10; alignZ = align;
h = waitbar(0, 'STA: acceleration to CIN unit spike times');
for x = 1:length(sub)
    st = sub(x).st; % Aligning signal to spike times
    if isempty(sub(x).burst); continue; end
    st = sub(x).burst.Windows(:,1); % Aligning signal to burst onset times
    
    idx = find(strcmp({beh.rec},sub(x).rec));
    sig = [beh(idx).vel(1); diff(movmean(beh(idx).vel,10))]; % Acceleration
%     [mat,~,matZ] = getSTA(sig, st, Fs, [time(1), time(end)]);
    
    stSub = cell(1,2);
    for y = 1:length(beh(idx).on)
        logicalIndexes = st < beh(idx).off(y)/Fs & st > beh(idx).on(y)/Fs;
        stSub{1} = [stSub{1}; st(logicalIndexes)];
    end
    for y = 1:length(beh(idx).onRest)
        logicalIndexes = st < beh(idx).offRest(y)/Fs & st > beh(idx).onRest(y)/Fs;
        stSub{2} = [stSub{2}; st(logicalIndexes)];
    end
    for y = 1:2
        [mat,~,matZ] = getSTA(sig, stSub{y}, Fs, [time(1), time(end)]);
        align{x,y} = mat; alignZ{x,y} = matZ;
    end
    
%     ev_new = shiftST(st, nShuff, 1/nShuff); %shuffle spike times n times
%     mat_shuff = [];
%     for z = 1:nShuff
%         tmp_shuff = getSTA(sig, ev_new{z}, Fs, [time(1), time(end)]);
%         mat_shuff(:,z) = nanmean(tmp_shuff,2); 
%     end
%     mu = nanmean(nanmean(mat_shuff,2)); sigma = nanmean(nanstd(mat_shuff,[],2));
%     mat_z = (mat - nanmean(mu))./nanmean(sigma);
%     align{x} = mat; alignZ{x} = matZ;
    waitbar(x/length(sub),h);
end
close(h); fprintf('Done! \n');

%%
align_mat = []; alignZ_mat = [];
for x = 1:length(align); y = 1;
    if isempty(align{x}); align_mat(:,x) = nan(length(time),1); continue; end
    align_mat(:,x) = nanmean(align{x},2); 
    alignZ_mat(:,x) = nanmean(alignZ{x},2);
end
%
figure; hold on
shadederrbar(time, nanmean(alignZ_mat,2), SEM(alignZ_mat,2), 'k'); 
ylabel('Acceleration (z-score)'); xlabel('Latency (s)'); grid on; xlim([-1 1]);
title(sprintf('Acceleration aligned to CIN burst (n = %d units)',size(align_mat,2)));

%%
figure; plm = floor(sqrt(length(align))); pln = ceil(length(align)/plm);
for x = 1:length(align)
    if isempty(align{x}); continue; end
    sp(x) = subplot(plm,pln,x); 
    shadederrbar(time, nanmean(alignZ{x},2), SEM(alignZ{x},2), 'k'); 
    xlabel('Latency to CIN st (s)'); grid on; xlim([-1 1]);
    title(sprintf('%s-#%d',sub(x).rec,sub(x).n));
end
