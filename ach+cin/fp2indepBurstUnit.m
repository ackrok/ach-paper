%% (1) Load data
load('fullcin_Feb21_ACh+DA+PF.mat'); load('beh_wt_ACh+DA+PF.mat')
%sub = cinwt; beh = behwt;
sub = cinwt(find(strcmp({cinwt.rec},'IV066_rec01')));
beh = behwt(find(strcmp({behwt.rec},sub(1).rec)));

%% (2) PETH
gen = struct; 
gen.bin = 0.05; gen.window = [-0.025 0.025]; %CHANGE: window for PETH
%gen.bin = 0.02; gen.window = [-0.01 0.01];
%gen.bin = 0.01; gen.window = [-0.005 0.005];

mat = struct; % Initialize structure
h = waitbar(0, 'PETH: CIN spikes to FP peak times');
for x = 1:length(sub)
    %% Extract spike times
    idx = find(strcmp({sub.rec},sub(x).rec)); idx(idx == x) = [];
    if isempty(idx); continue; end
    st = sub(x).burst.Windows(:,1) + ((sub(x).burst.Windows(:,2) - sub(x).burst.Windows(:,1))/2); % Extract center of burst
    st_other = {sub(idx).st}; % Extract spike times of all other units from this recording
    fr_other = [sub(idx).fr];
    %% PETH
    peth = getClusterPETH(st_other, st, gen); % PETH: spike times aligned to spike times
    %%
    prob = []; 
    t_0 = 1; cts0 = []; cts = {};
    for y = 1:length(idx)
        cts{y} = peth.cts{y}; 
        prob(:,y) = (sum(peth.cts{y},2))./length(st_other{y}); % Probability of firing = cts/st
        cts0(y,:) = peth.cts{y}(t_0,:);
    end  
    %% STA
    fp = beh(1).FP{1}; Fs = 50;
    %fp = [beh(1).vel(1); diff(movmean(beh(1).vel,10))]; Fs = 50;
%     fp = beh(1).vel; Fs = 50;
    st_forSTA = sub(x).burst.Windows(:,1);
    [sta_fp,sta_time] = getSTA(fp, st_forSTA, Fs, [-1, 1]);
    %% Load into output structure
    mat(x).rec = sub(x).rec; 
    mat(x).n = sub(x).n; mat(x).m = [sub(idx).n];
    mat(x).cts0 = cts0; mat(x).cts = cts; mat(x).prob = prob; 
    mat(x).sta = sta_fp;
    %% REST vs MVMT
    st_mvmt = extractEventST(st, beh(1).on/Fs, beh(1).off/Fs, 0); % Event times during movement
    st_rest = extractEventST(st, beh(1).onRest/Fs, beh(1).offRest/Fs, 0); % Event times during movement
    peth_mvmt = getClusterPETH(st_other, st_mvmt, gen); % PETH: spike times aligned to spike times
    peth_rest = getClusterPETH(st_other, st_rest, gen); % PETH: spike times aligned to spike times
    cts_mvmt = {}; prob_mvmt = []; cts0_mvmt = [];
    cts_rest = {}; prob_rest = []; cts0_rest = [];
    for y = 1:length(idx)
        cts_mvmt{y} = peth_mvmt.cts{y};
        cts_rest{y} = peth_rest.cts{y};
        prob_mvmt(:,y) = (sum(peth_mvmt.cts{y},2))./length(st_other{y}); % Probability of firing = cts/st
        prob_rest(:,y) = (sum(peth_rest.cts{y},2))./length(st_other{y}); % Probability of firing = cts/st
        cts0_mvmt(y,:) = peth_mvmt.cts{y}(t_0,:);
        cts0_rest(y,:) = peth_rest.cts{y}(t_0,:);
    end  
    
    st_mvmtSTA = extractEventST(st_forSTA, beh(1).on/Fs, beh(1).off/Fs, 0); % Event times during movement
    st_restSTA = extractEventST(st_forSTA, beh(1).onRest/Fs, beh(1).offRest/Fs, 0); % Event times during movement
    sta_mvmt = getSTA(fp, st_mvmtSTA, Fs, [-1, 1]); 
    sta_rest = getSTA(fp, st_restSTA, Fs, [-1, 1]);
    
    mat(x).cts_mvmt = cts_mvmt; mat(x).cts_rest = cts_rest; 
    mat(x).prob_mvmt = prob_mvmt; mat(x).prob_rest = prob_rest;
    mat(x).cts0_mvmt = cts0_mvmt; mat(x).cts0_rest = cts0_rest;
    mat(x).sta_mvmt = sta_mvmt; mat(x).sta_rest = sta_rest;
    %%
    waitbar(x/length(sub),h);
end
close(h); fprintf('Done: aligning spike times to spike times \n');
time = peth.time; 

%% PLOT: STA to BURST UNITS
figure; 
a = cell(length(mat),length(mat)); 
clr = {'k','m','b','c','g'}; % clr = {'k','m','c','g'}; % %clr = {'k','r','m','b','c','g'}; %
for x = 1:length(mat)
    b = mat(x).cts0_rest; b(b > 1) = 1; b = sum(b,1);
    a{1,x} = find(b == 0); a{2,x} = find(b == 1); 
    a{3,x} = find(b == 2); a{4,x} = find(b >= 3);
%     a{5,x} = find(b >= 4); %a{6,x} = find(b >= 5);
    sp(x) = subplot(2,3,x); 
    for y = 1:size(a,1)
        shadederrbar(sta_time, nanmean(mat(x).sta_rest(:,a{y,x}),2), SEM(mat(x).sta_rest(:,a{y,x}),2), clr{y}); hold on
    end
    xlabel('Latency to CIN spike (s)'); ylabel('ACh (dF/F)'); 
    title(sprintf('%s-#%d',mat(x).rec,mat(x).n),'Interpreter','none');
end
linkaxes(sp,'y');

%% PLOT ALL UNITS: STA to synchST
figure; 
a = cell(length(mat),length(mat)); 
for x = 1:length(mat)
    sp(x) = subplot(2,3,x);
        
    mat = mat_IV066rec01_burst;
    b = mat(x).cts0; b(b > 1) = 1; b = sum(b,1);
    a{1,x} = find(b == 0); a{2,x} = find(b == 1); 
    a{3,x} = find(b == 2); a{4,x} = find(b >= 3);
    a{5,x} = find(b >= 4); %a{6,x} = find(b >= 5);
    y = 1; % Bursts where there is no overlap with other spikes
    shadederrbar(sta_time, nanmean(mat(x).sta(:,a{y,x}),2), SEM(mat(x).sta(:,a{y,x}),2), 'k'); hold on
    
    mat = mat_IV066rec01;
    b = mat(x).cts0; b(b > 1) = 1; b = sum(b,1);
    a{1,x} = find(b == 0); a{2,x} = find(b == 1); 
    a{3,x} = find(b == 2); a{4,x} = find(b >= 3);
    a{5,x} = find(b >= 4); %a{6,x} = find(b >= 5);
    y = 5; % Spikes where 4/4 units are spiking
    shadederrbar(sta_time, nanmean(mat(x).sta(:,a{y,x}),2), SEM(mat(x).sta(:,a{y,x}),2), 'g'); hold on
    
    xlabel('Latency (s)'); ylabel('ACh (dF/F)'); 
    title(sprintf('%s-#%d',mat(x).rec,mat(x).n),'Interpreter','none');
end
