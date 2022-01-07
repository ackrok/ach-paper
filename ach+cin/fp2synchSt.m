%% (1) Load data
% load('C:\Users\Anya\Desktop\IV_LOCAL\data_comb\fullcin+beh_ACh.mat')
sub = cinACh(find(strcmp({cinACh.rec},'IV066_rec01')));
beh = behACh(find(strcmp({behACh.rec},sub(1).rec)));

%% (2) PETH
%bin = 0.05; window = [-0.075 0.075]; %CHANGE: window for PETH
bin = 0.05; window = [-0.025 0.025];
%bin = 0.01; window = [-0.005 0.005];

mat = struct; % Initialize structure
h = waitbar(0, 'PETH: CIN spikes to FP peak times');
for x = 1:length(sub)
    %% Extract spike times
    idx = find(strcmp({sub.rec},sub(x).rec)); idx(idx == x) = [];
    if isempty(idx); continue; end
    st = [sub(x).st]; % Extract spike times of this unit
    st_other = {sub(idx).st}; % Extract spike times of all other units from this recording
    fr_other = [sub(idx).fr];
    %% PETH
    peth = getClusterPETH(st_other, st, bin, window); % PETH: spike times aligned to spike times
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
    % fp = [beh(1).vel(1); diff(movmean(beh(1).vel,10))]; Fs = 50;
    [sta_fp,sta_time] = getSTA(fp, st, Fs, [-1, 1]);
    %% Load into output structure
    mat(x).rec = sub(x).rec; 
    mat(x).n = sub(x).n; mat(x).m = [sub(idx).n];
    mat(x).cts0 = cts0; mat(x).cts = cts; mat(x).prob = prob; 
    mat(x).sta = sta_fp;
    %% REST vs MVMT
    idxSub = cell(1,2);
    for y = 1:length(beh(1).on)
        logicalIndexes = st < beh(1).off(y)/Fs & st > beh(1).on(y)/Fs;
        idxSub{1} = [idxSub{1}; find(logicalIndexes)];
    end
    for y = 1:length(beh(1).onRest)
        logicalIndexes = st < beh(1).offRest(y)/Fs & st > beh(1).onRest(y)/Fs;
        idxSub{2} = [idxSub{2}; find(logicalIndexes)];
    end
    mat(x).cts0_mvmt = mat(x).cts0(:,idxSub{1}); 
    mat(x).cts0_rest = mat(x).cts0(:,idxSub{2}); 
    mat(x).sta_mvmt = mat(x).sta(:,idxSub{1}); 
    mat(x).sta_rest = mat(x).sta(:,idxSub{2}); 
    %%
    waitbar(x/length(sub),h);
end
close(h); fprintf('Done: aligning spike times to spike times \n');
time = peth.time; 

%% PLOT ALL UNITS: STA to synchST
figure; plm = floor(sqrt(length(mat))); pln = ceil(length(mat)/plm);
a = cell(2, length(mat)); % a = cell(length(mat),length(mat)); 
clr = {'m','r','b','c','g'}; % clr = {'k','m','c','g'}; % %clr = {'k','r','m','b','c','g'}; %
% clr = {'k','g'};
for x = 1:length(mat)
    b = mat(x).cts0_rest; % CHANGE
    b(b > 1) = 1; b = sum(b,1);
    % a{1,x} = find(b == 0); a{2,x} = find(b == max(b)); 
     a{1,x} = find(b == 0); a{2,x} = find(b == 1); 
     a{3,x} = find(b == 2); a{4,x} = find(b == 3);
     a{5,x} = find(b == 4); %a{6,x} = find(b >= 5);
    sp(x) = subplot(plm,pln,x); 
    sig = mat(x).sta_rest; % CHANGE
    for y = 1:size(a,1)
        shadederrbar(sta_time, nanmean(sig(:,a{y,x}),2), SEM(sig(:,a{y,x}),2), clr{y}); hold on
    end
    xlabel('Latency to CIN spike (s)'); ylabel('Photometry (% dF/F)'); grid on
    title(sprintf('%s-#%d',mat(x).rec,mat(x).n),'Interpreter','none');
end
linkaxes(sp,'y');

%% FULL vs MOV vs REST
x = 5;
fig = figure; fig.Position(3) = 1800;
a = cell(length(mat),length(mat)); 
% clr = {'k','m','c','g'}; %clr = {'k','r','m','b','c','g'}; %
z = 1; sp(z) = subplot(1,3,z);
b = sum(mat(x).cts0,1);
a{1,z} = find(b == 0); a{2,z} = find(b == 1); 
a{3,z} = find(b == 2); a{4,z} = find(b == 3);
a{5,z} = find(b >= 4); %a{6,z} = find(b >= 5);
for y = 1:size(a,1)
    shadederrbar(sta_time, nanmean(mat(x).sta(:,a{y,z}),2), SEM(mat(x).sta(:,a{y,z}),2), clr{y}); hold on
end
xlabel('Latency to CIN spike (s)'); ylabel('ACh (dF/F)');
title(sprintf('%s-#%d',mat(x).rec,mat(x).n),'Interpreter','none');

z = 2; sp(z) = subplot(1,3,z);
b = sum(mat(x).cts0_mvmt,1);
a{1,z} = find(b == 0); a{2,z} = find(b == 1); 
a{3,z} = find(b == 2); a{4,z} = find(b == 3);
a{5,z} = find(b >= 4); %a{6,z} = find(b >= 5);
for y = 1:size(a,1)
    shadederrbar(sta_time, nanmean(mat(x).sta_mvmt(:,a{y,z}),2), SEM(mat(x).sta_mvmt(:,a{y,z}),2), clr{y}); hold on
end
xlabel('Latency to CIN spike (s)'); ylabel('ACh (dF/F)');
title(sprintf('%s-#%d - MOV',mat(x).rec,mat(x).n),'Interpreter','none');

z = 3; sp(z) = subplot(2,3,z);
b = sum(mat(x).cts0_rest,1);
a{1,z} = find(b == 0); a{2,z} = find(b == 1); 
a{3,z} = find(b == 2); a{4,z} = find(b == 3);
a{5,z} = find(b == 4); a{6,z} = find(b >= 5);
for y = 1:size(a,1)
    shadederrbar(sta_time, nanmean(mat(x).sta_rest(:,a{y,z}),2), SEM(mat(x).sta_rest(:,a{y,z}),2), clr{y}); hold on
end
xlabel('Latency to CIN spike (s)'); ylabel('ACh (dF/F)');
title(sprintf('%s-#%d - REST',mat(x).rec,mat(x).n),'Interpreter','none');

linkaxes(sp,'y');

%% PIE GRAPH
b = cellfun('length',a); % Return length of elements in cell array
b_prop = [];
for x = 1:length(mat); b_prop(:,x) = b(:,x)./size(mat(x).cts0,2); end
nanmean(b_prop.*100,2)

figure; bar(b_prop'.*100,'stacked');
