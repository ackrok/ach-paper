mat = struct;
beh = behwt(1:19); beh(10:12) = []; % ACh DLS recordings only
% x = 7; sub = cinwt(find(strcmp({cinwt.rec},beh(x).rec)));

for x = 1:length(beh)
    sub = cinwt(find(strcmp({cinwt.rec},beh(x).rec)));
    mat(x).rec = beh(x).rec;
    mat(x).nUnits = length(sub);
    if isempty(sub); continue; end
    tic
%%
st = {sub.st};
timeEnd = beh(x).time(end);
bin = 0.02; % Bin size, in seconds
timeBin = [0:bin:timeEnd]; timeVec = timeBin(1:end-1) + bin/2; % Bin edges, bin midpoints

stBin = zeros(length(timeBin)-1, length(st)); % Initialize binned spike times
stBinRand = zeros(length(timeBin)-1, length(st)); % Initialize 
for y = 1:length(st)
    stBin(:,y) = histcounts(st{y}, timeBin); % Bin spike times
    stRand = poissonSpikeGen(sub(y).fr, timeEnd, 1); % Generate random poisson spike train using firing rate of this unit
    % stRand = shuffleST(st{y}, 1);
    stBinRand(:,y) = histcounts(stRand{1}, timeBin); % Bin spike times
end
stBin (stBin > 1) = 1; % Binarize spiking in bin to 0's or 1's
stSynch = sum(stBin, 2); % Number of cells co-firing in each bin. Range: 0:length(sub)
stBinRand (stBinRand > 1) = 1; stSynchRand = sum(stBinRand, 2); 

%% Bin photometry
fp = beh(x).FP{1}; % Extract photometry
fpBin = downsampleTLab(fp, length(fp)/length(timeVec), 2); % Bin photometry
fpStrat = zeros(1+max(stSynch), 2); fpStrat(fpStrat == 0) = nan; % Initialize matrix
fpStratSEM = fpStrat;
fpAll = cell(1+max(stSynch), 2);
for y = 1:1+max(stSynch)
    idx = find(stSynch == y-1); % Find bin indices that satisfy x,y
    if isempty(idx); fpStrat(y,1) = 0; continue; end
    fpStrat(y,1) = mean(fpBin(idx)); % Mean dF/F
    fpStratSEM(y,1) = SEM(fpBin(idx),1); % SEM dF/F
    
    idxRand = find(sum(stSynchRand, 2) == y-1); % Find bin indices that satisfy x,y
    if isempty(idxRand); fpStrat(y,2) = 0; continue; end
    fpStrat(y,2) = mean(fpBin(idxRand)); % Mean dF/F
    fpStratSEM(y,2) = SEM(fpBin(idxRand),1); % SEM dF/F
    
    fpAll{y,1} = fpBin(idx); fpAll{y,2} = fpBin(idxRand);
end
norm = fpStrat(2,1); % Average ACh when 1 unit is active
fpStratNorm = fpStrat./norm; fpStratSEM = fpStratSEM./norm; % Normalize FP dF/F to ACh when #co-active units = 1

%% Save into output structure
    mat(x).fpStratNorm = fpStratNorm;
    mat(x).fpStratSEM = fpStratSEM;
    mat(x).fpStrat = fpStrat;
toc
end
%% (3) PLOT
figure; hold on
% h = heatmap([1:2],[0:max(stSynch)],[fpStrat,fpRand],'ColorMap',parula,'CellLabelColor','none');
% h = heatmap([1:2],[0:max(stSynch)],[fpStratNorm,fpRandNorm],'ColorMap',parula,'CellLabelColor','none');
% ylabel('# co-active units');
% plot([0:1+max(stSynch)],fpStratNorm(:,2),'-ok'); plot([0:4],fpStratNorm(:,1),'-ob');
shadederrbar([0:max(stSynch)], fpStratNorm(:,2), fpStratSEM(:,2), 'k'); hold on
shadederrbar([0:max(stSynch)], fpStratNorm(:,1), fpStratSEM(:,1), 'g');
xlabel('# co-active units'); ylabel('ACh (normalized dF/F)');
title(sprintf('%s - ACh',beh(x).rec));

%%
mat2 = mat; mat([mat.nUnits] < 2) = [];
figure;
for x = 1:length(mat)
    sp(x) = subplot(3,3,x); hold on
%     plot([0:mat(x).nUnits],mat(x).fpRandNorm,'o-k'); plot([0:mat(x).nUnits],mat(x).fpStratNorm,'o-b');
    shadederrbar([0:mat(x).nUnits], mat(x).fpStratNorm(:,2), mat(x).fpStratSEM(:,2), 'k'); hold on
    shadederrbar([0:mat(x).nUnits], mat(x).fpStratNorm(:,1), mat(x).fpStratSEM(:,1), 'g');
    xlabel('# co-active units'); ylabel('ACh (normalized dF/F)');
    title(sprintf('%s - ACh',mat(x).rec));
end; % linkaxes(sp,'y');

%% adjust x-axis to be 0-100% units
fpStratNorm_all = []; fpStratNorm_rand = [];
for x = 1:length(mat)
    fpStratNorm_all(x,1) = mat(x).fpStratNorm(1,1);
    fpStratNorm_all(x,2) = mat(x).fpStratNorm(2,1);
    fpStratNorm_all(x,3) = mat(x).fpStratNorm(end,1);
    
    fpStratNorm_rand(x,1) = mat(x).fpStratNorm(1,2);
    fpStratNorm_rand(x,2) = mat(x).fpStratNorm(2,2);
    fpStratNorm_rand(x,3) = mat(x).fpStratNorm(end,2);
end
figure;
shadederrbar([0:2], nanmean(fpStratNorm_rand,1), SEM(fpStratNorm_rand,1), 'k'); hold on
shadederrbar([0:2], nanmean(fpStratNorm_all,1), SEM(fpStratNorm_all,1), 'g');


