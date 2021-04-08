%% (1) Load data
load('fullcin_Feb21_ACh+DA+PF.mat'); load('beh_wt_ACh+DA+PF.mat')
sub = cinwt; beh = behwt([5:9, 13:19]);
%Select ACh recordings

%% (2) BIN data
bin = 0.02; % bin size = 50ms

mat = struct;
h = waitbar(0, 'binning data');
for x = 1:length(beh)
    st = {sub(strcmp({sub.rec},beh(x).rec)).st};
    if length(st) < 2; continue; end
    
    %% (2a) BIN spike times
    timeEnd = beh(x).time(end);
    timeBin = [0:bin:timeEnd]; timeVec = timeBin(1:end-1) + bin/2; % Bin edges, bin midpoints
    stBin = zeros(length(timeBin)-1, length(st)); % Initialize binned spike times
    for y = 1:length(st)
        stBin(:,y) = histcounts(st{y}, timeBin); % Bin spike times
    end
    stBin (stBin > 1) = 1; % Binarize spiking in bin to 0's or 1's
    stSynch = sum(stBin, 2); % Number of cells co-firing in each bin. Range: 0:length(sub)
    clc
    %% (2b) Stratify stSynch
    synchLow = timeVec(stSynch < 1); % 0 units spiking in bin
    synchHigh = timeVec(stSynch == max(stSynch)); % ALL units spiking in bin
    
    %% (2c) STA stratified stSynch
    fp = beh(x).FP{1}; Fs = 50; % Extract photometry
    [staLow, staTime, staLowZ] = getSTA(fp, synchLow, Fs, [-1, 1]);
    [staHigh, ~, staHighZ] = getSTA(fp, synchHigh, Fs, [-1, 1]);
    
    %%
    mat(x).rec = beh(x).rec;
    mat(x).synchSt = stSynch;
    mat(x).staLow = staLow; mat(x).staLowZ = staLowZ; 
    mat(x).staHigh = staHigh; mat(x).staHighZ = staHighZ;
    mat(x).check = x;
    waitbar(x/length(beh),h);
end
close(h);
mat = mat([mat.check]); mat = rmfield(mat,'check');
        
%% (3) PLOT
figure; 
for x = 1:length(mat)
    subplot(3,3,x);
    shadederrbar(staTime, nanmean(mat(x).staLowZ,2), SEM(mat(x).staLowZ,2), 'k'); hold on
    shadederrbar(staTime, nanmean(mat(x).staHighZ,2), SEM(mat(x).staHighZ,2), 'g');
    xlabel('Latency to st (s)'); ylabel('ACh (z-score)');
    title(sprintf('%s - ACh to low/high synchSt',mat(x).rec));
end
