%% (1) Load data
load('fullcin_Feb21_ACh+DA+PF.mat'); load('beh_wt_ACh+DA+PF.mat')
sub = cinwt; beh = behwt([5:9, 13:19]);
%Select ACh recordings

%% (2) BIN data
bin = 0.02; % bin size = 40ms or 20ms (data needs to be divisible)

mat = struct;
h = waitbar(0, 'binning data');
for x = 1:length(beh)
    st = {sub(find(strcmp({sub.rec},beh(x).rec))).st};
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
    %% (2b) BIN velocity
    vel = beh(x).vel; % Extract velocity
    acc = [vel(x); diff(movmean(vel,10))]; % Extract acceleration
    behBin = downsampleTLab(vel, length(vel)/length(timeVec), 2); % Bin acceleration
    beh_max = round(ceil(max(behBin)), -1); % Round maximum
    behStrat = behBin;
    behStrat( behBin < 0.25 ) = 0; behStrat( behBin >= 0.25 & behBin < 10) = 1;
    for y = 2:beh_max/10
        behStrat( behBin >= (y-1)*10 & behBin < y*10 ) = y; % Stratify into low to high velocity
    end
    behStrat ( behBin >= y*10 ) = y+1;
    clc
    %% (2c) BIN photometry
    fp = beh(x).FP{1}; % Extract photometry
    fpBin = downsampleTLab(fp, length(fp)/length(timeVec), 2); % Bin photometry
    fpStrat = zeros(1+max(stSynch), 1+max(behStrat)); fpStrat(fpStrat == 0) = nan; % Initialize matrix
    for y = 1:1+max(stSynch)
        for z = 1:1+max(behStrat)
            idx = find(stSynch == y-1 & behStrat == z-1); % Find bin indices that satisfy x,y
            if isempty(idx); fpStrat(y,z) = 0; continue; end
            fpStrat(y,z) = mean(fpBin(idx));
        end
    end
    clc
    %%
    mat(x).rec = beh(x).rec;
    mat(x).synchSt = stSynch;
    mat(x).behStrat = behStrat;
    mat(x).fpStrat = fpStrat;
    waitbar(x/length(beh),h);
end
close(h);
        
%% (3) PLOT
figure; 
for x = 1:length(mat)
    subplot(3,3,x);
    h = heatmap([0:max(mat(x).behStrat)],[0:max(mat(x).synchSt)],mat(x).fpStrat,'ColorMap',parula,'CellLabelColor','none');
    xlabel('Velocity Binned (Low - High)'); ylabel('#Units co-active');
    title(sprintf('%s - ACh',mat(x).rec));
end
