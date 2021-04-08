%% Compute ACG
st = {data.clusters(data.lbl.cin).spikeTimes}; %extract spike times from data structure
%st = {sub.st};

acgLin = []; %initialize matrixes
for x = 1:length(st)
    [xLin, nLin, ~, ~]= myACG(st{x},[],[]); %compute autocorrelogram with linear bins
    nLin = [flip(nLin);nLin]; %symmetrical ACG
    xLin = [-1.*flip(xLin);xLin]; xLin = xLin.*1000; %x-values for ACG in milliseconds
    acgLin = [acgLin, nLin]; %add to matrix
end

%% PLOT ACGs in subplots
figure;
plM = floor(sqrt(length(st))); plN = ceil(length(st)/plM); c = hsv(length(st)); 
for x = 1:length(st)
    subplot(plM,plN,x); hold on
    bar(xLin, acgLin(:,x),'FaceColor',c(x,:),'FaceAlpha',0.5);
    xlim([-250 250]); xlabel('Lag (ms)');
    title(sprintf('ACG unit#%d',data.lbl.cin(x))); 
end

%% Plot overlapping ACGs
figure; hold on
c = hsv(length(st)); 
for x = 1:length(st)
    plot(xLin, acgLin(:,x),'Color',c(x,:),'DisplayName',sprintf('unit#%d',data.lbl.cin(x)));
end
xlim([-250 250]); xlabel('Lag (ms)');
title(sprintf('%s-%s: ACG',data.mouse,data.date)); 

%% Plot all CCGs from mat structuer
mat = struct; %NAME BASED ON COMB MAT STRUCTURE

figure; 
plM = floor(sqrt(length(mat))); plN = ceil(length(mat)/plM); c = hsv(length(mat)); 
for x = 1:length(mat)
    [xLin, nLin, ~, ~]= myACG(mat(x).st,[],[]); %compute autocorrelogram
    nLin = [flip(nLin);nLin]; %symmetrical ACG
    xLin = [-1.*flip(xLin);xLin]; xLin = xLin.*1000; %x-values for ACG in milliseconds
    
    subplot(plM,plN,x); 
    bar(xLin, nLin, 'FaceColor',c(x,:), 'FaceAlpha',0.5);
    xlim([-250 250]); xlabel('Lag (ms)');
end

