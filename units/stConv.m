

%% Input Variables
params.total_duration = data.gen.totalTime; %in seconds
params.set_length = data.gen.setLength; %in seconds
params.Fs = data.gen.acqFs; %sampling frequency 
params.dsRate = params.Fs/data.gen.Fs; newFs = params.Fs/params.dsRate; %downsample
params.tau_rise = 30; %in milliseconds
params.tau_decay = 450; %in milliseconds

%% Iterate over all CINs
st = {data.clusters(data.lbl.cin).spikeTimes}; 

fprintf('\n st conv: tau_decay = %ds, dsRate = %d, nUnits = %d \n', ...
    params.tau_decay/1000, params.dsRate, length(st));

stSim = []; %clear/initialize matrix
for x = 1:length(st)
    spike_times = st{x}; %spike times, in seconds

    tic
    signal = st_conv(spike_times, params);
    toc

    signal = signal(1:params.set_length*newFs); %truncate to match length of analog signals
    stSim(:,x) = signal;
end
time = [1 : length(signal)]./data.gen.Fs; time = time(:); 
    
%% PLOT: spike times convolution + bouts
figure; hold on
bouts = zeros(length(data.final.time),1);
for ii = 1:data.final.mov.numBouts
    bouts(data.final.mov.onsets(ii):data.final.mov.offsets(ii)) = 1;
end
area(data.final.time,max(max(stSim)).*bouts,'FaceColor','g','FaceAlpha',.1,'EdgeAlpha',.5,'ShowBaseLine','off','DisplayName','bouts')
plot(data.final.time, data.final.vel, 'k', 'DisplayName', 'vel');

for x = 1:size(stSim,2)
    %subplot(size(stSim,2),1,x); 
    plot(time, stSim(:,x),'DisplayName',sprintf('#%d',data.lbl.cin(x)))
end


xlabel('Time (s)'); xlim([0 params.set_length]); legend
title(sprintf('%s-%s: convolved CIN spike times (tau = %ds)',data.mouse,data.date,params.tau_decay/1000));

%% Combined CIN spike times
stComb = [];
for x = 1:length(st)
    stComb = [stComb; st{x}]; 
end
tic
stCombSim = st_conv(stComb, total_duration, ...
    Fs, dsRate, tau_rise, tau_decay);
toc
stCombSim = stCombSim(1:set_length*newFs); 
%% PLOT: convolved combined CIN spike times + bouts
figure; hold on
area(data.final.time,max(stCombSim).*bouts,'FaceColor','g','FaceAlpha',.1,'EdgeAlpha',.5,'ShowBaseLine','off','DisplayName','bouts')
plot(time, stCombSim,'DisplayName','combST')
xlabel('Time (s)'); legend
title(sprintf('%s-%s: convolved CIN spike times (tau = %ds)',data.mouse,data.date,tau_decay/1000));

%% Plotting ROC curves
pred = stSim; 
resp = bouts > 0; 

figure; hold on
title(sprintf('ROC glm - CIN - %s-%s',data.mouse,data.date)); 
xlabel('False Positive Rate'); ylabel('True Positive Rate'); legend('Location','southeast')
plot([0:1],[0:1],'--k','DisplayName','chance'); %plot chance line

AUC = []; 
for x = 1:size(pred,2)
    mdl = fitglm(pred(:,x),resp,'Distribution','binomial','Link','logit'); %fit generalized logistic linear regression model to data
    scores = mdl.Fitted.Probability; %probability estimates from GLM as scores
    [X,Y,~,AUC(x)] = perfcurve(resp, scores, 'true'); %compute standard ROC curve using probabilities for scores
    plot(X,Y,'DisplayName',sprintf('#%d: %1.2f',data.lbl.cin(x),AUC(x))); 
end

%Plotting ROC curves for combined ST
mdl = fitglm(predComb,resp,'Distribution','binomial','Link','logit'); %fit generalized logistic linear regression model to data
scores = mdl.Fitted.Probability; %probability estimates from GLM as scores
[X,Y,~,AUCcomb] = perfcurve(resp, scores, 'true'); %compute standard ROC curve using probabilities for scores
plot(X,Y,'LineWidth',2,'DisplayName',sprintf('comb: %1.2f',AUCcomb));

fprintf('CIN convolved with exp kernel (tau = %ds)\n',params.tau_decay/1000);

%% 
tBefore = 4; tAfter = 4;
data.final.nUnits = length(data.lbl.cin);
for x = 1:data.final.nUnits
    data.final.stconv{x} = stSim(:,x);
end
data = alignStconv2Mov(data, tBefore, tAfter);
time = [-1*tBefore:0.02:tAfter]';
%%
figure; x = 1;
%plot(data.mov.onset.avgST{x}); 
shadederrbar(time,mean(data.mov.onset.allST{x},2),SEM(data.mov.onset.allST{x},2),'b');
