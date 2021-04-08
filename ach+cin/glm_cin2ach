%% Extract Variables
units = data.lbl.cin;
st = {data.clusters(units).spikeTimes}; 

dt = 0.001; %Bins of 1ms width - spike times become discrete probablity distribution
edges = 0:dt:data.gen.setLength; %Vector of time bin edges (for histogram)

tau_decay = 270; %Set the decay constant, in ms
tau_rise = 50; %Set the rise time, in ms
kernel = getExpKer(tau_rise,tau_decay,1/dt); %Get the kernel
kernel = kernel(1:100:end); %Downsample

tic
stConv = zeros(length(edges)-1,length(st)); 
for n = 1:length(units)
    temp = histcounts(st{n}, edges); %Bin spike times into 1ms bins = discrete probability distribution
    temp( temp > 1 ) = 1; %Correction for if have multiple spikes per bin (should not occur unless have ISI violations)
    fr = conv(temp, kernel ,'same'); %Smooth binned spikeTimes using kernel
    stConv (:,n) = fr; %Load into matrix
end
toc 

% rawFP = openDat;
% lpCut = 10; filtOrder = 8;
% rawFs = 30000; dsRate = 30; dsType = 2;
% interpType = 'linear'; fitType = 'interp'; basePrc = 5; winSize = 10; winOv = 0;
% FPbase = baselineFP(rawFP,interpType,fitType,basePrc,winSize,winOv,rawFs);   
% FPfilt = filterFP(FPbase,rawFs,lpCut,filtOrder,'lowpass');
% FPfinal = downsampleTLab(FPfilt,dsRate,dsType); Fs = rawFs/dsRate;
% FPset = FPfinal(1:length(stConv));
fp = FPset; %open /Users/akrok/Desktop/ACh Paper/glm/IV043rec02_fp.mat

%% Input Variables: GLM Parameters
options = optimoptions('fminunc',...
    'Algorithm','trust-region',...    
    'OptimalityTolerance',1e-6,...     
    'StepTolerance',1e-6,...                      
    'SpecifyObjectiveGradient', true);
    options.MaxIterations = 100;
    options.Display = 'off'; %'iter-detailed';
   
repeats = 10; %CHANGE - how many times to iterate GLM
numSamp = size(stConv,1); %length(respY);
testIdx = reshape(1:numSamp, numSamp/repeats, repeats); %Cross validation segments

ratepred_const = numSamp/numSamp; 
LL_0 = numSamp.*log(ratepred_const) - numSamp.*ratepred_const;

%% Run GLM
tic
out = struct; 

respY = fp; % ResponseY is photometry signal
predX = stConv; %PredictorX are convolved spike times

fprintf('\n Running GLM ...');
for rep = 1:repeats % Repeat so each segment is tested once
    thisSegment = testIdx(:,rep); 
    testY = respY(thisSegment,:);  % Set aside 1 segment for testing
    trainY = respY; trainY(thisSegment,:) = []; % Remaining segments will be used to train GLM
    testX = predX(thisSegment,:);
    testX = [ones(size(testX,1), 1) testX]; % Pad predictor matrix with 1's, this is your constant
    trainX = predX; trainX(thisSegment,:) = [];
    trainX = [ones(size(trainX,1), 1) trainX];

    b = rand(size(trainX,2), 1); % Generate random beta values
    costfun_train = @(x) nll_glm(x,trainX,trainY,dt); % Generate cost-function based on training data set, using negative log-likelihood computation
    costfun_test = @(x) nll_glm(x,testX,testY,dt); % Same, for testing data set

    [b_opt] = fminunc(costfun_train, b, options); % **Optimization of beta's using cost-function
    toc; fprintf('Optimization complete for repeat #%d of %d. ',rep,repeats);

    f_test = costfun_test(b_opt); % Return negative log-likelihood for optimized beta's
    
    out(rep).rep = rep;
    out(rep).b_opt = b_opt; % Save optimized beta's
    out(rep).LL = -f_test; % Save log-likelihood
    out(rep).essi = mean((-f_test) - LL_0)/length(respY)/log(2); %Empirical single-spike information 
end

fprintf('done! \n');

%% PLOT predX, respY
bOpt_avg = mean([out.b_opt],2); % Average optimized beta's
fp_glm = exp( stConv*bOpt_avg(2:end) ) .* dt; % Evaluate GLM model

dsRate = 10; dsType = 2; fsNew = (1/dt)/dsRate; 
t_vec = [1/fsNew:1/fsNew:data.gen.setLength];
fp_2 = downsampleTLab(respY,dsRate,dsType); % Downsample data
fp_glm2 = downsampleTLab(fp_glm,dsRate,dsType); % Downsample predicted fp
stConv_2 = []; for x = 1:length(st); stConv_2(:,x) = downsampleTLab(stConv(:,x),dsRate,dsType); end % Downsample convolved spike times

figure; hold on
yyaxis left; plot(t_vec,fp_2,'g'); ylabel('ACh (dF/F)'); xlabel('Time (s)');
yyaxis right; 
plot(t_vec,fp_glm2,'r');
plot(t_vec,stConv_2,'k');
