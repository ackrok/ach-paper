function [out] = runglm_AChCIN(data)
% Generalized Linear Model (GLM) of CIN spiking and ACh photometry 
%   ACh = b0 + b1*CIN + error
%
% [out] = runglm_AChCIN(data)
% 
% Summary: Generalized linear model regression with response variable of
% ACh photometry and predictor variables of binned spike times
% Cross-validation by dividing time bins at random into N segments, then 
%   training GLM on N-1, and testing on the last segment. N repeats so each 
%   segment is tested once. 
%
%% Variables
    Fs = data.gen.Fs; %Sampling rate of photometry signal in samples/sec
    numSamp = data.gen.setLength*Fs; %Length of recording
    units = data.lbl.cin;
    
%% Response: ACh photometry
    fullY = data.final.FP{1}; %Processed photometry signal is response variable
    varY = data.final.FPnames;

%% Predictor: CIN spike times
    fullX = []; varX = {'b0'};
    for x = 1:length(units) 
        st_cin = data.clusters(units(x)).spikeTimes.*Fs; %Spike times, in samples matching that of photometry signal
        st_cin = st_cin(st_cin <= data.gen.setLength); %Retain only spikes that are within time range of photometry signal 
        bST_cin = hist(st_cin,numSamp); %Bin spike times
        fullX(:,x) = bST_cin(:); %Add to predictor matrix
        varX{x+1} = sprintf('unit #%d',units(x)); %Save unit number
    end
    
%% Response: Spike Times, Binned
    out = struct;
    x = 1; %for x = 1:length(units)
        
        %runX = fullX; 
        %runX(:,x) = []; %Remove column that matches current unit from predictor variables matrix

    %% Cross-Validation
    % Divide time bins at random into N segments. GLM is trained on N-1
    % segments and tested on the last segment. Repeat so each segment is
    % tested once. Allows for comparison of GLMs containing different numbers
    % of predictors because penalizes model overfitting. From Sjulson 2018
        temp = struct; 
        testIdx = randperm(numSamp); %Divide time bins at random into segments
        numRep = 10; %Number of segments 
        testIdx = reshape(testIdx, numSamp/numRep, numRep);
        h = waitbar(0,'Running GLM: ACh ~ CIN st');
        for rep = 1:numRep %Repeat 5x so each segment is tested once

            thisSegment = testIdx(:,rep); 
            testY = fullY(thisSegment,:);  %Set aside 1 segment (10%) for testing
            trainY = fullY; trainY(thisSegment,:) = []; %Remaining 9 segments (90%) will be used to train GLM
            testX = fullX(thisSegment,:);
            trainX = fullX; trainX(thisSegment,:) = [];

            [b,~,stats] = glmfit(trainX,trainY,'normal'); %Fit GLM to training set
            temp(rep).b = b; 
            temp(rep).p = stats.p;

            % Compute log-likelihood for poisson GLMs to compare performance. 
            % Let s be spike count in a bin and r is the predicted spike rate 
            % (known as "conditional intensity") in units of spikes/bin
            %   Poisson likelihood:         P(s|r) = r^s/s! exp(-r)
            %   giving log-likelihood:  log P(s|r) = s log r - r
            ratepred = exp(b(1) + testX*b(2:end)); %Evaluate GLM with testing set
            LL_glm = testY'*log(ratepred) - sum(ratepred);
            temp(rep).LL = LL_glm;

            %Calculate Root Mean Squared Error
            rmse = sqrt(mean(testY - ratepred).^2);
            temp(rep).rmse = rmse;
            
            waitbar(rep/numRep,h);
        end
        close (h);
        
    %% Log-likelihood
%     % Compute rate under homogenous Poisson model that asssumes a constant
%     % firing rate with correct mean spike count
%         ratepred_const = length(fullY)/numSamp; 
%         LL_0 = length(st).*log(ratepred_const) - numSamp.*ratepred_const;
% 
%     % Empirical single-spike information
%         SSinfo_glm = (mean([temp.LL]) - LL_0)/length(st)/log(2);
%         %fprintf('empirical single-spike information: GLM: %.2f bits/sp\n',SSinfo_glm);
% 
%     % P-values for coefficients
%         meanp = mean([temp.p],2);
%         %fprintf('GLM R(t) for unit #%d \n p-values of variables are: \n   b0: %.5f \n   b1 (velocity): %.5f \n',n,meanp(1),meanp(2));
%         if length(meanp) > 2
%             for y = 3:length(meanp)
%                 %fprintf('   b%d (CIN #%d): %.5f \n',(y-1),variables{y-2},meanp(y));
%             end
%         end

    %% Save into output structure
        out(x).mouse = data.mouse;
        out(x).date = data.date;
        out(x).fullY = fullY;
        out(x).fullX = fullX; 
        out(x).numRep = numRep;
        out(x).b = mean([temp.b],2);
        out(x).p = mean([temp.p],2);
        out(x).LL = mean([temp.LL],2);
        %out(x).LL_0 = LL_0;
        out(x).rmse = mean([temp.rmse],2);

    %end
    
end
