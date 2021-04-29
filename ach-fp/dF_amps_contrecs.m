%% Proper Rest Threshold
velThres = 2.5;
Fs = 50;
minRestTime = 4*Fs; minRunTime = 1*Fs;
timeThres = 4*Fs; timeShift = 0.5*Fs;
for x = 1:length(beh)
    try
    beh(x).onRest = []; beh(x).offRest = [];
    vel = beh(x).vel; vel = abs(vel);
    [onsetInd,offsetInd] = getOnsetOffset(-vel,-velThres,minRunTime,minRestTime,1);
    [onsetInd,offsetInd] = adjOnsetOffset(onsetInd,offsetInd,timeThres,vel);
    onsetInd = onsetInd+timeShift; offsetInd = offsetInd-timeShift;
    beh(x).onRest = onsetInd; beh(x).offRest = offsetInd;
    end
end; fprintf('Done %s: new rest threshold = 0.25.\n',beh(1).FPnames{1});

%% Extract ACh dF/F values for REST vs MOV
mat = struct;
for x = 1:length(beh)
    mat(x).rec = beh(x).rec;
    mat(x).FPnames = beh(x).FPnames;
    fp = beh(x).FP{1}; % Extract photometry signal for this sweep
    % fp = raw(x).fp;
    mat(x).fp = fp;
    mat(x).fp_mvmt = []; mat(x).fp_rest = []; % Initialize output vectors for this recording
    t_mr = cell(2,2); t_mr{1,1} = beh(x).on; t_mr{2,1} = beh(x).off; % Load mov/rest onsets/offsets into cell structure
    t_mr{1,2} = beh(x).onRest; t_mr{2,2} = beh(x).offRest; % Load mov/rest onsets/offsets into cell structure
    tmp_sig = cell(1,2);
    for a = 1:2 % Repeat over movement, rest
        tmp_vec = []; % Clear/initialize vector
        for z = 1:length(t_mr{1,a}) % Iterate over all mov/rest periods
            range = [t_mr{1,a}(z), t_mr{2,a}(z)];  % Sample range
            tmp_vec = [tmp_vec; fp(range(1):range(2))]; % Extract values within range
        end
        tmp_sig{a} = tmp_vec;
    end
    mat(x).fp_mvmt = [mat(x).fp_mvmt; tmp_sig{1}]; % Concatenate to output matrix
    mat(x).fp_rest = [mat(x).fp_rest; tmp_sig{2}]; % Concatenate to output matrix
end
fprintf('Done %s: extracted REST + MOV photometry signal.\n',beh(1).FPnames{1});

fp_rest_all = []; fp_mvmt_all = []; % Initialize output vector for concatenated photometry dF/F
for x = 1:length(mat)
    fp_rest_all = [fp_rest_all; mat(x).fp_rest]; % Concatenate to output vector
    fp_mvmt_all = [fp_mvmt_all; mat(x).fp_mvmt]; % Concatenate to output vector
end
