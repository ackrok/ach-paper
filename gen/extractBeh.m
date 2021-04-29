%%
fPath = 'R:\Grant Z\Fiber Photometry\'; 
[fName,fPath] = uigetfile('*.mat','MultiSelect','On');
%select .mat files you want to add to summary data structure

%%
% beh = struct;
for y = 1:length(fName)
    load(fullfile(fPath,fName{y})); 
    [an,b] = strtok(fName{y},'_'); day = strtok(b,'_'); task = 'openField';
    x = 1+length(beh);
    beh(x).rec = [an,'-',day]; beh(x).site = 'GFP-DLS';
    beh(x).task = task;
    
    %% Photometry
    beh(x).Fs = data.gen.Fs; % Sampling frequency, in Hz
    beh(x).time = data.final.time; % Time vector
    beh(x).FP = data.final.FP; % Photometry signal(s)
    beh(x).FPnames = data.final.FPnames; % Names of photometry signal(s)
    
    %% Movement
    if isfield(data.final,'vel') % If movement data exists
        beh(x).vel = data.final.vel; % Velocity signal
        beh(x).on = data.final.mov.onsets; beh(x).off = data.final.mov.offsets;                 % Movement onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
        beh(x).onRest = data.final.mov.onsetsRest; beh(x).offRest = data.final.mov.offsetsRest; % Rest onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
    else % Open field
        if isfield(data.acq,'opto'); signal = data.acq.opto.trace; end
        sigEdge = data.gen.params.FP.sigEdge; 
        rawFs = data.gen.params.acqFs; dsRate = data.gen.params.dsRate;
        if sigEdge ~= 0
            signal = signal((sigEdge*rawFs)+1:end-(sigEdge*rawFs));
        end
        camOn = getPulseOnsetOffset (signal, 0.5);
        camOn_Fs = round(camOn/dsRate);
        beh(x).cam = camOn_Fs;
    end
    
    %% Lick/Reward
    if isfield(data.acq,'rew')
        data = processReward(data, data.gen.params);
        beh(x).task = 'Pavlovian';
        beh(x).lick = data.final.lick.onset;        % Lick times in sampling freq (data.gen.Fs), NOT in seconds
        beh(x).cue = data.final.rew.cue;            % Cue onset times in sampling freq (data.gen.Fs), NOT in seconds
        beh(x).reward = data.final.rew.delivery;    % Reward delivery time in sampling freq (data.gen.Fs), NOT in seconds
    end
    
    %%
    fprintf('Extracted from %s\n',fName{y});
end

%%
%raw = struct;
for y = 1:length(fName)
    load(fullfile(fPath,fName{y})); 
    [an,b] = strtok(fName{y},'_'); day = strtok(b,'_'); task = 'wheel';
    x = 1+length(raw);
    raw(x).rec = [an,'-',day]; raw(x).Fs = 2000;
    raw(x).FPnames{1} = 'GFP';
    raw(x).rawfp = data.acq.FP{1};
end