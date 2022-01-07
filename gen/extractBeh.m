%%
fPath = 'R:\tritsn01labspace\'; 
[fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
if ~iscell(fName); fName = {fName}; end
%select .mat files you want to add to summary data structu
%%
if ~exist('beh','var'); beh = struct; end
for y = 1:length(fName) 
    load(fullfile(fPath,fName{y})); 
    [an,b] = strtok(fName{y},'_'); day = strtok(b,'_');
    x = 1+length(beh);
    beh(x).rec = [an,'_',day]; 
    beh(x).site = 'DLS';
    beh(x).task = 'wheel';
    
    %% Photometry
    beh(x).Fs = data.gen.Fs; % Sampling frequency, in Hz
    beh(x).time = data.final.time; % Time vector
    beh(x).FP = data.final.FP;
    beh(x).nbFP = data.final.nbFP; % Photometry signal(s)
    beh(x).FPnames = data.final.FPnames; % Names of photometry signal(s)
    
    %% Movement
    if isfield(data.final,'vel') % If movement data exists
        % data.final.vel = -1*data.final.vel;
        % data = processOnsetOffset(data,data.gen.params);
        % data = processRestOnsetOffset(data,data.gen.params);
        % save(fullfile(fPath,fName{y}),'data');
        beh(x).vel = data.final.vel; % Velocity signal
        beh(x).on = data.final.mov.onsets; beh(x).off = data.final.mov.offsets;                 % Movement onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
        beh(x).onRest = data.final.mov.onsetsRest; beh(x).offRest = data.final.mov.offsetsRest; % Rest onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
    else % Open field
        try 
            beh(x).cam = data.final.cam.on; 
        catch
            signal = data.acq.opto.trace;
            sigEdge = data.gen.params.FP.sigEdge; 
            rawFs = data.gen.params.acqFs; dsRate = data.gen.params.dsRate;
            if sigEdge ~= 0
                signal = signal((sigEdge*rawFs)+1:end-(sigEdge*rawFs));
            end
            camOn = getPulseOnsetOffset (signal, 0.5);
            camOn_Fs = round(camOn/dsRate);
            beh(x).cam = camOn_Fs;
        end
    end
    
    %% Lick/Reward
    if isfield(data.acq,'rew')
        % data = processReward(data, data.gen.params);
        beh(x).task = 'reward';
        % beh(x).cue = data.final.rew.cue;            % Cue onset times in sampling freq (data.gen.Fs), NOT in seconds
        beh(x).reward = data.final.rew.onset;    % Reward delivery time in sampling freq (data.gen.Fs), NOT in seconds
        beh(x).lick = data.final.lick.onset;        % Lick times in sampling freq (data.gen.Fs), NOT in seconds
        beh(x).lickVec = data.final.lick.trace;        % Lick trace
    end
    
    %%
    fprintf('Extracted from %s\n',fName{y});

%%
%     x = 1+length(raw);
%     raw(x).rec = [an,'-',day]; raw(x).Fs = data.gen.acqFs;
%     raw(x).FPnames{1} = 'GFP';
%     raw(x).rawfp = data.acq.FP{1};
%     raw(x).on{1} = data.final.mov.onsets;
%     raw(x).off{1} = data.final.mov.offsets;
%     raw(x).onRest{1} = data.final.mov.onsetsRest;
%     raw(x).offRest{1} = data.final.mov.offsetsRest;
end
if isempty(beh(1).Fs); beh(1) = []; end

%% Add recs from one structure to another
s1 = beh; s2 = modACh; % Names of structure. Adding s1 to s2

for x = 1:length(s1)
    y = 1+length(s2);
    f = fieldnames(s2);
    for z = 1:length(f)
        if isfield(s1,f{z})
            s2(y).(f{z}) = s1(x).(f{z});
        end
    end
end

%%
raw = struct;
for y = 1:length(fName)
    load(fullfile(fPath,fName{y})); 
    [an,b] = strtok(fName{y},'_'); day = strtok(b,'_'); task = 'wheel';
    x = 1+length(raw);
    raw(x).rec = [an,'-',day]; raw(x).Fs = 2000;
    raw(x).FPnames{1} = 'DA'; raw(x).acq = 'ACh';
    raw(x).rawfp = data.acq.FP{1};
end