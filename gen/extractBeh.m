%%
fPath = 'R:\Grant Z\Fiber Photometry\'; 
[fName,fPath] = uigetfile('*.mat','MultiSelect','On');
%select .mat files you want to add to summary data structure

%%
beh = struct;
for x = 1:length(fName)
    load(fullfile(fPath,fName{x})); 
    [an,b] = strtok(fName{x},'_'); [site,b] = strtok(b,'_'); day = strtok(b,'_'); task = '';
    beh(x).rec = [an,'-',day]; 
    beh(x).site = site;
    beh(x).task = task;
    
    %% Photometry
    beh(x).Fs = data.gen.Fs; % Sampling frequency, in Hz
    beh(x).time = data.final.time; % Time vector
    beh(x).FP = data.final.FP; % Photometry signal(s)
    beh(x).FPnames = data.final.FPnames; % Names of photometry signal(s)
    
    %% Movement
    beh(x).vel = data.final.vel; % Velocity signal
    beh(x).on = data.final.mov.onsets; beh(x).off = data.final.mov.offsets;                 % Movement onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
    beh(x).onRest = data.final.mov.onsetsRest; beh(x).offRest = data.final.mov.offsetsRest; % Rest onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
    
    %% Lick/Reward
    if isfield(data.acq,'rew')
        data = processReward(data, data.gen.params);
        beh(x).task = 'Pavlovian';
        beh(x).lick = data.final.lick.onset;        % Lick times in sampling freq (data.gen.Fs), NOT in seconds
        beh(x).cue = data.final.rew.cue;            % Cue onset times in sampling freq (data.gen.Fs), NOT in seconds
        beh(x).reward = data.final.rew.delivery;    % Reward delivery time in sampling freq (data.gen.Fs), NOT in seconds
    end
    %%
    fprintf('Extracted from %s\n',fName{x});
end
cd(fPath)