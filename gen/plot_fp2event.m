function [fig,align,time] = plot_fp2event(beh,varargin)
%Align photometry to acceleration peaks and plot in subplots
%
%   plot_fp2event(beh)
%   [fig,align,time] = plot_fp2event(beh)
%   [fig,align,time] = plot_fp2event(beh,window)
%
%   Description: This function is for aligning photometry signal to
%   discrete events (acceleration peaks) for multiple recordings using
%   input structure 'beh' for multiple photometry + behavior recordings.
%   The output of this function is figure with nSubplots = nRecordings and 
%   an optional cell array with the data used to generate each subplot.
%
%   Input:
%   - beh - A data structure containing extracts photometry, behavior
%   information for multiple recordings
%   - optional inputs:
%       - window - A vector with window for extraction, default is [-1:1]
%
%   Output:
%   - fig - Figure handle
%   - align - A cell array with photometry aligned to acceleration peaks
%   for each recording analyzed
%   - time - A vector of time points for plotting purposes
%
%   Author: Anya Krok, February 2020

%% Input Variables
evName = menu('Choose Event to Align To:','Movement Onset','Rest Onset','Movement Offset','Rest Offset','Acceleration Peaks','Cue Onset','Reward Onset','Lick');

%% Extract Acceleration Peaks
if evName == 5
    for x = 1:length(beh)
        if all(logical(~rem(beh(x).on,1))); diffFs = 1; else; diffFs = beh(x).Fs; end % Adjust for movement onset times being in samples or seconds
        vel = beh(x).vel; % Extract velocity signal during this mov bout
        vel_sm = fliplr(movmean(fliplr(movmean(vel,10)),10)); % Smooth velocity, flip left-right, smooth again, flip back
        acc = [vel(1); diff(vel_sm)]; % Acceleration vector is diff of smoothed velocity vector
        [~,locs] = findpeaks(acc,'MinPeakProminence',0.5,'MinPeakDistance',0.5); % Location of peaks, using findpeaks function
        locs = beh(x).time(locs); % Convert peak locations to seconds
        beh(x).acc_locs = locs; clc
    end; clc
end

%% Align photometry to acceleration times
Fs = beh(1).Fs; 
time = [-5:1/Fs:5]; %CHANGE: window for STA
align = cell(length(beh),length(beh(1).FP)); nShuff = 10; 
h = waitbar(0, 'STA: signal to acc pks');
for x = 1:length(beh)
    if all(logical(~rem(beh(x).on,1))); diffFs = 1; else; diffFs = 50; end
    switch evName
        case 1; ev = beh(x).on/(Fs/diffFs); % Extract movement onset times, adjusting event times to be in seconds if necessary
        case 2; ev = beh(x).onRest/(Fs/diffFs); % Extract rest onset times, adjusting event times to be in seconds if necessary
        case 3; ev = beh(x).off/(Fs/diffFs); % Extract movement offset times, adjusting event times to be in seconds if necessary
        case 4; ev = beh(x).offRest/(Fs/diffFs); % Extract rest offset times, adjusting event times to be in seconds if necessary
        case 5; ev = beh(x).acc_locs; % Extract acceleration peak times, adjusting event times to be in seconds if necessary
        case 6; ev = beh(x).cue/(Fs/diffFs); % Extract cue times, adjusting event times to be in seconds if necessary
        case 7; ev = beh(x).reward/(Fs/diffFs); % Extract reward delivery times, adjusting event times to be in seconds if necessary
        case 8; ev = beh(x).lick/(Fs/diffFs); % Extract lick times, adjusting event times to be in seconds if necessary
    end
    if isempty(beh(x).FP); continue; end %if no photometry, continue to next recording
    for y = 1:length(beh(x).FP)
        sig = beh(x).FP{y}; %CHANGE signal
        [mat,~,mat_z] = getSTA(sig, ev, Fs, [time(1), time(end)]); % STA: aligning photometry to event times
%         ev_new = shiftST(ev, nShuff, 1/nShuff); % Shuffle event times n times
%         mat_shuff = [];
%         for z = 1:nShuff
%             tmp_shuff = getSTA(sig, ev_new{z}, Fs, [time(1), time(end)]); % STA: aligning photometry to shuffled event times
%             mat_shuff = [mat_shuff, tmp_shuff]; 
%         end
%         mu = nanmean(nanmean(mat_shuff,2)); sigma = nanmean(nanstd(mat_shuff,[],2)); % Use mu, sigma of shuffled STA to z-score 
        align{x,y} = mat_z; % Load STA into output cell array
    end
    waitbar(x/length(beh),h);
end
close(h); fprintf('Done aligning photometry to events! \n');

%% Plot for each recording
fig = figure; % Figure handle
plm = floor(sqrt(size(align,1))); pln = ceil(size(align,1)/plm); % Subplot size depending on number of recordings
clr = {'g','m','b','r'};
switch evName
    case 1; lbl = 'Movement Onset'; case 2; lbl = 'Rest Onset'; % Asign label based on event time aligning to
    case 3; lbl = 'Movement Offset'; case 4; lbl = 'Rest Offset';
    case 5; lbl = 'Acceleration Peak';
    case 6; lbl = 'Cue Onset'; case 7; lbl = 'Reward Delivery'; case 8; lbl = 'Lick';
end
for x = 1:size(align,1) % Iterate over each recording
    sp(x) = subplot(plm,pln,x); 
    for y = 1:size(align,2) % Iterate over each photometry signal
        if isempty(align{x,y}); continue; end
        shadederrbar(time, nanmean(align{x,y},2), SEM(align{x,y},2), clr{y}); hold on % Plot average STA for each photometry signal
    end
    xlabel(sprintf('Latency to %s (s)',lbl)); 
    ylabel('FP (z-score)'); grid on; xlim([time(1) time(end)]);
    title(sprintf('%s - %s',beh(x).rec,beh(x).site)); 
end
%linkaxes(sp,'y'); 
linkaxes(sp,'x'); % Link x,y axes

end
