Fs = beh(1).Fs; 
time = [-2:1/Fs:2]; %CHANGE: window for STA
align = cell(length(beh),2); nShuff = 10; 
h = waitbar(0, 'STA: signal to acc pks');
for x = 1:length(beh)
    for y = 1:length(beh(x).FP)
        sig = beh(x).FP{y}; %CHANGE signal
        peakProm = 0.1; peakDist = 0.5; % Parameters for FP peaks
        switch y
            case 1; fp = beh(x).FP{2}; case 2; fp = beh(x).FP{1}; 
        end
        fp_norm = (fp - min(fp)) / (max(fp) - min(fp)); % Min Max Normalization
        [pks,locs] = findpeaks(fp_norm,'MinPeakProminence',peakProm,'MinPeakDistance',peakDist); 
        locs = beh(x).time(locs);
        %%
        ev = locs; % Photometry peaks
        ev = extractEventST(locs, beh(x).on/Fs, beh(x).off/Fs, 1); % Event times during movement
        ev = extractEventST(locs, beh(x).onRest/Fs, beh(x).offRest/Fs, 1); % Event times during rest
        ev = extractEventST(locs, beh(x).reward/Fs, beh(x).reward/Fs + 0.4, 1);
        if isempty(ev); continue; end
        %%
        [mat,~,mat_z] = getSTA(sig, ev, Fs, [time(1), time(end)]); % STA: aligning photometry to event times
        align{x,y} = mat_z; % Load STA into output cell array
    end
    waitbar(x/length(beh),h);
end
close(h); fprintf('Done aligning photometry to events! \n');

%%
figure;
plm = floor(sqrt(size(align,1))); pln = ceil(size(align,1)/plm); % Subplot size depending on number of recordings
clr = {'g','r','b'}; 
lbl = 'Reward Delivery'; lbl = 'Movement Onset'; lbl = 'Acceleration Peak'; lbl = 'ACh peak';
for x = 1:length(align)
    sp(x) = subplot(plm,pln,x); y = 2;
    %x = 13; y = 1;
    %plot(time, align{x,y}, ':'); hold on % Plot average STA for each photometry signal
    shadederrbar(time, nanmean(align{x,y},2), SEM(align{x,y},2), clr{3}); hold on % Plot average STA for each photometry signal
    %shadederrbar(time, nanmean(align_mov{x,y},2), SEM(align_mov{x,y},2), clr{1}); hold on % Plot average STA for each photometry signal
    %shadederrbar(time, nanmean(align_rest{x,y},2), SEM(align_rest{x,y},2), clr{2}); hold on % Plot average STA for each photometry signal
    xlabel(sprintf('Latency to %s (s)',lbl)); 
    ylabel('DA (z-score)'); grid on; xlim([time(1) time(end)]);
    title(sprintf('%s - %s',beh(x).rec,beh(x).site)); 
    xlim([-1 1])
end; linkaxes(sp,'y'); 

%%
align_GZ3 = cell(1,3); y = 2;
align_GZ3{1} = [align_mov{7,y}, align_mov{8,y}, align_mov{9,y}];
align_GZ3{2} = [align_rest{7,y}, align_rest{8,y}, align_rest{9,y}];
align_GZ3{3} = [align_rew{2,y}, align_rew{5,y}, align_rew{8,y}, align_rew{11,y}, align_rew{14,y}];

% align_GZ2 = cell(1,3); y = 2;
% align_GZ2{1} = [align_mov{4,y}, align_mov{5,y}, align_mov{6,y}];
% align_GZ2{2} = [align_rest{4,y}, align_rest{5,y}, align_rest{6,y}];
% align_GZ2{3} = [align_rew{1,y}, align_rew{4,y}, align_rew{7,y}, align_rew{10,y}, align_rew{13,y}];

figure; clr = {'g','r','b'}; 
for x = 1:3
    shadederrbar(time, nanmean(align_GZ3{x},2), SEM(align_GZ3{x},2), clr{x}); hold on % Plot average STA for each photometry signal
end
% ylabel('DA (z-score)');
% yyaxis right; x = 3; shadederrbar(time, nanmean(align_GZ2{x},2), SEM(align_GZ2{x},2), clr{x}); 
lbl = 'ACh peak'; xlabel(sprintf('Latency to %s (s)',lbl)); 
ylabel('DA (z-score)'); grid on; xlim([time(1) time(end)]);
title('GZ003 - DLS DA to ACh peaks');
xlim([-1 1])
