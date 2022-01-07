%% Align photometry data to reward
good_rew = [10:16,20,22:29,32:35,40:42,44:47];
beh = modAChDA(good_rew);
[align_full, t] = plot_fp2event(beh, [-6 2], 0);

%% n = X mice, concatenate all trials for each animal
an = {}; for x = 1:length(beh); an{x} = strtok(beh(x).rec,'-'); end
uni = unique(an);

align_adj = {}; % adjust for baseline
for x = 1:length(beh); for y = 1:2; align_adj{x,y} = align_full{x,y} - nanmean(beh(x).FP{y}); end; end

align_an = {}; % Cell array with concatenated trials for each animal
for x = 1:length(uni) % Iterate over each unique animal
    for y = 1:2 % Repeat for ACh, DA
        align_an{x,y} = [];
        matchAN = find(strcmp(an,uni{x})); % Find all recordings that match animal ID
        for z = 1:length(matchAN)
            align_an{x,y} = [align_an{x,y}, align_adj{matchAN(z),y}]; % Concatenate trials from recordings with same animal ID
        end
    end
end

%% n = X mice, average all trials for each animals
align_an_avg = cell(1,2); 
align_an_avg_norm = cell(1,2);
for x = 1:size(align_an,1)
    for y = 1:size(align_an,2)
        tmp = align_an{x,y};
        align_an_avg{1,y}(:,x) = nanmean(tmp,2); % Average over all trials for each animal
        tmp = nanmean(tmp,2);
        bb = nanmean(tmp([find(t==-6):find(t==-1)],:),1);
        
        if y == 1 
            mm = min(tmp([find(t==0):find(t==1)],:));
            norm = (tmp - mm)./(bb - mm);
        elseif y == 2
            mm = max(tmp([find(t==0):find(t==1)],:));
            norm = (tmp - bb)./(mm - bb); end
        align_an_avg_norm{1,y}(:,x) = norm; % normalize(align_an_avg{1,y}(:,x),'range'); % Normalize so range is [0 1]
    end
end
 
%% PLOT: %dF/F 
fig = figure;
plot([0 0],[-5 11],'--k'); hold on
shadederrbar(t, nanmean(align_an_avg{2},2), SEM(align_an_avg{2},2), 'm'); hold on
shadederrbar(t, nanmean(align_an_avg{1},2), SEM(align_an_avg{1},2), 'g');
ylabel('FP (% dF/F)'); xlabel('Latency to reward delivery (s)');
xlim([-1 2]);
axis('square')
title(sprintf('reward: ACh3.0 + rDA1m (n = %d mice)',size(align_an,1)));
 
%PLOT: %dF/F normalized range [0 1]
fig = figure;
plot([0 0],[-1.1 1.1],'--k'); hold on
shadederrbar(t, nanmean(align_an_avg_norm{2},2), SEM(align_an_avg_norm{2},2), 'm'); hold on
shadederrbar(t, nanmean(align_an_avg_norm{1}-1,2), SEM(align_an_avg_norm{1},2), 'g');
ylabel('FP (a.u.)'); xlabel('Latency to reward delivery (s)');
xlim([-0.2 0.8]); xticks([0 0.5]);
ylim([-1.1 1.1]); yticks([-1:1]);
axis('square')
title(sprintf('reward: ACh3.0 + rDA1m (n = %d mice)',size(align_an,1))); 

%% PLOT: ACh, DA separately %dF/F for each animal
figure;
plot([0 0],[0 17],'--k'); hold on
plot(t, align_an_avg{1}); xlim([-1 2]);
ylabel('ACh3.0 (% dF/F)'); xlabel('Latency to reward delivery (s)');
title(sprintf('reward: ACh3.0 (n = %d mice)',size(align_an,1)));
 
figure; 
plot([0 0],[-2 14],'--k'); hold on
plot(t, align_an_avg{2}); xlim([-1 2]);
ylabel('rDA1m (% dF/F)'); xlabel('Latency to reward delivery (s)');
title(sprintf('reward: rDA1m (n = %d mice)',size(align_an,1)));
 
%% n = X mice, find ACh pauses and DA peaks
save_fit = []; r2 = [];
% figure;
for x = 1:size(align_an,1) % Iterate over each animal
    y = 1; ach_min = min(align_an{x,y}([301:340],:)); % ACh pause: minimum value within defined time range: +0.0 +0.8 seconds after reward delivery
    y = 2; da_max = max(align_an{x,y}([301:340],:)); % DA peak: maximal value within defined time range: +0.0 +0.8 seconds after reward delivery
    x_vec = da_max(~isnan(da_max)); % Remove NaNs from vector of maxima
    y_vec = ach_min(~isnan(ach_min)); % Remove NaNs from vector of minma
    lm = fitlm(x_vec, y_vec); % Simple linear regression
    r2(x) = lm.Rsquared.Ordinary; % R-squared coefficient of determination
    da_range = [5:0.5:30]; % Defined DA peak range to fit line to
    ach_fit = polyval([lm.Coefficients{2,1}, lm.Coefficients{1,1}], da_range);  % Evaluate fitted line over defined range
    save_fit(:,x) = ach_fit; % Store fit values in matrix
 
%     sp(x) = subplot(4,3,x); hold on
%     plot(x_vec, y_vec, '.b'); % Plot points: ACh pause, DA peak for non-nan trials
%     plot([0:0.5:30], ach_fit, '-r'); % Plot fit line
%     title(sprintf('%s: %1.2f',uni{x},r2(x)));
%     xlabel('DA peak'); ylabel('ACh pause');
end
% linkaxes(sp,'x'); linkaxes(sp,'y');
 
%% PLOT: ACh pause ~ b0 + b1*DA peak?
figure; hold on;
plot(da_range, save_fit, '--');
shadederrbar(da_range, nanmean(save_fit,2), SEM(save_fit,2), 'k');
xlabel('DA peak values for fit (% dF/F)'); xlim([0 35]);
ylabel('Linear regression fit for ACh pause (% dF/F)'); ylim([-5 2]);
title(sprintf('Regression: ACh pause ~ DA peak (n = %d mice)',size(save_fit,2)));
figure;
violinplot(r2); title('Regression: R-squared'); 
ylim([-0.2 1]); yticks([0:0.5:1]); ylabel('R-squared');
xlim([0.2 1.8]); xticklabels({''}); grid on

