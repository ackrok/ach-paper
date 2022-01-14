a = [7:12,14:16];
beh = allDMS(a); % Extract recordings with reward
[align_full, time, ev] = plot_fp2event(beh,[-6 2],0); % Align photometry to events
align_adj = align_full;
for x = 1:length(beh)
    align_adj{x,1} = align_adj{x,1} - nanmean(beh(x).FP{1}); end

% n = X mice
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
align_an = cell(length(uni),1); align_an_avg = [];
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    align_an_avg(:,x) = nanmean([align_adj{ii,1}],2);
    align_an{x} = [align_adj{ii,1}];
end

%% Population AVERAGE
figure;
shadederrbar(time, nanmean(align_an_avg,2), SEM(align_an_avg,2), 'g');
xlabel('Latency to reward (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('ACh3.0 fluorescence'); ylim([-3 1]); yticks([-2:0]);
title(sprintf('ACh3.0 DMS Population Average (n = %d mice)',size(align_an_avg,2)));

%% Population HEATMAP
figure;
h = heatmap(align_an_avg([find(time == -1):find(time == 2)],:)');
h.Title = 'ACh3.0 DMS Population Average';
h.XLabel = 'Latency to reward (s)'; h.YLabel = 'Animal';
h.Colormap = gray;

%% Example session INDIVIDUAL TRIALS
x = 2; % example AK193
z = sort(randi(size(align_an{x},2),[100,1]));
figure; hold on
plot(time, align_an{x}(:,z), 'Color', [0.05 0.75 0.45 0.5]);
shadederrbar(time, nanmean(align_an{x},2), SEM(align_an{x},2), 'k');
xlabel('Latency to reward (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('ACh3.0 fluorescence'); ylim([-5 7]); yticks([-4:2:10]);
title('ACh3.0 DMS Example (n = 100 trials)');

%% Licking
align_lick = cell(length(beh),1);
for x = 1:length(beh)
    lick = beh(x).lick/beh(x).Fs; % Licks, in seconds
    peth = getClusterPETH(lick, ev{x}(~isnan(ev{x})), 1/20, [-1 2]);
    align_lick{x} = peth.cts{1};
end
tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
lick_an = cell(length(uni),1); lick_an_avg = [];
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    lick_an_avg(:,x) = nanmean([align_lick{ii,1}],2);
    lick_an{x} = [align_lick{ii,1}];
end

%% LICK population AVERAGE
figure;
shadederrbar(peth.time, nanmean(lick_an_avg,2), SEM(lick_an_avg,2), 'k');
xlabel('Latency to reward (s)'); % xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('Lick (Hz)'); ylim([0 1]); yticks([0:0.2:1]);
title(sprintf('Lick Population Average (n = %d mice)',size(lick_an_avg,2)));

%% LICK population HEATMAP
figure;
h = heatmap(lick_an_avg');
h.Title = 'Lick Population Average';
h.XLabel = 'Latency to reward (s)'; 
h.YLabel = 'Animal';
h.Colormap = gray; caxis([0 1.1])

%% LICK example session RASTER
% use x from example above
zz = sort(randi(size(lick_an{x},2),[100,1]));
figure; hold on
plot([0 0],[0 100],'-k');
for ii = 1:length(zz)
    plot(peth.time, ii.*(lick_an{x}(:,zz(ii)) > 0), '.k');
end
xlabel('Latency to reward (s)'); xlim([-1 2]); xticks([-1:0.5:2]);
ylabel('Trial Number');
title('Lick Example Session (n = 100 trials)');