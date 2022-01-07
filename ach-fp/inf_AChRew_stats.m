%% create structure with infusion photometry data
s(1).s = acsf; s(2).s = d1d2; s(3).s = glu;
s(1).inf = 'aCSF'; s(2).inf = 'D1R/D2R'; s(3).inf = 'NMDA/AMPA';
s(1).win = [20 40; 20 40; 20 40; 20 40]; s(2).win = s(1).win; 
s(3).win = s(1).win + 10; s(3).win(4,:) = [20 40];

[s(1).a, t, s(1).ev] = plot_fp2event(acsf, [-6 2]);
[s(2).a, ~, s(2).ev] = plot_fp2event(d1d2, [-6 2]);
[s(3).a, ~, s(3).ev] = plot_fp2event(glu, [-6 2]);

%% mu, std of rest ACh distribution during infusions
mu = []; sigma = []; % initialize/clear matrix
a_cell = {}; 
rew_min = {}; rew_lag = {}; rew_base = {};

z = 1; % ACh or DA
for y = 1:length(s) % iterate over infusion
    beh = s(y).s;
    for x = 1:length(beh) % iterate over animal
        ii = find(strcmp({cannula.rec},beh(x).rec)); ii = ii(1); % identify matching baseline recording from cannula larger structure
        fp_base = cannula(ii).FP{z}; % extract full baseline photometry signal
        fp_base_mu = nanmean(fp_base); % mean of baseline photometry recording
        fp_base = fp_base - fp_base_mu; % subtract mean of baseline recording to center around 0%
        fp_base_rest = fp_base(extractEventST([1:length(fp_base)]', cannula(ii).onRest, cannula(ii).offRest, 1)); % rest baseline
        
        a_rew = s(y).a{x,z}; % extract reward-aligned signal
        a_rew = a_rew - fp_base_mu; % subtract baseline%
        ev = s(y).ev{x}; % extract event times
        inf_win = s(y).win(x,:).*60;
        idx_inf = ev > inf_win(1) & ev < inf_win(2); % exclude events outside of infusion window
        a_rew = a_rew(:,idx_inf); % exclude events outside of infusion window
        a_rew_base = nanmean(a_rew(1:find(t==-1),:)); % baseline of pause
        [m,i] = min(a_rew(find(t == 0):find(t == 1.5),:));

        m = a_rew_base - m; % adjust for magnitude of pause
        i = i + find(t == 0) - 1; % adjust for start index
        
        a_cell{x,y} = a_rew;
        rew_min{x,y} = m(:); % reward pause magnitude (%dF/F)
        rew_lag{x,y} = t(i(:)); % reward pause latency (s)
        rew_base{x,y} = a_rew_base(:); % baseline (%dF/F)
    end
end

%%
% COLORS: teal - [0.05 0.75 0.45] || purple - [0.5 0.2 0.55] || orange - [0.85 0.35 0.1]
figure; hold on
y = 1; % sp(y) = subplot(1,2,1); hold on
plot([0 0],[-5 3],'--k');
avg_mat = []; for x = 1:4; avg_mat(:,x) = nanmean(a_cell{x,y},2); end
shadederrbar(t, nanmean(avg_mat,2), SEM(avg_mat,2), [0.05 0.75 0.45]); hold on
for x = 1:4; plot(t, nanmean(a_cell{x,y},2), 'Color', [0.05 0.75 0.45 0.2]); end
xlabel('Latency to Reward (s)'); ylabel('relative change ACh3.0 (%dF/F');
title(sprintf('ACh - %s',s(y).inf));
xlim([-1 2]);

y = 2; % sp(y) = subplot(1,2,2); hold on
plot([0 0],[-5 3],'--k');
avg_mat = []; for x = 1:4; avg_mat(:,x) = nanmean(a_cell{x,y},2); end
shadederrbar(t, nanmean(avg_mat,2), SEM(avg_mat,2), [0.5 0.2 0.55]); hold on
for x = 1:4; plot(t, nanmean(a_cell{x,y},2), 'Color', [0.5 0.2 0.55 0.2]); end
xlabel('Latency to Reward (s)'); ylabel('relative change ACh3.0 (%dF/F');
title(sprintf('ACh - %s',s(y).inf));
xlim([-1 2]);
linkaxes(sp,'y');

%%
y = [1 2]; % infusion
figure;
subplot(2,2,1); hold on; 
plot_me = rew_base; 
mu = cellfun(@nanmean, plot_me);
errorbar(nanmean(mu(:,y))',SEM(mu(:,y),1)','.k')
plot([mu(:,y)'],'.:k');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({s(y).inf});
ylabel('baseline (%dF/F)');
[~,p] = ttest(mu(:,y(1)),mu(:,y(2)));
title(sprintf('ACh - baseline (p = %1.3f)',p));

subplot(2,2,2); hold on
plot_me = rew_min; 
mu = cellfun(@nanmean, plot_me);
errorbar(nanmean(mu(:,y))',SEM(mu(:,y),1)','.k')
plot([mu(:,y)'],'.:k');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({s(y).inf}); 
ylabel('magnitude (%dF/F)');
[~,p] = ttest(mu(:,y(1)),mu(:,y(2)));
title(sprintf('ACh - pause magnitude (p = %1.3f)',p));

subplot(2,2,3); hold on
plot_me = rew_lag; 
mu = cellfun(@nanmean, plot_me);
errorbar(nanmean(mu(:,y))',SEM(mu(:,y),1)','.k')
plot([mu(:,y)'],'.:k');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({s(y).inf});
ylabel('latency (s)'); ylim([0 1]);
[~,p] = ttest(mu(:,y(1)),mu(:,y(2)));
title(sprintf('ACh - pause latency (p = %1.3f)',p));

% subplot(2,2,4); hold on
% plot_me = rew_lag; 
% mu = cellfun(@nanmean, plot_me);
% errorbar(nanmean(mu(:,y))',SEM(mu(:,y),1)','.k')
% plot([mu(:,y)'],'.:k');
% xlim([0.5 2.5]); xticks([1 2]); xticklabels({s(y).inf});
% ylabel('ACh pause duration (s)'); ylim([0 1]);
% title(sprintf('ACh - pause duration (p = %1.3f)',signrank(mu(:,y(1)),mu(:,y(2)))));

%% PLOT violinplot STD
figure;
violinplot(sigma); xticklabels({s.inf});
title('STD of ACh rest distribution (%dF/F)'); 
ylabel('standard deviation (%dF/F)'); ylim([0 2.5]);
 
%% PLOT histograms
figure;
bin = 0.2;
for x = 1:size(a_cell,1); for y = 1:size(a_cell,2)
    z = x + (4*(y-1)); sp(z) = subplot(3,4,z); hold on
    % histogram(fp_cell{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',[0.05 0.75 0.45],'FaceAlpha',0.2,'EdgeAlpha',0.1);
    % histogram(fp_2{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeAlpha',0.1);
    % histogram(fp_1{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',[0.05 0.75 0.45],'FaceAlpha',0.2,'EdgeAlpha',0.1);
    title(sprintf('%s - %s',strtok(s(y).s(x).rec,'_'),s(y).inf));
    end; end
linkaxes(sp,'x');

%% EXAMPLES: AK190
x = 2; % AK190
clr = [0.05 0.75 0.45; 0.5 0.2 0.55; 0.85 0.35 0.1]; % teal, purple, orange
for y = 1:3
    figure; hold on
    histogram(fp_2{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeAlpha',0.1);
    histogram(fp_1{x,y},'BinWidth',bin,'Normalization','probability','FaceColor',clr(y,:),'FaceAlpha',0.2,'EdgeAlpha',0.1); 
    legend({'pre-infusion',s(y).inf});
    title(sprintf('%s - %s',strtok(s(y).s(x).rec,'_'),s(y).inf));
    xlabel('ACh3.0 fluorescence (dF/F %)'); xlim([-5 25]);
end

