%% PLOT: STATS -- DA to ACh pause -- comparing peak magnitude, peak latency
group = [1*ones(nAn,1);2*ones(nAn,1);3*ones(nAn,1)];
[~,~,stats] = anova1(max_time(:),group,'off'); [c_time] = multcompare(stats,'Display','off');
[~,~,stats] = anova1(max_val(:),group,'off'); [c_val] = multcompare(stats,'Display','off');

fig = figure; fig.Position([3 4]) = [879 852];
subplot(2,2,1); hold on
r = [1 3]; a = max_time(:,r);
errorbar(nanmean(a)',SEM(a,1)','.m','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),a','.:k','MarkerSize',20);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'immobility','reward'});
ylabel('Latency to ACh pause (s)'); ylim([-0.4 0]); yticks([-0.4:0.1:0]);
[~, p] = ttest(a(:,1),a(:,2));
title(sprintf('peak latency (imm/rew: %1.3f)',p)); axis('square');

subplot(2,2,2); hold on
r = [1 3]; a = max_val(:,r);
errorbar(nanmean(a)',SEM(a,1)','.m','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),a','.:k','MarkerSize',20);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'immobility','reward'});
ylabel('rDA1m (%dF/F)'); ylim([-5 15]); yticks([0:5:15]);
[~, p] = ttest(a(:,1),a(:,2));
title(sprintf('peak amplitude (imm/rew: %1.3f)',p)); axis('square');

movegui(gcf,'center');

% PLOT: AVERAGES for all groups
clr = {'r','g','b'};
subplot(2,2,3);
plot([0 0],[-2 8],'--k'); hold on; 
for y = [1 3]
shadederrbar(sta_time, nanmean(da2pause_an{y},2), SEM(da2pause_an{y},2), clr{y}); hold on
end
xlabel('Latency to ACh pause (s)'); ylabel('rDA1m (%dF/F)'); xlim([-1 1]);
axis('square')
title(sprintf('DA to ACh pause (n = %d mice)',nAn));

subplot(2,2,4);
plot([0 0],[-6 3],'--k'); hold on; 
for y = [1 3]
shadederrbar(sta_time, nanmean(ach2pause_an{y},2), SEM(ach2pause_an{y},2), clr{y}); hold on
end
xlabel('Latency to ACh pause (s)'); ylabel('ACh3.0 (%dF/F)'); xlim([-1 1]);
axis('square')
title(sprintf('ACh to ACh pause (n = %d mice)',nAn));
movegui(gcf,'center');