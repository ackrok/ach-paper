z = 1; % PAUSE(1) or PEAK(2)

switch z; case 1; lbl = 'Pause'; case 2; lbl = 'Peak'; end
%% PLOT
fig = figure; fig.Position(3) = 1375; fig.Position(4) = 872;

%IMM - frequency
freq_mat = [];
for y = 1:3; freq_mat(:,y) = s(y).freq(:,z); end
freq_mat(isnan(freq_mat)) = 0;
% group = [1*ones(4,1);2*ones(4,1);3*ones(4,1)];
% xx = freq_mat; [~,~,stats] = anova1(xx(:),group,'off'); [c] = multcompare(stats,'Display','off');
subplot(2,3,1); hold on
r = [1 2]; a = freq_mat(:,r)./freq_mat(:,1);
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a'),'.:k','MarkerSize',20);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'aCSF','D1R/D2R antag'});
ylabel(sprintf('%s Frequency (Hz)',lbl)); ylim([0 3]); yticks([0:0.5:3]);
[~, p] = ttest(a(:,1),a(:,2));
title(sprintf('IMM: %s Freq (a/da: %1.3f)',lbl,p)); axis('square');

subplot(2,3,4); hold on
r = [1 3]; a = freq_mat(:,r)./freq_mat(:,r(1));
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a'),'.:k','MarkerSize',20);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'aCSF','NMDA/AMPA antag'});
ylabel(sprintf('%s Frequency (Hz)',lbl)); ylim([-0.5 1.5]); yticks([0:0.5:1]);
[~, p] = ttest(a(:,1),a(:,2));
title(sprintf('IMM: %s Freq (a/glu: %1.3f)',lbl,p)); axis('square');

%IMM - duration
dur_mat = [];
for y = 1:3; for x = 1:4; dur_mat(x,y) = nanmean(s(y).dur{x,z}); end; end
dur_mat(isnan(dur_mat)) = 0;
subplot(2,3,2); hold on
r = [1 2]; a = (dur_mat(:,r).*(1000/Fs))./(dur_mat(:,1).*(1000/Fs));
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a'),'.:k','MarkerSize',20);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'aCSF','D1R/D2R antag'});
ylabel(sprintf('%s Duration (ms)',lbl)); ylim([0 2]); yticks([0:0.5:2]); 
[~, p] = ttest(a(:,1),a(:,2));
title(sprintf('IMM: %s Dur (a/da: %1.3f)',lbl,p)); axis('square');

%IMM - amplitude
amp_mat = [];
for y = 1:3; for x = 1:4; amp_mat(x,y) = nanmean(s(y).mag{x,z}); end; end
amp_mat(isnan(amp_mat)) = 0;
switch z; case 1; amp_mat = -1*amp_mat; end
subplot(2,3,3); hold on
r = [1 2]; a = amp_mat(:,r)./amp_mat(:,1);
errorbar(abs(nanmean(a)'),abs(SEM(a,1)'),'.k','MarkerSize',20);
plot([1.15; 1.85].*ones(2,nAn),abs(a'),'.:k','MarkerSize',20);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'aCSF','D1R/D2R antag'});
ylabel(sprintf('%s Amplitude (df/f)',lbl)); 
switch z; case 1; ylim([0 2]); yticks([0:0.5:2]); case 2; ylim([0 3]); yticks([0:3]); end
[~, p] = ttest(a(:,1),a(:,2));
title(sprintf('IMM: %s Amp (a/da: %1.3f)',lbl,p)); axis('square');

movegui(gcf,'center');