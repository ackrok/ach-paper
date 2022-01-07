z = 1; % PAUSE(1) or PEAK(2)

switch z; case 1; lbl = 'Pause'; case 2; lbl = 'Peak'; end
%% PLOT
fig = figure; fig.Position(3) = 1375;

%IMM - frequency
freq_mat = [];
for y = 1:3; freq_mat(:,y) = s(y).freq(:,z); end
freq_mat(isnan(freq_mat)) = 0;
% group = [1*ones(4,1);2*ones(4,1);3*ones(4,1)];
% xx = freq_mat; [~,~,stats] = anova1(xx(:),group,'off'); [c] = multcompare(stats,'Display','off');
subplot(1,3,1);
violinplot(freq_mat);
xticklabels({s.inf}); xlim([0.5 3.5]); 
ylabel(sprintf('%s Frequency (Hz)',lbl)); ylim([-0.02 0.5]); yticks([0:0.1:0.5]);
axis('square');
title(sprintf('IMM: %s Freq (a/d: %1.3f | a/n: %1.3f)',lbl,signrank(freq_mat(:,1),freq_mat(:,2)),signrank(freq_mat(:,1),freq_mat(:,3))));

%IMM - duration
dur_mat = [];
for y = 1:3; for x = 1:4; dur_mat(x,y) = nanmean(s(y).dur{x,z}); end; end
dur_mat(isnan(dur_mat)) = 0;
% group = [1*ones(4,1);2*ones(4,1);3*ones(4,1)];
% xx = dur_mat; [~,~,stats] = anova1(xx(:),group,'off'); [c] = multcompare(stats,'Display','off');
subplot(1,3,2);
violinplot(dur_mat.*(1000/Fs));
xticklabels({s.inf}); xlim([0.5 3.5]); 
ylabel(sprintf('%s Duration (ms)',lbl)); ylim([-20 500]); yticks([0:100:500]);
axis('square');
title(sprintf('IMM: %s Dur (a/d: %1.3f | a/n: %1.3f)',lbl,signrank(dur_mat(:,1),dur_mat(:,2)),signrank(dur_mat(:,1),dur_mat(:,3))));

%IMM - amplitude
amp_mat = [];
for y = 1:3; for x = 1:4; amp_mat(x,y) = nanmean(s(y).mag{x,z}); end; end
amp_mat(isnan(amp_mat)) = 0;
switch z; case 1; amp_mat = -1*amp_mat; end
% group = [1*ones(4,1);2*ones(4,1);3*ones(4,1)];
% xx = amp_mat; [~,~,stats] = anova1(xx(:),group,'off'); [c] = multcompare(stats,'Display','off');
subplot(1,3,3);
violinplot(amp_mat);
xticklabels({s.inf}); xlim([0.5 3.5]); 
ylabel(sprintf('%s Amplitude (df/f)',lbl)); 
switch z; case 1; ylim([-0.2 6]); yticks([0:5]); case 2; ylim([-0.2 15]); yticks([0:5:15]); end
axis('square');
title(sprintf('IMM: %s Amp (a/d: %1.3f | a/n: %1.3f)',lbl,signrank(amp_mat(:,1),amp_mat(:,2)),signrank(amp_mat(:,1),amp_mat(:,3))));

movegui(gcf,'center');